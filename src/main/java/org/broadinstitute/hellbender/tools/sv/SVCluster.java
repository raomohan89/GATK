package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.codecs.SVCallRecordCodec;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.PrintStream;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Clusters SVs with similar breakpoints based on coordinates and supporting evidence.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         Unclustered structural variant VCF
 *     </li>
 *     <li>
 *         PE evidence file
 *     </li>
 *     <li>
 *         SR evidence file
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Structural variant VCF
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVCluster
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Clusters structural variants",
        oneLineSummary = "Clusters structural variants",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public final class SVCluster extends GATKTool {
    public static final String SPLIT_READ_LONG_NAME = "split-reads-file";
    public static final String DISCORDANT_PAIRS_LONG_NAME = "discordant-pairs-file";
    public static final String MIN_SIZE_LONG_NAME = "min-size";

    @Argument(
            doc = "Split reads file",
            fullName = SPLIT_READ_LONG_NAME
    )
    private String splitReadsFile;

    @Argument(
            doc = "Discordant pairs file",
            fullName = DISCORDANT_PAIRS_LONG_NAME
    )
    private String discordantPairsFile;

    @Argument(
            doc = "Input file",
            fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME
    )
    private String inputFile;

    @Argument(
            doc = "Output file",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputFile;

    @Argument(
            doc = "Min event size",
            fullName = MIN_SIZE_LONG_NAME,
            minValue = 0,
            maxValue = Integer.MAX_VALUE,
            optional = true
    )
    private int minEventSize = 50;
    private SAMSequenceDictionary dictionary;

    private List<Tuple2<SimpleInterval,List<SVCallRecord>>> currentClusters = new LinkedList<>();
    private String currentContig = null;
    private int numProcessedClusters = 0;
    private FeatureDataSource<SVCallRecord> reader;
    private PrintStream writer;
    private SVCallRecordCodec callRecordCodec;

    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;

    private final double MIN_RECIPROCAL_OVERLAP_DEPTH = 0.8;
    private final double MIN_RECIPROCAL_OVERLAP_PESR = 0.1;
    private final int BREAKEND_CLUSTERING_WINDOW = 300;

    private final int SPLIT_READ_WINDOW = 25;
    private final int DISCORDANT_PAIR_WINDOW = 500;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        callRecordCodec = new SVCallRecordCodec();
        loadSplitReadEvidenceDataSource();
        loadDiscordantPairDataSource();
        reader = new FeatureDataSource<>(inputFile, "inputFile", 100000, SVCallRecord.class, getDefaultCloudPrefetchBufferSize(), getDefaultCloudIndexPrefetchBufferSize());
        writer = new PrintStream(BucketUtils.createFile(outputFile));
    }

    @Override
    public Object onTraversalSuccess() {
        reader.close();
        writer.close();
        return null;
    }

    private void loadSplitReadEvidenceDataSource() {
        splitReadSource = new FeatureDataSource<>(
                splitReadsFile,
                "splitReadsFile",
                100000,
                SplitReadEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void loadDiscordantPairDataSource() {
        discordantPairSource = new FeatureDataSource<>(
                discordantPairsFile,
                "discordantPairsFile",
                100000,
                DiscordantPairEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private List<SplitReadEvidence> getStartSplitReads(final Collection<SVCallRecord> cluster) {
        return getClusterEvidence(cluster, SVCallRecord::getStartAsInterval, splitReadSource, SPLIT_READ_WINDOW);
    }

    private List<SplitReadEvidence> getEndSplitReads(final Collection<SVCallRecord> cluster) {
        return getClusterEvidence(cluster, SVCallRecord::getEndAsInterval, splitReadSource, SPLIT_READ_WINDOW);
    }

    private List<DiscordantPairEvidence> getStartDiscordantPairs(final Collection<SVCallRecord> cluster) {
        return getClusterEvidence(cluster, SVCallRecord::getStartAsInterval, discordantPairSource, DISCORDANT_PAIR_WINDOW);
    }

    private List<DiscordantPairEvidence> getEndDiscordantPairs(final Collection<SVCallRecord> cluster) {
        return getClusterEvidence(cluster, SVCallRecord::getEndAsInterval, discordantPairSource, DISCORDANT_PAIR_WINDOW);
    }

    private <T extends Feature> List<T> getClusterEvidence(final Collection<SVCallRecord> cluster,
                                                                final Function<SVCallRecord,SimpleInterval> intervalFunction,
                                                                final FeatureDataSource<T> dataSource,
                                                                final int padding) {
        final List<SimpleInterval> intervalList = cluster.stream()
                .map(r -> intervalFunction.apply(r))
                .collect(Collectors.toList());
        final SimpleInterval interval = IntervalUtils.getSpanningInterval(intervalList).expandWithinContig(padding, dictionary);
        return Lists.newArrayList(dataSource.query(interval));
    }

    private boolean clusterTogether(final SVCallRecord a, final SVCallRecord b) {
        // TODO does not handle BNDs / intrachromosomal events
        if (!a.getContig().equals(b.getContig())) return false;
        final SimpleInterval intervalA = new SimpleInterval(a.getContig(), a.getStart(), a.getEnd());
        final SimpleInterval intervalB = new SimpleInterval(b.getContig(), b.getStart(), b.getEnd());
        //TODO min or max here?
        final double reciprocalOverlap = Math.max(getReciprocalOverlap(a), getReciprocalOverlap(b));
        if (IntervalUtils.isReciprocalOverlap(intervalA, intervalB, reciprocalOverlap)) return true;
        final int window = Math.min(getBreakendWindow(a), getBreakendWindow(b));
        return Math.abs(intervalA.getStart() - intervalB.getStart()) <= window &&
                Math.abs(intervalA.getEnd() - intervalB.getEnd()) <= window;
    }

    private SimpleInterval getClusteringStartInterval(final SVCallRecord variant, final SimpleInterval clusterMinStartInterval) {
        if (clusterMinStartInterval != null && !clusterMinStartInterval.getContig().equals(variant.getContig())) {
            throw new IllegalArgumentException("Attempted to get clustering start interval for variant on " + variant.getContig() + " but the given cluster start interval is " + clusterMinStartInterval);
        }
        if (!currentContig.equals(variant.getContig())) {
            throw new IllegalArgumentException("Attempted to get clustering start interval for variant on " + variant.getContig() + " but the current contig is " + currentClusters);
        }
        if (dictionary.getSequence(variant.getContig()) == null) {
            throw new IllegalArgumentException("Variant contig " + variant.getContig() + " not found in sequence dictionary.");
        }
        final int start = variant.getStart();
        final double min_reciprocal_overlap = getReciprocalOverlap(variant);
        final int minStart = Math.max(start - BREAKEND_CLUSTERING_WINDOW, (int) (start - min_reciprocal_overlap * (variant.getEnd() - start)));
        final int maxStart = Math.max(start + BREAKEND_CLUSTERING_WINDOW, (int) (start + min_reciprocal_overlap * (variant.getEnd() - start)));
        if (clusterMinStartInterval == null) {
            return IntervalUtils.trimIntervalToContig(currentContig, minStart, maxStart, dictionary.getSequence(currentContig).getSequenceLength());
        }
        final int newMinStart = Math.min(minStart, clusterMinStartInterval.getStart());
        final int newMaxStart = Math.min(maxStart, clusterMinStartInterval.getEnd());
        return IntervalUtils.trimIntervalToContig(currentContig, newMinStart, newMaxStart, dictionary.getSequence(currentContig).getSequenceLength());
    }

    private double getReciprocalOverlap(final SVCallRecord variant) {
        return isDepthOnlyCall(variant) ? MIN_RECIPROCAL_OVERLAP_DEPTH : MIN_RECIPROCAL_OVERLAP_PESR;
    }

    private int getBreakendWindow(final SVCallRecord variant) {
        return isDepthOnlyCall(variant) ? 0 : BREAKEND_CLUSTERING_WINDOW;
    }

    private boolean isDepthOnlyCall(final SVCallRecord variant) {
        final List<String> algorithms = variant.getAlgorithms();
        for (final String alg : algorithms) {
            if (!alg.equals("depth")) return true;
        }
        return false;
    }

    private void processBuffer() {
        while (!currentClusters.isEmpty()) {
            processCluster(0);
        }
    }

    private void seedCluster(final SVCallRecord seed) {
        if (!currentContig.equals(seed.getContig())) {
            throw new IllegalArgumentException("Attempted to seed new cluster with variant on contig " + seed.getContig() + " but the current contig is " + currentContig);
        }
        final List<SVCallRecord> newCluster = new ArrayList<>(1);
        newCluster.add(seed);
        currentClusters.add(new Tuple2<>(getClusteringStartInterval(seed, null), newCluster));
    }

    private void addToCluster(final int clusterIndex, final SVCallRecord variant) {
        if (!currentContig.equals(variant.getContig())) {
            throw new IllegalArgumentException("Attempted to add new variant on contig " + variant.getContig() + " but the current contig is " + currentContig);
        }
        if (clusterIndex >= currentClusters.size()) {
            throw new IllegalArgumentException("Specified cluster index " + clusterIndex + " is greater than the largest index.");
        }
        Tuple2<SimpleInterval, List<SVCallRecord>> cluster = currentClusters.get(clusterIndex);
        final SimpleInterval clusterInterval = cluster._1;
        final List<SVCallRecord> clusterVariants = cluster._2;
        clusterVariants.add(variant);
        final SimpleInterval clusteringStartInterval = getClusteringStartInterval(variant, clusterInterval);
        if (clusteringStartInterval.getStart() != clusterInterval.getStart() || clusteringStartInterval.getEnd() != clusterInterval.getEnd()) {
            currentClusters.remove(clusterIndex);
            currentClusters.add(clusterIndex, new Tuple2<>(clusteringStartInterval, clusterVariants));
        }
    }

    private void processCluster(final int i) {
        if (i < 0 || i >= currentClusters.size()) {
            throw new IllegalArgumentException("Specified cluster index " + i + " is out of range.");
        }
        final Tuple2<SimpleInterval, List<SVCallRecord>> cluster = currentClusters.remove(i);
        final List<SVCallRecord> clusterVariants = cluster._2;

        final List<SplitReadEvidence> startSplitReads = getStartSplitReads(clusterVariants);
        final List<DiscordantPairEvidence> startDiscordantPairs = getStartDiscordantPairs(clusterVariants);

        final List<SplitReadEvidence> endSplitReads = getEndSplitReads(clusterVariants);
        final List<DiscordantPairEvidence> endDiscordantPairs = getEndDiscordantPairs(clusterVariants);


        //writeCluster(clusterVariants, Integer.toString(numProcessedClusters));

        numProcessedClusters++;
        progressMeter.update(cluster._1);
    }

    private List<SplitReadEvidence> getSplitReadSupport(final List<SVCallRecord> calls, final List<SplitReadEvidence> evidenceList) {
        final Set<String> carrierSamples = calls.stream().map(SVCallRecord::getSample).collect(Collectors.toSet());
        final int finalMaxPosition = getPositionWithMostSplitReadSupport(evidenceList, carrierSamples);
        return evidenceList.stream().filter(e -> e.getStart() == finalMaxPosition).collect(Collectors.toList());
    }

    private int getPositionWithMostSplitReadSupport(final List<SplitReadEvidence> evidenceList, final Set<String> carrierSamples) {
        int position = 0;
        int carrierCount = 0;
        int backgroundCount = 1;
        int maxPosition = 0;
        double maxRatio = 0;
        for (final SplitReadEvidence e : evidenceList) {
            if (e.getStart() != position) {
                final double ratio = carrierCount / (double) backgroundCount;
                if (ratio > maxRatio) {
                    maxPosition = position;
                    maxRatio = ratio;
                }
                position = e.getStart();
                carrierCount = 0;
                backgroundCount = 1;
            }
            if (carrierSamples.contains(e.getSample())) {
                carrierCount += e.getCount();
            } else {
                backgroundCount += e.getCount();
            }
        }
        return maxPosition;
    }

    private void writeCluster(final List<SVCallRecord> clusterVariants, final String clusterID) {
        for (final SVCallRecord variant : clusterVariants) {
            writer.print(callRecordCodec.encode(variant));
            writer.println(SVCallRecordCodec.COL_DELIMITER + clusterID);
        }
    }

    @Override
    public void traverse() {
        StreamSupport.stream(Spliterators.spliteratorUnknownSize(reader.iterator(), 0), false).forEachOrdered(v -> apply(v));
    }

    public void apply(final SVCallRecord variant) {
        if (!isValidSize(variant, minEventSize)) return;

        //TODO : only DEL for now
        if (!variant.getType().equals(StructuralVariantType.DEL)) return;

        final String contig = variant.getContig();
        if (currentContig == null || !contig.equals(currentContig)) {
            processBuffer();
            currentContig = contig;
            seedCluster(variant);
            return;
        }

        int clusterIndex = 0;
        final List<Integer> clustersToProcess = new ArrayList<>();
        boolean hasExistingCluster = false;
        int matchingClusterIndex = -1;
        for (final Tuple2<SimpleInterval, List<SVCallRecord>> cluster : currentClusters) {
            final SimpleInterval clusterInterval = cluster._1;
            final List<SVCallRecord> clusterVariants = cluster._2;
            if (variant.getStart() > clusterInterval.getEnd()) {
                clustersToProcess.add(clusterIndex);
            } else {
                boolean hasClusteringVariant = false;
                for (final SVCallRecord other : clusterVariants) {
                    if (clusterTogether(variant, other)) {
                        hasClusteringVariant = true;
                        matchingClusterIndex = clusterIndex;
                        break;
                    }
                }
                hasExistingCluster = hasExistingCluster || hasClusteringVariant;
            }
            clusterIndex++;
        }
        if (hasExistingCluster) {
            addToCluster(matchingClusterIndex, variant);
        } else {
            seedCluster(variant);
        }
        for (int i = clustersToProcess.size() - 1; i >= 0; i--) {
            processCluster(i);
        }
    }

    private boolean isValidSize(final Locatable loc, final int min) {
        return loc.getLengthOnReference() >= min;
    }


}
