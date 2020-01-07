package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import com.google.common.collect.Ordering;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalTree;
import htsjdk.samtools.util.OverlapDetector;
import htsjdk.variant.variantcontext.*;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.vcf.*;
import org.apache.commons.io.IOUtils;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.spark.sv.utils.GATKSVVCFConstants;
import org.broadinstitute.hellbender.utils.*;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import scala.Tuple2;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Paths;
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
 *         Unclustered structural variants from
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
    public static final String SAMPLE_COVERAGE_LONG_NAME = "sample-coverage";
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
            doc = "Sample coverage tsv",
            fullName = SAMPLE_COVERAGE_LONG_NAME
    )
    private String sampleCoverageFile;

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

    private FeatureDataSource<SVCallRecord> reader;
    private VariantContextWriter writer;

    private SVClusterEngine clusterEngine;
    private SplitReadEvidenceProcessor splitReadEvidenceProcessor;
    private SimpleInterval splitReadCacheInterval = null;
    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private SimpleInterval discordantPairCacheInterval = null;

    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;

    private Map<String,Double> sampleCoverageMap;
    private List<String> samplesList;
    private Set<String> samplesSet;

    private final Map<String,IntervalTree> whitelistedIntervalTreeMap = new HashMap<>();

    public static String END_CONTIG_ATTRIBUTE = "CHR2";
    public static String END_POS_ATTRIBUTE = VCFConstants.END_KEY;
    public static String SVLEN_ATTRIBUTE = GATKSVVCFConstants.SVLEN;
    public static String SVTYPE_ATTRIBUTE = VCFConstants.SVTYPE;
    public static String STRANDS_ATTRIBUTE = "STRANDS";
    public static String ALG_ATTRIBUTE = "ALG";
    public static String SPLIT_READ_START_COUNT_ATTRIBUTE = "SRS";
    public static String SPLIT_READ_END_COUNT_ATTRIBUTE = "SRE";
    public static String DISCORDANT_PAIR_COUNT_ATTRIBUTE = "PE";

    private final int SPLIT_READ_PADDING = 50;
    private final int SPLIT_READ_WINDOW = (SPLIT_READ_PADDING * 2) + 1;
    private final int MAX_INSERTION_SPLIT_READ_CROSS_DISTANCE = 20;
    private final int DISCORDANT_PAIR_PADDING = 500;

    private final int SPLIT_READ_QUERY_LOOKAHEAD = 0;
    private final int DISCORDANT_PAIR_QUERY_LOOKAHEAD = 0;
    private final int INPUT_QUERY_LOOKAHEAD = 10000;

    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
        loadSampleCoverage();
        loadSplitReadEvidenceDataSource();
        loadDiscordantPairDataSource();
        initializeWhitelistedIntervalTreeMap();

        clusterEngine = new SVClusterEngine(dictionary);
        splitReadEvidenceProcessor = new SplitReadEvidenceProcessor();
        reader = new FeatureDataSource<>(inputFile, "inputFile", INPUT_QUERY_LOOKAHEAD, SVCallRecord.class, getDefaultCloudPrefetchBufferSize(), getDefaultCloudIndexPrefetchBufferSize());
        reader.setIntervalsForTraversal(getTraversalIntervals());
        progressMeter.setRecordsBetweenTimeChecks(100);
        writer = createVCFWriter(Paths.get(outputFile));
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
                SPLIT_READ_QUERY_LOOKAHEAD,
                SplitReadEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void loadDiscordantPairDataSource() {
        discordantPairSource = new FeatureDataSource<>(
                discordantPairsFile,
                "discordantPairsFile",
                DISCORDANT_PAIR_QUERY_LOOKAHEAD,
                DiscordantPairEvidence.class,
                cloudPrefetchBuffer,
                cloudIndexPrefetchBuffer);
    }

    private void loadSampleCoverage() {
        try {
            sampleCoverageMap = IOUtils.readLines(BucketUtils.openFile(sampleCoverageFile), Charset.defaultCharset()).stream()
                    .map(line -> line.split("\t"))
                    .collect(Collectors.toMap(tokens -> tokens[0], tokens -> Double.valueOf(tokens[1])));
            samplesList = IOUtils.readLines(BucketUtils.openFile(sampleCoverageFile), Charset.defaultCharset()).stream()
                    .map(line -> line.split("\t")[0])
                    .collect(Collectors.toList());
            samplesSet = sampleCoverageMap.keySet();
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(sampleCoverageFile, e);
        }
    }


    private void initializeWhitelistedIntervalTreeMap() {
        for (final SimpleInterval interval : getTraversalIntervals()) {
            final String contig = interval.getContig();
            whitelistedIntervalTreeMap.putIfAbsent(contig, new IntervalTree<>());
            final IntervalTree tree = whitelistedIntervalTreeMap.get(contig);
            tree.put(interval.getStart(), interval.getEnd(), null);
        }
    }

    @Override
    public void traverse() {
        logger.info("Clustering raw calls...");
        final List<SVCallRecordWithEvidence> clusteredCalls = clusterCalls();
        logger.info("Processing start positions of " + clusteredCalls.size() + " clusters...");
        final List<SVCallRecordWithEvidence> startRefined = processStartPositions(clusteredCalls);
        logger.info("Processing end positions of " + clusteredCalls.size() + " clusters...");
        final List<SVCallRecordWithEvidence> endRefined = processEndPositions(startRefined);
        final List<SVCallRecordWithEvidence> finalCalls = SVClusterEngine.deduplicateCalls(endRefined, dictionary);
        logger.info("Writing output...");
        writeOutput(finalCalls);
    }

    private List<SVCallRecordWithEvidence> clusterCalls() {
        StreamSupport.stream(Spliterators.spliteratorUnknownSize(reader.iterator(), Spliterator.ORDERED), false)
                .filter(this::isValidSize)
                .filter(this::isWhitelisted)
                .map(SVCallRecordWithEvidence::new)
                .forEach(clusterEngine::addVariant);
        return clusterEngine.getOutput();
    }

    private List<SVCallRecordWithEvidence> processStartPositions(final List<SVCallRecordWithEvidence> calls) {
        final OverlapDetector splitReadStartIntervalOverlapDetector = getEvidenceOverlapDetector(calls, this::getStartSplitReadInterval);
        final OverlapDetector discordantPairIntervalOverlapDetector = getEvidenceOverlapDetector(calls, this::getDiscordantPairStartInterval);
        return calls.stream().map(c -> refineStartPosition(c, splitReadStartIntervalOverlapDetector, discordantPairIntervalOverlapDetector))
                .collect(Collectors.toList());
    }

    private List<SVCallRecordWithEvidence> processEndPositions(final List<SVCallRecordWithEvidence> calls) {
        final OverlapDetector splitReadEndIntervalOverlapDetector = getEvidenceOverlapDetector(calls, this::getEndSplitReadInterval);
        return calls.stream().sorted(Comparator.comparing(c -> c.getEndAsInterval(), IntervalUtils.getDictionaryOrderComparator(dictionary)))
                .map(c -> refineEndPosition(c, splitReadEndIntervalOverlapDetector))
                .collect(Collectors.toList());
    }

    private void writeOutput(final List<SVCallRecordWithEvidence> calls) {
        writeVCFHeader();
        calls.stream()
                .sorted(Comparator.comparing(c -> c.getStartAsInterval(), IntervalUtils.getDictionaryOrderComparator(dictionary)))
                .map(this::buildVariantContext)
                .forEachOrdered(writer::add);
    }

    private List<SplitReadEvidence> getSplitReads(final Function<SVCallRecord,SimpleInterval> intervalGetter,
                                                  final SVCallRecordWithEvidence call,
                                                  final OverlapDetector splitReadOverlapDetector) {
        final SimpleInterval interval = intervalGetter.apply(call);
        if (invalidCacheInterval(splitReadCacheInterval, interval)) {
            Utils.nonNull(splitReadOverlapDetector, "Split read cache missed but overlap detector is null");
            final Set<SimpleInterval> queryIntervalSet = splitReadOverlapDetector.getOverlaps(interval);
            if (queryIntervalSet.size() != 1) {
                throw new IllegalArgumentException("Call split read interval " + interval + " overlapped " + queryIntervalSet.size() + " query intervals");
            }
            splitReadCacheInterval = queryIntervalSet.iterator().next();
            splitReadSource.queryAndPrefetch(splitReadCacheInterval);
        }
        return splitReadSource.queryAndPrefetch(interval).stream()
                .filter(e -> e.getStrand() == call.getStartStrand())
                .collect(Collectors.toList());
    }

    private List<SplitReadEvidence> getStartSplitReads(final SVCallRecordWithEvidence call,
                                                       final OverlapDetector splitReadStartOverlapDetector) {
        return getSplitReads(this::getStartSplitReadInterval, call, splitReadStartOverlapDetector);
    }

    private SimpleInterval getStartSplitReadInterval(final SVCallRecord call) {
        return call.getStartAsInterval().expandWithinContig(SPLIT_READ_PADDING, dictionary);
    }

    private List<SplitReadEvidence> getEndSplitReads(final SVCallRecordWithEvidence call,
                                                     final OverlapDetector splitReadEndOverlapDetector) {
        return getSplitReads(this::getEndSplitReadInterval, call, splitReadEndOverlapDetector);
    }

    private SimpleInterval getEndSplitReadInterval(final SVCallRecord call) {
        final int lowerBound = getEndLowerBound(call);
        final SimpleInterval paddedInterval = call.getEndAsInterval().expandWithinContig(SPLIT_READ_PADDING, dictionary);
        return new SimpleInterval(paddedInterval.getContig(),
                Math.max(lowerBound, paddedInterval.getStart()),
                Math.max(lowerBound + 1, paddedInterval.getEnd()));
    }

    private boolean invalidCacheInterval(final SimpleInterval cacheInterval, final SimpleInterval queryInterval) {
        return cacheInterval == null
                || !queryInterval.getContig().equals(cacheInterval.getContig())
                || !queryInterval.spanWith(cacheInterval).equals(cacheInterval);
    }

    private List<DiscordantPairEvidence> getDiscordantPairs(final SVCallRecord call,
                                                            final OverlapDetector discordantPairStartOverlapDetector) {
        final SimpleInterval startInterval = getDiscordantPairStartInterval(call);
        if (invalidCacheInterval(discordantPairCacheInterval, startInterval)) {
            final Set<SimpleInterval> queryIntervalSet = discordantPairStartOverlapDetector.getOverlaps(startInterval);
            if (queryIntervalSet.size() != 1) {
                throw new IllegalArgumentException("Call end split read interval " + startInterval + " overlapped " + queryIntervalSet.size() + " query intervals");
            }
            discordantPairCacheInterval = queryIntervalSet.iterator().next();
            discordantPairSource.queryAndPrefetch(discordantPairCacheInterval);
        }
        final SimpleInterval endInterval = getDiscordantPairEndInterval(call);
        return discordantPairSource.queryAndPrefetch(startInterval).stream()
                .filter(e -> discordantPairOverlapsInterval(e, startInterval, endInterval))
                .filter(e -> e.getStartStrand() == call.getStartStrand() && e.getEndStrand() == call.getEndStrand())
                .collect(Collectors.toList());
    }

    private SimpleInterval getDiscordantPairStartInterval(final SVCallRecord call) {
        return call.getStartAsInterval().expandWithinContig(DISCORDANT_PAIR_PADDING, dictionary);
    }

    private SimpleInterval getDiscordantPairEndInterval(final SVCallRecord call) {
        return call.getEndAsInterval().expandWithinContig(DISCORDANT_PAIR_PADDING, dictionary);
    }

    private boolean discordantPairOverlapsInterval(final DiscordantPairEvidence evidence,
                                                   final SimpleInterval startInterval,
                                                   final SimpleInterval endInterval) {
        return evidence.getContig().equals(startInterval.getContig())
                && evidence.getStart() >= startInterval.getStart()
                && evidence.getStart() < startInterval.getEnd()
                && evidence.getEndContig().equals(endInterval.getContig())
                && evidence.getEnd() >= endInterval.getStart()
                && evidence.getEnd() < endInterval.getEnd();
    }


    private SVCallRecordWithEvidence refineStartPosition(final SVCallRecordWithEvidence call,
                                     final OverlapDetector splitReadStartIntervalOverlapDetector,
                                     final OverlapDetector discordantPairIntervalOverlapDetector) {
        final SVCallRecordWithEvidence refinedCall;
        if (SVClusterEngine.isDepthOnlyCall(call)) {
            refinedCall = new SVCallRecordWithEvidence(call, Collections.emptyList(), Collections.emptyList(), Collections.emptyList());
        } else {
            final List<SplitReadEvidence> startSplitReads = getStartSplitReads(call, splitReadStartIntervalOverlapDetector);
            final List<DiscordantPairEvidence> discordantPairs = getDiscordantPairs(call, discordantPairIntervalOverlapDetector);
            final Set<String> backgroundSamples = getBackgroundSamples(call);
            final SplitReadSite splitReadStartSite = splitReadEvidenceProcessor.process(call.getSamples(), backgroundSamples, startSplitReads, call.getStart(), sampleCoverageMap);
            final List<SplitReadSite> startSitesList = splitReadEvidenceProcessor.getSites();
            refinedCall = new SVCallRecordWithEvidence(
                    call.getContig(), splitReadStartSite.getPosition(), call.getStartStrand(), call.getEndContig(), call.getEnd(), call.getEndStrand(),
                    call.getType(), call.getLength(), call.getAlgorithms(), call.getSamples(), startSitesList, null, discordantPairs);
        }
        progressMeter.update(call.getStartAsInterval());
        return refinedCall;
    }

    private Set<String> getBackgroundSamples(final SVCallRecord call) {
        return samplesList.stream().filter(s -> call.getSamples().contains(s)).collect(Collectors.toSet());
    }

    // OverlapDetector may be null if the split reads are guaranteed to be cached
    private SVCallRecordWithEvidence refineEndPosition(final SVCallRecordWithEvidence call,
                                   final OverlapDetector splitReadEndIntervalOverlapDetector) {
        final SVCallRecordWithEvidence refinedCall;
        if (SVClusterEngine.isDepthOnlyCall(call)) {
            refinedCall = new SVCallRecordWithEvidence(call, Collections.emptyList(), Collections.emptyList(), Collections.emptyList());
        } else {
            final int endLowerBound = getEndLowerBound(call);
            final List<SplitReadEvidence> endSplitReads = getEndSplitReads(call, splitReadEndIntervalOverlapDetector);
            final int defaultEndPosition = Math.max(endLowerBound, call.getEnd());
            final Set<String> backgroundSamples = getBackgroundSamples(call);
            final SplitReadSite splitReadEndSite = splitReadEvidenceProcessor.process(call.getSamples(), backgroundSamples, endSplitReads, defaultEndPosition, sampleCoverageMap);
            final List<SplitReadSite> endSitesList = splitReadEvidenceProcessor.getSites();
            refinedCall = new SVCallRecordWithEvidence(
                    call.getContig(), call.getStart(), call.getStartStrand(), call.getEndContig(), splitReadEndSite.getPosition(), call.getEndStrand(),
                    call.getType(), call.getLength(), call.getAlgorithms(), call.getSamples(), call.getStartSplitReadSites(), endSitesList, call.getDiscordantPairs());
        }
        progressMeter.update(call.getEndAsInterval());
        return refinedCall;
    }

    private int getEndLowerBound(final SVCallRecord call) {
        return call.getType().equals(StructuralVariantType.INS) ?
                call.getStart() - MAX_INSERTION_SPLIT_READ_CROSS_DISTANCE :
                call.getStart() + 2; //TODO seems to need 2 instead of 1...
    }


    private final class SplitReadEvidenceProcessor {

        private final List<SplitReadSite> sites;
        private SplitReadSite selectedSite;

        public SplitReadEvidenceProcessor() {
            sites = new ArrayList<>(SPLIT_READ_WINDOW);
            reset();
        }

        void reset() {
            sites.clear();
            selectedSite = null;
        }

        private SplitReadSite process(final Set<String> carrierSamples,
                              final Set<String> backgroundSamples,
                              final List<SplitReadEvidence> splitReadEvidence,
                              final int defaultPosition,
                              final Map<String,Double> sampleCoverageMap) {
            if (!sampleCoverageMap.keySet().containsAll(carrierSamples)) {
                throw new IllegalArgumentException("One or more carrier samples not found in sample coverage map");
            }
            if (!sampleCoverageMap.keySet().containsAll(backgroundSamples)) {
                throw new IllegalArgumentException("One or more non-carrier samples not found in sample coverage map");
            }
            reset();
            computeSites(splitReadEvidence);
            selectedSite = backgroundSamples.isEmpty() ? testSitesByCount(defaultPosition, carrierSamples) : testSitesMannWhitney(defaultPosition, carrierSamples, backgroundSamples, sampleCoverageMap);
            if (selectedSite == null) {
                selectedSite = new SplitReadSite(defaultPosition, Collections.emptyMap());
            }
            return selectedSite;
        }

        private void computeSites(final List<SplitReadEvidence> evidenceList) {
            if (!Ordering.from(IntervalUtils.getDictionaryOrderComparator(dictionary)).isOrdered(evidenceList)) {
                throw new IllegalArgumentException("Evidence list is not dictionary sorted");
            }
            int position = 0;
            Map<String,Integer> sampleCounts = new HashMap<>();
            for (final SplitReadEvidence e : evidenceList) {
                if (e.getStart() != position) {
                    if (position > 0) {
                        sites.add(new SplitReadSite(position, sampleCounts));
                    }
                    position = e.getStart();
                    sampleCounts = new HashMap<>();
                }
                final String sample = e.getSample();
                sampleCounts.put(sample, e.getCount());
            }
        }

        private SplitReadSite testSitesByCount(final int defaultPosition,
                                      final Set<String> carrierSamples) {
            final OptionalInt maxCount = sites.stream().mapToInt(s -> s.getSampleCountSum(carrierSamples)).max();
            if (!maxCount.isPresent()) return null;
            return selectedSite = sites.stream()
                    .filter(s -> s.getSampleCountSum(carrierSamples) == maxCount.getAsInt())
                    .min(Comparator.comparingInt(s -> Math.abs(s.getPosition() - defaultPosition)))
                    .get();
        }

        private SplitReadSite testSitesMannWhitney(final int defaultPosition,
                               final Set<String> carrierSamples,
                               final Set<String> backgroundSamples,
                               final Map<String,Double> sampleCoverageMap) {
            final MannWhitneyU model = new MannWhitneyU();
            List<Tuple2<SplitReadSite,Double>> result = sites.stream()
                    .map(s -> new Tuple2<>(s, runTest(s, carrierSamples, backgroundSamples, sampleCoverageMap, model)))
                    .collect(Collectors.toList());
            final OptionalDouble maxP = result.stream().mapToDouble(Tuple2::_2).max();
            if (!maxP.isPresent()) return null;
            return selectedSite = result.stream()
                    .filter(t -> t._2 == maxP.getAsDouble())
                    .min(Comparator.comparingInt(t -> Math.abs(t._1.getPosition() - defaultPosition)))
                    .map(Tuple2::_1)
                    .get();

        }

        private double runTest(final SplitReadSite site,
                               final Set<String> carrierSamples,
                               final Set<String> backgroundSamples,
                               final Map<String,Double> sampleCoverageMap,
                               final MannWhitneyU model) {
            final double[] backgroundCounts = getNormalizedCounts(site, backgroundSamples, sampleCoverageMap);
            final double[] carrierCounts = getNormalizedCounts(site, carrierSamples, sampleCoverageMap);
            return model.test(backgroundCounts, carrierCounts, MannWhitneyU.TestType.FIRST_DOMINATES).getP();
        }

        private double[] getNormalizedCounts(final SplitReadSite site,
                                             final Set<String> samples,
                                             final Map<String,Double> sampleCoverageMap) {
            final double[] nonZeroCounts = site.getSampleCountsMap().entrySet().stream()
                    .filter(e -> samples.contains(e.getKey()))
                    .mapToDouble(e -> e.getValue() / sampleCoverageMap.get(e.getKey()))
                    .toArray();
            if (nonZeroCounts.length == samples.size()) return nonZeroCounts;
            final double[] counts = new double[samples.size()];
            for (int i = 0; i < nonZeroCounts.length; i++) {
                counts[i] = nonZeroCounts[i];
            }
            return counts;
        }

        public List<SplitReadSite> getSites() {
            return new ArrayList(sites);
        }

        public SplitReadSite getSelectedSite() {
            return selectedSite;
        }
    }

    private void writeVCFHeader() {
        final VCFHeader header = new VCFHeader(getDefaultToolVCFHeaderLines(), samplesList);
        header.setSequenceDictionary(dictionary);
        header.addMetaDataLine(new VCFInfoHeaderLine(END_CONTIG_ATTRIBUTE, 1, VCFHeaderLineType.String, "End contig"));
        header.addMetaDataLine(new VCFInfoHeaderLine(END_POS_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "End position"));
        header.addMetaDataLine(new VCFInfoHeaderLine(SVLEN_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Variant length"));
        header.addMetaDataLine(new VCFInfoHeaderLine(SVTYPE_ATTRIBUTE, 1, VCFHeaderLineType.String, "Variant type"));
        header.addMetaDataLine(new VCFInfoHeaderLine(ALG_ATTRIBUTE, 1, VCFHeaderLineType.String, "List of calling algorithms"));
        header.addMetaDataLine(new VCFFormatHeaderLine(SPLIT_READ_START_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at start of variant"));
        header.addMetaDataLine(new VCFFormatHeaderLine(SPLIT_READ_END_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Split read count at end of variant"));
        header.addMetaDataLine(new VCFFormatHeaderLine(DISCORDANT_PAIR_COUNT_ATTRIBUTE, 1, VCFHeaderLineType.Integer, "Discordant pair count"));
        writer.writeHeader(header);
    }

    public VariantContext buildVariantContext(final SVCallRecordWithEvidence call) {
        Utils.nonNull(call);
        Utils.nonNull(call.getStartSplitReadSites());
        Utils.nonNull(call.getEndSplitReadSites());
        Utils.nonNull(call.getDiscordantPairs());
        final Allele altAllele = Allele.create("<" + call.getType().name() + ">", false);
        final VariantContextBuilder builder = new VariantContextBuilder("", call.getContig(), call.getStart(), call.getEnd(),
                Lists.newArrayList(Allele.REF_N, altAllele));
        builder.attribute(END_CONTIG_ATTRIBUTE, call.getEndContig());
        builder.attribute(END_POS_ATTRIBUTE, call.getEnd());
        builder.attribute(SVLEN_ATTRIBUTE, call.getLength());
        builder.attribute(SVTYPE_ATTRIBUTE, call.getType());
        builder.attribute(ALG_ATTRIBUTE, call.getAlgorithms());
        final List<Genotype> genotypes = new ArrayList<>();
        final Map<String,Integer> startSplitReadCounts = getSplitReadCountsAtPosition(call.getStartSplitReadSites(), call.getStart());
        final Map<String,Integer> endSplitReadCounts = getSplitReadCountsAtPosition(call.getEndSplitReadSites(), call.getEnd());
        final Map<String,Integer> discordantPairCounts = getDiscordantPairCountsMap(call);
        for (final String sample : sampleCoverageMap.keySet()) {
            final GenotypeBuilder genotypeBuilder = new GenotypeBuilder(sample);
            genotypeBuilder.attribute(SPLIT_READ_START_COUNT_ATTRIBUTE, startSplitReadCounts.getOrDefault(sample, 0));
            genotypeBuilder.attribute(SPLIT_READ_END_COUNT_ATTRIBUTE, endSplitReadCounts.getOrDefault(sample, 0));
            genotypeBuilder.attribute(DISCORDANT_PAIR_COUNT_ATTRIBUTE, discordantPairCounts.getOrDefault(sample, 0));
            genotypes.add(genotypeBuilder.make());
        }
        builder.genotypes(genotypes);
        return builder.make();
    }

    private static Map<String,Integer> getSplitReadCountsAtPosition(final List<SplitReadSite> sites, final int pos) {
        return sites.stream()
                .filter(s -> s.getPosition() == pos)
                .map(SplitReadSite::getSampleCountsMap)
                .findAny()
                .orElse(Collections.emptyMap());
    }

    private Map<String,Integer> getDiscordantPairCountsMap(final SVCallRecordWithEvidence call) {
        return call.getDiscordantPairs().stream()
                .collect(Collectors.groupingBy(DiscordantPairEvidence::getSample,
                        Collectors.reducing(0, e -> 1, Integer::sum)));
    }

    private <T extends SVCallRecord> OverlapDetector getEvidenceOverlapDetector(final List<T> calls,
                                                       final Function<T,SimpleInterval> evidenceFunction) {
        final List<SimpleInterval> rawIntervals = calls.stream()
                .map(c -> evidenceFunction.apply(c))
                .sorted(IntervalUtils.getDictionaryOrderComparator(dictionary))
                .collect(Collectors.toList());
        final GenomeLocParser parser = new GenomeLocParser(dictionary);
        final List<GenomeLoc> rawLocs = IntervalUtils.genomeLocsFromLocatables(parser, rawIntervals);
        final List<GenomeLoc> mergedLocs = IntervalUtils.mergeIntervalLocations(rawLocs, IntervalMergingRule.ALL);
        final List<SimpleInterval> mergedIntervals = IntervalUtils.convertGenomeLocsToSimpleIntervals(mergedLocs);
        return OverlapDetector.create(mergedIntervals);
    }

    private boolean isValidSize(final SVCallRecord call) {
        return !((call.getType().equals(StructuralVariantType.DEL)
                || call.getType().equals(StructuralVariantType.DUP)
                || call.getType().equals(StructuralVariantType.INV))
                && call.getLength() < minEventSize);
    }

    private boolean isWhitelisted(final SVCallRecord call) {
        if (!whitelistedIntervalTreeMap.containsKey(call.getContig()) || !whitelistedIntervalTreeMap.containsKey(call.getEndContig())) {
            return false;
        }
        final IntervalTree startTree = whitelistedIntervalTreeMap.get(call.getContig());
        if (!startTree.overlappers(call.getStart(), call.getStart() + 1).hasNext()) {
            return false;
        }
        final IntervalTree endTree = whitelistedIntervalTreeMap.get(call.getEndContig());
        if (!endTree.overlappers(call.getEnd(), call.getEnd() + 1).hasNext()) {
            return false;
        }
        return true;
    }
}
