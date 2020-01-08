package org.broadinstitute.hellbender.tools.sv;

import com.google.common.collect.Lists;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.IntervalTree;
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
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Paths;
import java.util.*;
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

    private FeatureDataSource<SplitReadEvidence> splitReadSource;
    private FeatureDataSource<DiscordantPairEvidence> discordantPairSource;

    private SVClusterEngine clusterEngine;
    private BreakpointRefiner breakpointRefiner;
    private SVEvidenceCollector evidenceCollector;
    private Map<String,Double> sampleCoverageMap;
    private List<String> samplesList;

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
        breakpointRefiner = new BreakpointRefiner(sampleCoverageMap);
        evidenceCollector = new SVEvidenceCollector(splitReadSource, discordantPairSource, dictionary, progressMeter);

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
        logger.info("Collecting evidence for " + clusteredCalls.size() + " clusters...");
        final List<SVCallRecordWithEvidence> callsWithEvidence = evidenceCollector.collectEvidence(clusteredCalls);
        logger.info("Refining breakpoints of " + clusteredCalls.size() + " clusters...");
        final List<SVCallRecordWithEvidence> refinedCalls = callsWithEvidence.stream().map(breakpointRefiner::refineCalls).collect(Collectors.toList());
        final List<SVCallRecordWithEvidence> finalCalls = SVClusterEngine.deduplicateCalls(refinedCalls, dictionary);
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

    private void writeOutput(final List<SVCallRecordWithEvidence> calls) {
        writeVCFHeader();
        calls.stream()
                .sorted(Comparator.comparing(c -> c.getStartAsInterval(), IntervalUtils.getDictionaryOrderComparator(dictionary)))
                .map(this::buildVariantContext)
                .forEachOrdered(writer::add);
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
        final Map<String,Integer> discordantPairCounts = getDiscordantPairCountsMap(call.getDiscordantPairs());
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
        Utils.nonNull(sites);
        Utils.validateArg(pos > 0, "Non-positive position");
        if (sites.stream().map(SplitReadSite::getPosition).distinct().count() != sites.size()) {
            throw new IllegalArgumentException("Sites did not have unique positions");
        }
        return sites.stream()
                .filter(s -> s.getPosition() == pos)
                .map(SplitReadSite::getSampleCountsMap)
                .findAny()
                .orElse(Collections.emptyMap());
    }

    private Map<String,Integer> getDiscordantPairCountsMap(final Collection<DiscordantPairEvidence> discordantPairs) {
        Utils.nonNull(discordantPairs);
        return discordantPairs.stream()
                .collect(Collectors.groupingBy(DiscordantPairEvidence::getSample,
                        Collectors.reducing(0, e -> 1, Integer::sum)));
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
