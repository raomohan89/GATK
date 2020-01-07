package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.StructuralVariantDiscoveryProgramGroup;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.codecs.SVCallRecordCodec;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.Spliterators;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

/**
 * Creates sparsely formatted file of structural variants.
 *
 * <h3>Inputs</h3>
 *
 * <ul>
 *     <li>
 *         SV VCFs
 *     </li>
 * </ul>
 *
 * <h3>Output</h3>
 *
 * <ul>
 *     <li>
 *         Sparse variants file
 *     </li>
 * </ul>
 *
 * <h3>Usage example</h3>
 *
 * <pre>
 *     gatk SVCreateSparseVariants
 * </pre>
 *
 * @author Mark Walker &lt;markw@broadinstitute.org&gt;
 */

@CommandLineProgramProperties(
        summary = "Creates sparse structural variants file",
        oneLineSummary = "Creates sparse structural variants file",
        programGroup = StructuralVariantDiscoveryProgramGroup.class
)
@DocumentedFeature
public final class CreateSVCallRecords extends GATKTool {
    public static final String IGNORE_DICTIONARY_LONG_NAME = "ignore-dict";

    @Argument(
            doc = "Input VCFs",
            fullName = StandardArgumentDefinitions.VARIANT_LONG_NAME,
            shortName = StandardArgumentDefinitions.VARIANT_SHORT_NAME
    )
    private List<String> inputFiles;

    @Argument(
            doc = "Output VCF",
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private File outputFile;

    @Argument(
            doc = "Skip VCF sequence dictionary check",
            fullName = IGNORE_DICTIONARY_LONG_NAME,
            optional = true
    )
    private Boolean skipVcfDictionaryCheck = false;

    private List<SVCallRecord> records;
    private SAMSequenceDictionary dictionary;


    @Override
    public void onTraversalStart() {
        dictionary = getBestAvailableSequenceDictionary();
        if (dictionary == null) {
            throw new UserException("Reference sequence dictionary required");
        }
    }

    @Override
    public Object onTraversalSuccess() {
        return null;
    }

    @Override
    public void traverse() {
        records = new ArrayList<>();
        int index = 0;
        for (final String file : inputFiles) {
            records.addAll(StreamSupport.stream(Spliterators.spliteratorUnknownSize(getVariants(file, "vcf_" + index), 0), false)
                    .map(SVCallRecord::create).collect(Collectors.toList()));
        }
        records.sort(IntervalUtils.getDictionaryOrderComparator(dictionary));
        writeVariants();
    }

    private Iterator<VariantContext> getVariants(final String vcf, final String name) {
        final FeatureDataSource<VariantContext> featureDataSource = new FeatureDataSource<>(vcf, name, 100000, VariantContext.class, cloudPrefetchBuffer, cloudIndexPrefetchBuffer);
        if (!skipVcfDictionaryCheck) {
            featureDataSource.getSequenceDictionary().assertSameDictionary(dictionary);
        }
        featureDataSource.setIntervalsForTraversal(getTraversalIntervals());
        return featureDataSource.iterator();
    }

    private void writeVariants() {
        final SVCallRecordCodec callRecordCodec = new SVCallRecordCodec();
        try (final PrintStream writer = new PrintStream(outputFile)) {
            records.stream().forEachOrdered(r -> writer.println(callRecordCodec.encode(r)));
        } catch(final IOException e) {
            throw new GATKException("Error writing output file", e);
        }
    }
}
