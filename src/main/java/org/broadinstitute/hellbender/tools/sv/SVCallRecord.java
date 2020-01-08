package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.Arrays;
import java.util.List;
import java.util.Set;
import java.util.stream.Collectors;

public class SVCallRecord implements Feature {

    private final String startContig;
    private final int start;
    private final boolean startStrand;
    private final String endContig;
    private final int end;
    private final boolean endStrand;
    private final StructuralVariantType type;
    private int length;
    private final List<String> algorithms;
    private final Set<String> samples;

    private final static List<String> attributeStrings = Arrays.asList(
            SVCluster.END_CONTIG_ATTRIBUTE,
            SVCluster.ALG_ATTRIBUTE,
            SVCluster.STRANDS_ATTRIBUTE,
            SVCluster.SVLEN_ATTRIBUTE
    );

    public static SVCallRecord create(final VariantContext variant) {
        Utils.nonNull(variant);
        for (final String attr : attributeStrings) {
            if (!variant.hasAttribute(attr)) {
                throw new IllegalArgumentException("Attribute not found: " + attr);
            }
        }
        final String startContig = variant.getContig();
        final int start = variant.getStart();
        final String endContig = variant.getAttributeAsString(SVCluster.END_CONTIG_ATTRIBUTE, "NA");
        final int end = variant.getEnd();
        final StructuralVariantType type = variant.getStructuralVariantType();
        final List<String> algorithms = variant.getAttributeAsStringList(SVCluster.ALG_ATTRIBUTE, "NA");
        final String strands = variant.getAttributeAsString(SVCluster.STRANDS_ATTRIBUTE, "0");
        final char startStrandChar = strands.charAt(0);
        if (startStrandChar != '+' && startStrandChar != '-') {
            throw new IllegalArgumentException("Valid start strand not found");
        }
        final char endStrandChar = strands.charAt(1);
        if (endStrandChar != '+' && endStrandChar != '-') {
            throw new IllegalArgumentException("Valid end strand not found");
        }
        final boolean startStrand = strands.charAt(0) == '+';
        final boolean endStrand = strands.charAt(1) == '+';
        final int length = variant.getAttributeAsInt(SVCluster.SVLEN_ATTRIBUTE, 0);
        final Set<String> samples = variant.getGenotypes().stream()
                .filter(Genotype::isCalled)
                .map(Genotype::getSampleName)
                .collect(Collectors.toSet());
        return new SVCallRecord(startContig, start, startStrand, endContig, end, endStrand, type, length, algorithms, samples);
    }

    public SVCallRecord(final String startContig,
                        final int start,
                        final boolean startStrand,
                        final String endContig,
                        final int end,
                        final boolean endStrand,
                        final StructuralVariantType type,
                        final int length,
                        final List<String> algorithms,
                        final Set<String> samples) {
        Utils.nonNull(startContig);
        Utils.nonNull(endContig);
        Utils.nonNull(type);
        Utils.nonNull(algorithms);
        Utils.nonNull(samples);
        Utils.nonEmpty(algorithms);
        Utils.nonEmpty(samples);
        Utils.containsNoNull(algorithms, "Encountered null algorithm");
        Utils.containsNoNull(samples, "Encountered null sample");
        this.startContig = startContig;
        this.start = start;
        this.startStrand = startStrand;
        this.endContig = endContig;
        this.end = end;
        this.endStrand = endStrand;
        this.type = type;
        this.length = length;
        this.algorithms = algorithms;
        this.samples = samples;
    }

    @Override
    public String getContig() {
        return startContig;
    }

    @Override
    public int getStart() {
        return start;
    }

    public boolean getStartStrand() {
        return startStrand;
    }

    public String getEndContig() {
        return endContig;
    }

    @Override
    public int getEnd() {
        return end;
    }

    public boolean getEndStrand() {
        return endStrand;
    }

    public StructuralVariantType getType() {
        return type;
    }

    public int getLength() {
        return length;
    }

    public List<String> getAlgorithms() {
        return algorithms;
    }

    public Set<String> getSamples() {
        return samples;
    }

    public SimpleInterval getStartAsInterval() {
        return new SimpleInterval(startContig, start, start + 1);
    }

    public SimpleInterval getEndAsInterval() {
        return new SimpleInterval(endContig, end, end + 1);
    }
}
