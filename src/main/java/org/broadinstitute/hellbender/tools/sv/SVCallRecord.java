package org.broadinstitute.hellbender.tools.sv;

import htsjdk.tribble.Feature;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.StructuralVariantType;
import htsjdk.variant.variantcontext.VariantContext;
import org.broadinstitute.hellbender.utils.SimpleInterval;

import java.util.Collection;
import java.util.List;
import java.util.stream.Collectors;

public class SVCallRecord implements Feature {

    private final String startContig;
    private final int start;
    private final boolean startStrand;

    private final String endContig;
    private final int end;
    private final boolean endStrand;

    private final StructuralVariantType type;
    private final int length;
    private final List<String> algorithms;

    private final String sample;

    public static Collection<SVCallRecord> create(final VariantContext variant) {
        final String startContig = variant.getContig();
        final int start = variant.getStart();
        final String endContig = variant.getAttributeAsString("CHR2", "NA");
        final int end = variant.getEnd();
        final StructuralVariantType type = variant.getStructuralVariantType();
        final List<String> algorithms = variant.getAttributeAsStringList("ALGORITHMS", "NA");
        final String strands = variant.getAttributeAsString("STRANDS", "XX");
        final boolean startStrand = strands.charAt(0) == '+';
        final boolean endStrand = strands.charAt(1) == '+';
        final int length = variant.getAttributeAsInt("SVLEN", 0);
        return variant.getGenotypes().stream()
                .filter(Genotype::isCalled)
                .map(Genotype::getSampleName)
                .map(s -> new SVCallRecord(startContig, start, startStrand, endContig, end, endStrand, type, length, algorithms, s))
                .collect(Collectors.toList());
    }

    public SVCallRecord(String startContig, int start, boolean startStrand, String endContig, int end, boolean endStrand, StructuralVariantType type, int length, List<String> algorithms, String sample) {
        this.startContig = startContig;
        this.start = start;
        this.startStrand = startStrand;
        this.endContig = endContig;
        this.end = end;
        this.endStrand = endStrand;
        this.type = type;
        this.length = length;
        this.algorithms = algorithms;
        this.sample = sample;
    }

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

    public String getSample() {
        return sample;
    }

    public SimpleInterval getStartAsInterval() {
        return new SimpleInterval(startContig, start, start + 1);
    }

    public SimpleInterval getEndAsInterval() {
        return new SimpleInterval(endContig, end, end + 1);
    }
}
