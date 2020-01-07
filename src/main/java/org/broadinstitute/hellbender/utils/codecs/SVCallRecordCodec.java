package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.tools.sv.SVCallRecord;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class SVCallRecordCodec extends AsciiFeatureCodec<SVCallRecord> {

    public static final String COL_DELIMITER = "\t";

    public SVCallRecordCodec() {
        super(SVCallRecord.class);
    }

    @Override
    public SVCallRecord decode(final String line) {
        final String[] tokens = line.split(COL_DELIMITER);
        return new SVCallRecord(
                tokens[0],
                Integer.parseUnsignedInt(tokens[1]),
                tokens[2].equals("+"),
                tokens[3],
                Integer.parseUnsignedInt(tokens[4]),
                tokens[5].equals("+"),
                StructuralVariantType.valueOf(tokens[6]),
                Integer.parseInt(tokens[7]),
                Arrays.asList(tokens[8]),
                Collections.singleton(tokens[9])
        );
    }

    @Override
    public TabixFormat getTabixFormat() {
        return TabixFormat.BED;
    }

    @Override
    public boolean canDecode(final String path) {
        return path.endsWith(".gz");
    } // TODO

    @Override
    public Object readActualHeader(final LineIterator reader) { return null; }

    public String encode(final SVCallRecord record) {
        final List<String> data = Arrays.asList(
                record.getContig(),
                Integer.toString(record.getStart()),
                record.getStartStrand() ? "+" : "-",
                record.getEndContig(),
                Integer.toString(record.getEnd()),
                record.getEndStrand() ? "+" : "-",
                record.getType().name(),
                Integer.toString(record.getLength()),
                String.join(",", record.getAlgorithms()),
                String.join(",", record.getSamples())
        );
        return String.join(COL_DELIMITER, data);
    }
}
