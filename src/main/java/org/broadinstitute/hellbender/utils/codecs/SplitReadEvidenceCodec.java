package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.tools.sv.SplitReadEvidence;

public class SplitReadEvidenceCodec extends AsciiFeatureCodec<SplitReadEvidence> {

    public static final String COL_DELIMITER = "\t";

    public SplitReadEvidenceCodec() {
        super(SplitReadEvidence.class);
    }

    @Override
    public SplitReadEvidence decode(final String line) {
        final String[] tokens = line.split(COL_DELIMITER);
        final String contig = tokens[0];
        final int position = Integer.parseUnsignedInt(tokens[1]);
        final boolean strand = tokens[2].equals("right");
        final int count = Integer.parseUnsignedInt(tokens[3]);
        final String sample = tokens[4];
        return new SplitReadEvidence(sample, contig, position, count, strand);
    }


    @Override
    public TabixFormat getTabixFormat() {
        return TabixFormat.BED;
    }

    @Override
    public boolean canDecode(final String path) {
        return true;
    } // TODO

    @Override
    public Object readActualHeader(final LineIterator reader) { return null; }
}
