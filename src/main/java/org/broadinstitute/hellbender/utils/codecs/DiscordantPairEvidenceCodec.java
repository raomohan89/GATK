package org.broadinstitute.hellbender.utils.codecs;

import htsjdk.tribble.AsciiFeatureCodec;
import htsjdk.tribble.index.tabix.TabixFormat;
import htsjdk.tribble.readers.LineIterator;
import org.broadinstitute.hellbender.tools.sv.DiscordantPairEvidence;

public class DiscordantPairEvidenceCodec extends AsciiFeatureCodec<DiscordantPairEvidence> {

    public static final String COL_DELIMITER = "\t";

    public DiscordantPairEvidenceCodec() {
        super(DiscordantPairEvidence.class);
    }

    @Override
    public DiscordantPairEvidence decode(final String line) {
        final String[] tokens = line.split(COL_DELIMITER);
        final String startContig = tokens[0];
        final int start = Integer.parseUnsignedInt(tokens[1]);
        final boolean startStrand = tokens[2].equals("+");
        final String endContig = tokens[3];
        final int end = Integer.parseUnsignedInt(tokens[4]);
        final boolean endStrand = tokens[5].equals("+");
        final String sample = tokens[6];
        return new DiscordantPairEvidence(sample, startContig, start, startStrand, endContig, end, endStrand);
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
