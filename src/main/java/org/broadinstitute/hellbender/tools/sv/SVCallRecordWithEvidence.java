package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.StructuralVariantType;

import java.util.List;
import java.util.Set;

public class SVCallRecordWithEvidence extends SVCallRecord {

    List<SplitReadSite> startSplitReadSites;
    List<SplitReadSite> endSplitReadSites;
    List<DiscordantPairEvidence> discordantPairs;

    public SVCallRecordWithEvidence(final SVCallRecord record) {
        super(record.getContig(), record.getStart(), record.getStartStrand(), record.getEndContig(), record.getEnd(),
                record.getEndStrand(), record.getType(), record.getLength(), record.getAlgorithms(), record.getSamples());
        this.startSplitReadSites = null;
        this.endSplitReadSites = null;
        this.discordantPairs = null;
    }

    public SVCallRecordWithEvidence(final SVCallRecord record,
                                    List<SplitReadSite> startSplitReadSites,
                                    List<SplitReadSite> endSplitReadSites,
                                    List<DiscordantPairEvidence> discordantPairs) {
        super(record.getContig(), record.getStart(), record.getStartStrand(), record.getEndContig(), record.getEnd(),
                record.getEndStrand(), record.getType(), record.getLength(), record.getAlgorithms(), record.getSamples());
        this.startSplitReadSites = startSplitReadSites;
        this.endSplitReadSites = endSplitReadSites;
        this.discordantPairs = discordantPairs;
    }

    public SVCallRecordWithEvidence(String startContig,
                        int start,
                        boolean startStrand,
                        String endContig,
                        int end,
                        boolean endStrand,
                        StructuralVariantType type,
                        int length,
                        List<String> algorithms,
                        Set<String> samples,
                        List<SplitReadSite> startSplitReadSites,
                        List<SplitReadSite> endSplitReadSites,
                        List<DiscordantPairEvidence> discordantPairs) {
        super(startContig, start, startStrand, endContig, end, endStrand, type, length, algorithms, samples);
        this.startSplitReadSites = startSplitReadSites;
        this.endSplitReadSites = endSplitReadSites;
        this.discordantPairs = discordantPairs;
    }

    public void setDiscordantPairs(final List<DiscordantPairEvidence> discordantPairs) {
        this.discordantPairs = discordantPairs;
    }

    public List<DiscordantPairEvidence> getDiscordantPairs() {
        return discordantPairs;
    }

    public List<SplitReadSite> getStartSplitReadSites() {
        return startSplitReadSites;
    }

    public List<SplitReadSite> getEndSplitReadSites() {
        return endSplitReadSites;
    }

    public void setStartSplitReadSites(final List<SplitReadSite> startSplitReadSites) {
        this.startSplitReadSites = startSplitReadSites;
    }

    public void setEndSplitReadSites(final List<SplitReadSite> endSplitReadSites) {
        this.endSplitReadSites = endSplitReadSites;
    }
}
