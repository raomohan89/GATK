package org.broadinstitute.hellbender.tools.sv;

import java.util.Map;
import java.util.Set;

final class SplitReadSite {
    private final int position;
    private final Map<String,Integer> sampleCountsMap;

    public SplitReadSite(final int position, final Map<String,Integer> sampleCountsMap) {
        this.position = position;
        this.sampleCountsMap = sampleCountsMap;
    }

    public int getPosition() {
        return position;
    }

    public Map<String,Integer> getSampleCountsMap() {
        return sampleCountsMap;
    }

    public int getSampleCountSum(final Set<String> samples) {
        return sampleCountsMap.entrySet().stream()
                .filter(e -> samples.contains(e.getKey()))
                .mapToInt(e -> e.getValue())
                .sum();
    }

    public int getCountSum() {
        return sampleCountsMap.values().stream().mapToInt(Integer::intValue).sum();
    }
}
