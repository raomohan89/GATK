package org.broadinstitute.hellbender.tools.sv;

import htsjdk.variant.variantcontext.StructuralVariantType;
import org.apache.commons.math3.util.FastMath;
import org.apache.commons.math3.util.MathUtils;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;

public class BreakpointRefiner {

    private final Map<String,Double> sampleCoverageMap;
    private int maxInsertionSplitReadCrossDistance;
    private double minSiteProbability;

    public static final int DEFAULT_MAX_INSERTION_CROSS_DISTANCE = 20;
    public static final double DEFAULT_MIN_SITE_PROBABILITY = 0.99;
    private final static double[] FACTORIALS = {1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880, 3628800, 39916800, 479001600};

    public BreakpointRefiner(final Map<String,Double> sampleCoverageMap) {
        this.sampleCoverageMap = sampleCoverageMap;
        this.maxInsertionSplitReadCrossDistance = DEFAULT_MAX_INSERTION_CROSS_DISTANCE;
        minSiteProbability = DEFAULT_MIN_SITE_PROBABILITY;
    }

    public void setMaxInsertionSplitReadCrossDistance(final int distance) {
        maxInsertionSplitReadCrossDistance = distance;
    }

    public void setMinSiteProbability(final double p) {
        minSiteProbability = p;
    }

    public SVCallRecordWithEvidence refineCalls(final SVCallRecordWithEvidence call) {
        Utils.nonNull(call);
        Utils.nonNull(call.getStartSplitReadSites());
        final SVCallRecordWithEvidence refinedCall;
        if (SVClusterEngine.isDepthOnlyCall(call)) {
            refinedCall = call;
        } else {
            final Set<String> backgroundSamples = getBackgroundSamples(call);
            final SplitReadSite refinedStartSite = getRefinedSite(call.getStartSplitReadSites(), call.getSamples(), backgroundSamples, call.getStart());
            final int endLowerBound = getEndLowerBound(call.getType(), call.getContig(), refinedStartSite.getPosition(), call.getEndContig());
            final int defaultEndPosition = Math.max(endLowerBound, call.getEnd());
            final List<SplitReadSite> validSites = getValidEndSplitReadSites(call, endLowerBound);
            final SplitReadSite splitReadEndSite = getRefinedSite(validSites, call.getSamples(), backgroundSamples, defaultEndPosition);
            refinedCall = new SVCallRecordWithEvidence(
                    call.getContig(), refinedStartSite.getPosition(), call.getStartStrand(),
                    call.getEndContig(), splitReadEndSite.getPosition(), call.getEndStrand(),
                    call.getType(), call.getLength(), call.getAlgorithms(), call.getSamples(),
                    call.getStartSplitReadSites(), call.getEndSplitReadSites(), call.getDiscordantPairs());
        }
        return refinedCall;
    }

    private List<SplitReadSite> getValidEndSplitReadSites(final SVCallRecordWithEvidence call, final int endLowerBound) {
        return call.getEndSplitReadSites().stream()
                .filter(s -> s.getPosition() >= endLowerBound)
                .collect(Collectors.toList());
    }

    private Set<String> getBackgroundSamples(final SVCallRecord call) {
        return sampleCoverageMap.keySet().stream().filter(s -> !call.getSamples().contains(s)).collect(Collectors.toSet());
    }

    private int getEndLowerBound(final StructuralVariantType type, final String startContig, final int startPosition, final String endContig) {
        if (!startContig.equals(endContig)) return 0;
        return type.equals(StructuralVariantType.INS) ?
                startPosition - maxInsertionSplitReadCrossDistance :
                startPosition + 2; //TODO seems to need 2 instead of 1...
    }

    private SplitReadSite getRefinedSite(final List<SplitReadSite> sites,
                                         final Set<String> carrierSamples,
                                         final Set<String> backgroundSamples,
                                         final int defaultPosition) {
        if (!sampleCoverageMap.keySet().containsAll(carrierSamples)) {
            throw new IllegalArgumentException("One or more carrier samples not found in sample coverage map");
        }
        if (!sampleCoverageMap.keySet().containsAll(backgroundSamples)) {
            throw new IllegalArgumentException("One or more non-carrier samples not found in sample coverage map");
        }
        return testSitesByLikelihood(sites, defaultPosition, carrierSamples, backgroundSamples);
    }

    private SplitReadSite testSitesByRateExcess(final List<SplitReadSite> sites,
                                                final int defaultPosition,
                                                final Set<String> carrierSamples,
                                                final Set<String> backgroundSamples) {
        if (sites.isEmpty() || carrierSamples.isEmpty()) return new SplitReadSite(defaultPosition, Collections.emptyMap());
        final double carrierCoverage = carrierSamples.stream().mapToDouble(sampleCoverageMap::get).sum();
        final double backgroundCoverage = backgroundSamples.stream().mapToDouble(sampleCoverageMap::get).sum();
        List<Tuple2<SplitReadSite,Double>> siteScorePairs = sites.stream()
                .map(s -> new Tuple2<>(s, calculateRateExcess(s, carrierSamples, backgroundSamples, carrierCoverage, backgroundCoverage)))
                .collect(Collectors.toList());
        final double maxScore = siteScorePairs.stream().mapToDouble(Tuple2::_2).max().getAsDouble();
        SplitReadSite result = siteScorePairs.stream()
                .filter(p -> p._2 == maxScore)
                .min(Comparator.comparingInt(s -> Math.abs(s._1.getPosition() - defaultPosition)))
                .get()._1;
        return result;
    }

    private double calculateRateExcess(final SplitReadSite site,
                                       final Set<String> carrierSamples,
                                       final Set<String> backgroundSamples,
                                       final double carrierCoverage,
                                       final double backgroundCoverage) {
        final double carrierRate = site.getSampleCountSum(carrierSamples) / carrierCoverage;
        final double backgroundRate = backgroundSamples.isEmpty() ? 0. : site.getSampleCountSum(backgroundSamples) / backgroundCoverage;
        return Math.max(carrierRate - backgroundRate, 0.);
    }

    private SplitReadSite testSitesByCount(final List<SplitReadSite> sites,
                                           final int defaultPosition,
                                           final Set<String> carrierSamples) {
        final OptionalInt maxCount = sites.stream().mapToInt(s -> s.getSampleCountSum(carrierSamples)).max();
        if (!maxCount.isPresent()) return null;
        return sites.stream()
                .filter(s -> s.getSampleCountSum(carrierSamples) == maxCount.getAsInt())
                .min(Comparator.comparingInt(s -> Math.abs(s.getPosition() - defaultPosition)))
                .get();
    }

    private SplitReadSite testSitesByLikelihood(final List<SplitReadSite> sites,
                                                final int defaultPosition,
                                                final Set<String> carrierSamples,
                                                final Set<String> backgroundSamples) {
        if (sites.isEmpty() || carrierSamples.isEmpty()) return new SplitReadSite(defaultPosition, Collections.emptyMap());
        final double epsilon = sites.stream().mapToDouble(s -> s.getNormalizedCountSum(sampleCoverageMap)).sum() / (sites.size() * (carrierSamples.size() + backgroundSamples.size()));
        List<Tuple2<SplitReadSite,Double>> siteScorePairs = sites.stream()
                .map(s -> new Tuple2<>(s, calculateJointLikelihood(s, carrierSamples, backgroundSamples, epsilon)))
                .collect(Collectors.toList());
        final double maxScore = siteScorePairs.stream().mapToDouble(Tuple2::_2).max().getAsDouble();
        final double maxProbability = 1. / siteScorePairs.stream().mapToDouble(p -> FastMath.exp(p._2 - maxScore)).sum();
        if (maxProbability < minSiteProbability) return new SplitReadSite(defaultPosition, Collections.emptyMap());
        return siteScorePairs.stream()
                .filter(p -> p._2 == maxScore)
                .min(Comparator.comparingInt(s -> Math.abs(s._1.getPosition() - defaultPosition)))
                .get()._1;
    }

    private double calculateJointLikelihood(final SplitReadSite site,
                                            final Set<String> carrierSamples,
                                            final Set<String> backgroundSamples,
                                            final double epsilon) {
        final double phi = 0.35;
        final int N = 2;
        final Map<String,Integer> sampleCountMap = site.getSampleCountsMap();
        double logLik = 1;
        for (final String sample : carrierSamples) {
            double pSample = 0;
            for (int j = 1; j <= N; j++) {
                final double rate = sampleCoverageMap.get(sample) * j * phi;
                final int count = sampleCountMap.containsKey(sample) ? sampleCountMap.get(sample) : 0;
                pSample += FastMath.exp(logPoissonProbability(rate, count));
            }
            logLik += FastMath.log(pSample);
        }
        for (final String sample : backgroundSamples) {
            final double rate = sampleCoverageMap.get(sample) * epsilon;
            final int count = sampleCountMap.containsKey(sample) ? sampleCountMap.get(sample) : 0;
            logLik += logPoissonProbability(rate, count);
        }
        return logLik;
    }

    private static double logPoissonProbability(double mean, int x) {
        Utils.validateArg(mean > 0, "Non-positive mean");
        Utils.validateArg(x >= 0, "Negative count");
        if (x == 0) {
            return -mean;
        } else if (x < FACTORIALS.length) {
            return x * FastMath.log(mean) - mean - FastMath.log(FACTORIALS[x]);
        }
        return x * FastMath.log(mean) - mean - 0.5 * FastMath.log(MathUtils.TWO_PI) - 0.5 * FastMath.log(x) - x * FastMath.log(x) + x;
    }
}
