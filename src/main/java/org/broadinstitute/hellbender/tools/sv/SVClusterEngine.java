package org.broadinstitute.hellbender.tools.sv;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.variant.variantcontext.StructuralVariantType;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import scala.Tuple2;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class SVClusterEngine {

    private final double MIN_RECIPROCAL_OVERLAP_DEPTH = 0.8;
    private final double BREAKEND_CLUSTERING_WINDOW_FRACTION = 0.5;
    private final int MIN_BREAKEND_CLUSTERING_WINDOW = 50;
    private final int MAX_BREAKEND_CLUSTERING_WINDOW = 300;

    private final SAMSequenceDictionary dictionary;
    private final List<Tuple2<SimpleInterval,List<Integer>>> currentClusters;
    private final Map<Integer,SVCallRecordWithEvidence> idToVariantMap;
    private final List<SVCallRecordWithEvidence> outputBuffer;
    private int currentVariantId;
    private String currentContig;

    public SVClusterEngine(final SAMSequenceDictionary dictionary) {
        this.dictionary = dictionary;
        this.currentClusters = new LinkedList<>();
        this.idToVariantMap = new HashMap<>();
        this.outputBuffer = new ArrayList<>();
        currentVariantId = 0;
        currentContig = null;
    }

    public List<SVCallRecordWithEvidence> getOutput() {
        processBuffer();
        return deduplicateCalls(outputBuffer, dictionary);
    }

    public void addVariant(final SVCallRecordWithEvidence variant) {
        // Keep track of a unique id for each variant
        idToVariantMap.put(currentVariantId, variant);

        // Start a new cluster if on a new contig
        if (!variant.getContig().equals(currentContig)) {
            processBuffer();
            currentContig = variant.getContig();
            seedCluster(currentVariantId);
            return;
        }

        final List<Integer> clusterIdsToProcess = cluster(variant);
        processFinalizedClusters(clusterIdsToProcess);
        deleteRedundantClusters();
        currentVariantId++;
    }

    private List<Integer> cluster(final SVCallRecordWithEvidence variant) {
        // Get list of variant IDs from active clusters that cluster with this variant
        final Set<Integer> linkedVariantIds = idToVariantMap.entrySet().stream()
                .filter(other -> other.getKey().intValue() != currentVariantId && clusterTogether(variant, other.getValue()))
                .map(Map.Entry::getKey)
                .collect(Collectors.toSet());

        // Find clusters to which this variant belongs, and which active clusters we're definitely done with
        int clusterIndex = 0;
        final List<Integer> clusterIdsToProcess = new ArrayList<>();
        final List<Integer> fullyMatchingClusters = new ArrayList<>();
        final List<Integer> partiallyMatchingClusters = new ArrayList<>();
        for (final Tuple2<SimpleInterval, List<Integer>> cluster : currentClusters) {
            final SimpleInterval clusterInterval = cluster._1;
            final List<Integer> clusterVariantIds = cluster._2;
            if (variant.getStart() > clusterInterval.getEnd()) {
                clusterIdsToProcess.add(clusterIndex);
            } else {
                final int n = (int) clusterVariantIds.stream().filter(linkedVariantIds::contains).count();
                if (n == clusterVariantIds.size()) {
                    fullyMatchingClusters.add(clusterIndex);
                } else if (n > 0) {
                    partiallyMatchingClusters.add(clusterIndex);
                }
            }
            clusterIndex++;
        }

        // Add to variant cliques
        for (final int index : fullyMatchingClusters) {
            addToCluster(index, currentVariantId);
        }
        // Create new cliques
        for (final int index : partiallyMatchingClusters) {
            seedWithExistingCluster(currentVariantId, index, linkedVariantIds);
        }
        // If there weren't any matches, create a new singleton clique
        if (fullyMatchingClusters.isEmpty() && partiallyMatchingClusters.isEmpty()) {
            seedCluster(currentVariantId);
        }
        return clusterIdsToProcess;
    }

    private void processFinalizedClusters(final List<Integer> clusterIdsToProcess) {
        final Set<Integer> activeClusterIds = IntStream.range(0, currentClusters.size()).boxed().collect(Collectors.toSet());
        activeClusterIds.removeAll(clusterIdsToProcess);
        final Set<Integer> activeClusterVariantIds = activeClusterIds.stream().flatMap(i -> currentClusters.get(i)._2.stream()).collect(Collectors.toSet());
        final Set<Integer> finalizedVariantIds = clusterIdsToProcess.stream()
                .flatMap(i -> currentClusters.get(i)._2.stream())
                .filter(i -> !activeClusterVariantIds.contains(i))
                .collect(Collectors.toSet());
        for (int i = clusterIdsToProcess.size() - 1; i >= 0; i--) {
            processCluster(clusterIdsToProcess.get(i));
        }
        finalizedVariantIds.stream().forEach(idToVariantMap::remove);
    }

    private void deleteRedundantClusters() {
        final Set<Integer> redundantClusterSet = new HashSet<>();
        for (int i = 0; i < currentClusters.size(); i++) {
            final Set<Integer> clusterSetA = new HashSet<>(currentClusters.get(i)._2);
            for (int j = 0; j < i; j++) {
                final Set<Integer> clusterSetB = new HashSet<>(currentClusters.get(j)._2);
                if (clusterSetA.containsAll(clusterSetB)) {
                    redundantClusterSet.add(j);
                } else if (clusterSetA.size() != clusterSetB.size() && clusterSetB.containsAll(clusterSetA)) {
                    redundantClusterSet.add(i);
                }
            }
        }
        final List<Integer> redundantClustersList = new ArrayList<>(redundantClusterSet);
        redundantClustersList.sort(Comparator.naturalOrder());
        for (int i = redundantClustersList.size() - 1; i >= 0; i--) {
            currentClusters.remove((int)redundantClustersList.get(i));
        }
    }

    private void processCluster(final int clusterIndex) {
        if (clusterIndex < 0 || clusterIndex >= currentClusters.size()) {
            throw new IllegalArgumentException("Specified cluster index " + clusterIndex + " is out of range.");
        }
        final Tuple2<SimpleInterval, List<Integer>> cluster = currentClusters.remove(clusterIndex);
        final List<Integer> clusterVariantIds = cluster._2;
        if (clusterVariantIds.isEmpty()) {
            throw new IllegalArgumentException("Encountered empty cluster at index " + clusterIndex);
        }

        final List<SVCallRecordWithEvidence> clusterVariants = clusterVariantIds.stream().map(idToVariantMap::get).collect(Collectors.toList());
        final List<Integer> startPositions = clusterVariants.stream().map(SVCallRecordWithEvidence::getStart).sorted().collect(Collectors.toList());
        final List<Integer> endPositions = clusterVariants.stream().map(SVCallRecordWithEvidence::getEnd).sorted().collect(Collectors.toList());
        final int medianStart = startPositions.get(startPositions.size() / 2);
        final int medianEnd = endPositions.get(endPositions.size() / 2);
        final SVCallRecordWithEvidence exampleCall = clusterVariants.get(0);
        final int length = exampleCall.getContig().equals(exampleCall.getEndContig()) && !exampleCall.getType().equals(StructuralVariantType.INS) ? medianEnd - medianStart : exampleCall.getLength();
        final List<String> algorithms = clusterVariants.stream().flatMap(v -> v.getAlgorithms().stream()).distinct().collect(Collectors.toList());
        final Set<String> clusterSamples = clusterVariants.stream().flatMap(v -> v.getSamples().stream()).collect(Collectors.toSet());

        final int newStart;
        final int newEnd;
        if (exampleCall.getType().equals(StructuralVariantType.INS)) {
            // Insertions should be a single locus; also fixes case where end-supporting split reads are to the
            // left of start-supporting split reads
            final int mean = (medianStart + medianEnd) / 2;
            newStart = mean;
            newEnd = mean + 1;
        } else {
            newStart = medianStart;
            newEnd = medianEnd;
        }

        Utils.validate(!exampleCall.getContig().equals(exampleCall.getEndContig()) || newStart < newEnd,
                "Start position " + newStart + " comes after end position " + newEnd + " for " +
                        exampleCall.getType().name() + " cluster at approximately " + exampleCall.getContig() +
                        ":" + exampleCall.getStart() + ", " + exampleCall.getEndContig() + ":" + exampleCall.getEnd());
        outputBuffer.add(new SVCallRecordWithEvidence(exampleCall.getContig(), newStart, exampleCall.getStartStrand(),
                exampleCall.getEndContig(), newEnd, exampleCall.getEndStrand(), exampleCall.getType(), length, algorithms, clusterSamples,
                exampleCall.getStartSplitReadSites(), exampleCall.getEndSplitReadSites(), exampleCall.getDiscordantPairs()));
    }

    private boolean clusterTogether(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        if (!a.getType().equals(b.getType())) return false;
        final boolean depthOnlyA = isDepthOnlyCall(a);
        final boolean depthOnlyB = isDepthOnlyCall(b);
        if (depthOnlyA && depthOnlyB) {
            return clusterTogetherBothDepthOnly(a, b);
        } else if (!(depthOnlyA || depthOnlyB)) {
            return clusterTogetherBothWithEvidence(a, b);
        }
        // Mixed depth-only and non-depth-only
        return false;
    }

    private boolean clusterTogetherBothDepthOnly(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        if (!a.getContig().equals(a.getEndContig()) || !b.getContig().equals(b.getEndContig())) {
            throw new IllegalArgumentException("Attempted to cluster depth-only calls with endpoints on different contigs");
        }
        final SimpleInterval intervalA = new SimpleInterval(a.getContig(), a.getStart(), a.getEnd());
        final SimpleInterval intervalB = new SimpleInterval(b.getContig(), b.getStart(), b.getEnd());
        return IntervalUtils.isReciprocalOverlap(intervalA, intervalB, MIN_RECIPROCAL_OVERLAP_DEPTH);
    }

    private boolean clusterTogetherBothWithEvidence(final SVCallRecordWithEvidence a, final SVCallRecordWithEvidence b) {
        // Reject if one is intrachromosomal and the other isn't
        final boolean intrachromosomalA = a.getContig().equals(a.getEndContig());
        final boolean intrachromosomalB = b.getContig().equals(b.getEndContig());
        if (intrachromosomalA != intrachromosomalB) return false;

        // Matching endpoints
        final SimpleInterval intervalAStart =  getStartClusteringInterval(a);
        final SimpleInterval intervalAEnd =  getEndClusteringInterval(a);
        final SimpleInterval intervalBStart =  getStartClusteringInterval(b);
        final SimpleInterval intervalBEnd =  getEndClusteringInterval(b);
        return intervalAStart.overlaps(intervalBStart) && intervalAEnd.overlaps(intervalBEnd);
    }

    private SimpleInterval getStartClusteringInterval(final SVCallRecordWithEvidence call) {
        final int padding = getEndpointClusteringPadding(call);
        return call.getStartAsInterval().expandWithinContig(padding, dictionary);
    }

    private SimpleInterval getEndClusteringInterval(final SVCallRecordWithEvidence call) {
        final int padding =  getEndpointClusteringPadding(call);
        return call.getEndAsInterval().expandWithinContig(padding, dictionary);
    }

    private int getEndpointClusteringPadding(final SVCallRecordWithEvidence call) {
        return (int) Math.min(MAX_BREAKEND_CLUSTERING_WINDOW, Math.max(MIN_BREAKEND_CLUSTERING_WINDOW, BREAKEND_CLUSTERING_WINDOW_FRACTION * call.getLength()));
    }

    private SimpleInterval getClusteringStartInterval(final SVCallRecordWithEvidence variant, final SimpleInterval clusterMinStartInterval) {
        if (clusterMinStartInterval != null && !clusterMinStartInterval.getContig().equals(variant.getContig())) {
            throw new IllegalArgumentException("Attempted to get clustering start interval for variant on " + variant.getContig() + " but the given cluster start interval is " + clusterMinStartInterval);
        }
        if (!currentContig.equals(variant.getContig())) {
            throw new IllegalArgumentException("Attempted to get clustering start interval for variant on " + variant.getContig() + " but the current contig is " + currentClusters);
        }
        if (dictionary.getSequence(variant.getContig()) == null) {
            throw new IllegalArgumentException("Variant contig " + variant.getContig() + " not found in sequence dictionary.");
        }
        final int minStart;
        final int maxStart;
        if (isDepthOnlyCall(variant)) {
            minStart = (int) (variant.getStart() - (1.0 - MIN_RECIPROCAL_OVERLAP_DEPTH) * variant.getLength());
            maxStart = (int) (variant.getStart() + (1.0 - MIN_RECIPROCAL_OVERLAP_DEPTH) * variant.getLength());
        } else {
            minStart = variant.getStart() - MAX_BREAKEND_CLUSTERING_WINDOW;
            maxStart = variant.getStart() + MAX_BREAKEND_CLUSTERING_WINDOW;
        }
        if (clusterMinStartInterval == null) {
            return IntervalUtils.trimIntervalToContig(currentContig, minStart, maxStart, dictionary.getSequence(currentContig).getSequenceLength());
        }
        final int newMinStart = Math.max(minStart, clusterMinStartInterval.getStart());
        final int newMaxStart = Math.min(maxStart, clusterMinStartInterval.getEnd());
        return IntervalUtils.trimIntervalToContig(currentContig, newMinStart, newMaxStart, dictionary.getSequence(currentContig).getSequenceLength());
    }

    private void processBuffer() {
        while (!currentClusters.isEmpty()) {
            processCluster(0);
        }
    }

    private void seedCluster(final int seedId) {
        final SVCallRecordWithEvidence seed = idToVariantMap.get(seedId);
        if (seed == null) {
            throw new IllegalArgumentException("Variant id " + seedId + " not found in table");
        }
        if (!currentContig.equals(seed.getContig())) {
            throw new IllegalArgumentException("Attempted to seed new cluster with variant on contig " + seed.getContig() + " but the current contig is " + currentContig);
        }
        final List<Integer> newCluster = new ArrayList<>(1);
        newCluster.add(seedId);
        currentClusters.add(new Tuple2<>(getClusteringStartInterval(seed, null), newCluster));
    }

    private void seedWithExistingCluster(final int seedId, final int existingClusterIndex, final Set<Integer> clusteringIds) {
        final SVCallRecordWithEvidence seed = idToVariantMap.get(seedId);
        if (seed == null) {
            throw new IllegalArgumentException("Variant id " + seedId + " not found in table");
        }
        if (!currentContig.equals(seed.getContig())) {
            throw new IllegalArgumentException("Attempted to seed new cluster with variant on contig " + seed.getContig() + " but the current contig is " + currentContig);
        }
        final List<Integer> existingCluster = currentClusters.get(existingClusterIndex)._2;
        final List<Integer> validClusterIds = existingCluster.stream().filter(clusteringIds::contains).collect(Collectors.toList());
        final List<Integer> newCluster = new ArrayList<>(1 + existingCluster.size());
        newCluster.addAll(validClusterIds);
        newCluster.add(seedId);
        currentClusters.add(new Tuple2<>(getClusteringStartInterval(seed, currentClusters.get(existingClusterIndex)._1), newCluster));
    }

    private void addToCluster(final int clusterIndex, final int variantId) {
        final SVCallRecordWithEvidence variant = idToVariantMap.get(variantId);
        if (variant == null) {
            throw new IllegalArgumentException("Variant id " + variant + " not found in table");
        }
        if (!currentContig.equals(variant.getContig())) {
            throw new IllegalArgumentException("Attempted to add new variant on contig " + variant.getContig() + " but the current contig is " + currentContig);
        }
        if (clusterIndex >= currentClusters.size()) {
            throw new IllegalArgumentException("Specified cluster index " + clusterIndex + " is greater than the largest index.");
        }
        final Tuple2<SimpleInterval, List<Integer>> cluster = currentClusters.get(clusterIndex);
        final SimpleInterval clusterInterval = cluster._1;
        final List<Integer> clusterVariants = cluster._2;
        clusterVariants.add(variantId);
        final SimpleInterval clusteringStartInterval = getClusteringStartInterval(variant, clusterInterval);
        if (clusteringStartInterval.getStart() != clusterInterval.getStart() || clusteringStartInterval.getEnd() != clusterInterval.getEnd()) {
            currentClusters.remove(clusterIndex);
            currentClusters.add(clusterIndex, new Tuple2<>(clusteringStartInterval, clusterVariants));
        }
    }

    public static boolean isDepthOnlyCall(final SVCallRecordWithEvidence call) {
        for (final String alg : call.getAlgorithms()) {
            if (!alg.equals("depth")) return false;
        }
        return true;
    }

    public static List<SVCallRecordWithEvidence> deduplicateCalls(final List<SVCallRecordWithEvidence> calls,
                                                                  final SAMSequenceDictionary dictionary) {
        final List<SVCallRecordWithEvidence> sortedCalls = sortCallsByStart(calls, dictionary);
        final List<SVCallRecordWithEvidence> deduplicatedList = new ArrayList<>();
        int i = 0;
        while (i < sortedCalls.size()) {
            final SVCallRecordWithEvidence record = sortedCalls.get(i);
            int j = i + 1;
            final Collection<Integer> identicalCallIndexes = new ArrayList<>();
            while (j < sortedCalls.size() && record.getStartAsInterval().equals(sortedCalls.get(j).getStartAsInterval())) {
                final SVCallRecordWithEvidence other = sortedCalls.get(j);
                if (record.getEndAsInterval().equals(other.getEndAsInterval())
                        && record.getType().equals(other.getType())
                        && record.getStartStrand() == other.getStartStrand()
                        && record.getEndStrand() == other.getEndStrand()) {
                    identicalCallIndexes.add(j);
                }
                j++;
            }
            if (identicalCallIndexes.isEmpty()) {
                deduplicatedList.add(record);
                i++;
            } else {
                identicalCallIndexes.add(i);
                final List<SVCallRecordWithEvidence> identicalCalls = identicalCallIndexes.stream().map(sortedCalls::get).collect(Collectors.toList());
                final Set<String> samples = identicalCalls.stream()
                        .map(SVCallRecordWithEvidence::getSamples)
                        .flatMap(Collection::stream)
                        .collect(Collectors.toSet());
                final List<String> algorithms = identicalCalls.stream()
                        .map(SVCallRecordWithEvidence::getAlgorithms)
                        .flatMap(Collection::stream)
                        .distinct()
                        .collect(Collectors.toList());
                deduplicatedList.add(new SVCallRecordWithEvidence(
                        record.getContig(),
                        record.getStart(),
                        record.getStartStrand(),
                        record.getEndContig(),
                        record.getEnd(),
                        record.getEndStrand(),
                        record.getType(),
                        record.getLength(),
                        algorithms,
                        samples,
                        record.getStartSplitReadSites(),
                        record.getEndSplitReadSites(),
                        record.getDiscordantPairs()));
                i = j;
            }
        }
        return deduplicatedList;
    }

    public static <T extends SVCallRecord> List<T> sortCallsByStart(final Collection<T> calls,
                                                                    final SAMSequenceDictionary dictionary) {
        return calls.stream().sorted(IntervalUtils.getDictionaryOrderComparator(dictionary)).collect(Collectors.toList());
    }
}
