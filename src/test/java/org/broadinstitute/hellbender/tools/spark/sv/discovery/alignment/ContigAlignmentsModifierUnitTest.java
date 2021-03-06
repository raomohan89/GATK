package org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.hellbender.GATKBaseTest;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.spark.sv.StructuralVariationDiscoveryArgumentCollection;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SvCigarUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;
import scala.Tuple2;
import scala.Tuple3;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;

import static org.broadinstitute.hellbender.tools.spark.sv.discovery.alignment.AlignmentInterval.NO_NM;

public class ContigAlignmentsModifierUnitTest extends GATKBaseTest {

    @Test(groups = "sv")
    public void testCompactifyNeighboringSoftClippings() {
        Assert.assertEquals(new Cigar(SvCigarUtils.compactifyNeighboringSoftClippings(TextCigarCodec.decode("1H2S3S4M5D6M7I8M9S10S11H")
                        .getCigarElements())),
                TextCigarCodec.decode("1H5S4M5D6M7I8M19S11H"));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_OneInsertion() {

        final Cigar cigar = TextCigarCodec.decode("56S27M15I32M21S");
        final AlignmentInterval alignmentInterval = new AlignmentInterval(new SimpleInterval("1", 100, 158),
                57, 130, cigar, true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        final List<AlignmentInterval> generatedARList = Utils.stream(ContigAlignmentsModifier.splitGappedAlignment(alignmentInterval,
                1, cigar.getReadLength())).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentInterval(new SimpleInterval("1", 100, 126),
                57, 83, TextCigarCodec.decode("56S27M68S"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(1), new AlignmentInterval(new SimpleInterval("1", 127, 158),
                99, 130, TextCigarCodec.decode("98S32M21S"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));

    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_OneDeletion() {
        final Cigar cigar = TextCigarCodec.decode("2S205M2D269M77S");
        final AlignmentInterval alignmentInterval = new AlignmentInterval(new SimpleInterval("1", 100, 575),
                208, 476, cigar, true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        final List<AlignmentInterval> generatedARList = Utils.stream(ContigAlignmentsModifier.splitGappedAlignment(alignmentInterval,
                1, cigar.getReadLength())).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentInterval(new SimpleInterval("1", 100, 304),
                3, 207, TextCigarCodec.decode("2S205M346S"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(1), new AlignmentInterval(new SimpleInterval("1", 307, 575),
                208, 476, TextCigarCodec.decode("207S269M77S"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_Complex() {

        final Cigar cigar = TextCigarCodec.decode("397S118M2D26M6I50M7I26M1I8M13D72M398S");
        final AlignmentInterval alignmentInterval = new AlignmentInterval(new SimpleInterval("1", 100, 414),
                398, 711, cigar, true, 60, 65, 100, ContigAlignmentsModifier.AlnModType.NONE);

        final List<AlignmentInterval> generatedARList = Utils.stream(ContigAlignmentsModifier.splitGappedAlignment(alignmentInterval,
                1, cigar.getReadLength())).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 6);

        Assert.assertEquals(generatedARList.get(0), new AlignmentInterval(new SimpleInterval("1", 100, 217),
                398, 515, TextCigarCodec.decode("397S118M594S"),
                true, 60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(1), new AlignmentInterval(new SimpleInterval("1", 220, 245),
                516, 541, TextCigarCodec.decode("515S26M568S"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(2), new AlignmentInterval(new SimpleInterval("1", 246, 295),
                548, 597, TextCigarCodec.decode("547S50M512S"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(3), new AlignmentInterval(new SimpleInterval("1", 296, 321),
                605, 630, TextCigarCodec.decode("604S26M479S"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(4), new AlignmentInterval(new SimpleInterval("1", 322, 329),
                632, 639, TextCigarCodec.decode("631S8M470S"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(5), new AlignmentInterval(new SimpleInterval("1", 343, 414),
                640, 711, TextCigarCodec.decode("639S72M398S"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_GapSizeSensitivity() {

        final Cigar cigar = TextCigarCodec.decode("10M10D10M60I10M10I10M50D10M");
        final AlignmentInterval alignmentInterval = new AlignmentInterval(new SimpleInterval("1", 100, 209),
                1, 120, cigar, true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        final List<AlignmentInterval> generatedARList = Utils.stream(ContigAlignmentsModifier.splitGappedAlignment(alignmentInterval,
                StructuralVariationDiscoveryArgumentCollection.DiscoverVariantsFromContigsAlignmentsSparkArgumentCollection.GAPPED_ALIGNMENT_BREAK_DEFAULT_SENSITIVITY,
                cigar.getReadLength())).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 3);
        Assert.assertEquals(generatedARList.get(0), new AlignmentInterval(new SimpleInterval("1", 100, 129),
                1, 20, TextCigarCodec.decode("10M10D10M100S"),
                true, 60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(1), new AlignmentInterval(new SimpleInterval("1", 130, 149),
                81, 110, TextCigarCodec.decode("80S10M10I10M10S"),
                true, 60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(2), new AlignmentInterval(new SimpleInterval("1", 200, 209),
                111, 120, TextCigarCodec.decode("110S10M"),
                true, 60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_HardAndSoftClip() {

        final Cigar cigar = TextCigarCodec.decode("1H2S3M5I10M20D6M7S8H");
        final AlignmentInterval alignmentInterval = new AlignmentInterval(new SimpleInterval("1", 100, 138),
                4, 27, cigar, true, 60, 0,
                100, ContigAlignmentsModifier.AlnModType.NONE);

        final List<AlignmentInterval> generatedARList = Utils.stream(ContigAlignmentsModifier.splitGappedAlignment(alignmentInterval,
                1, cigar.getReadLength()+1+8)).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.get(0), new AlignmentInterval(new SimpleInterval("1", 100, 102),
                4, 6, TextCigarCodec.decode("1H2S3M28S8H"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(1), new AlignmentInterval(new SimpleInterval("1", 103, 112),
                12, 21, TextCigarCodec.decode("1H10S10M13S8H"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(2), new AlignmentInterval(new SimpleInterval("1", 133, 138),
                22, 27, TextCigarCodec.decode("1H20S6M7S8H"), true,
                60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_TerminalInsertionOperatorToSoftClip() {

        // beginning with 'I'
        Cigar cigar = TextCigarCodec.decode("10I10M5I10M");
        AlignmentInterval alignmentInterval = new AlignmentInterval(new SimpleInterval("1", 101, 120),
                11, 35, cigar, true, 60, 0,
                100, ContigAlignmentsModifier.AlnModType.NONE);

        List<AlignmentInterval> generatedARList = Utils.stream(ContigAlignmentsModifier.splitGappedAlignment(alignmentInterval, 1,
                cigar.getReadLength())).collect(Collectors.toList());

        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentInterval(
                new SimpleInterval("1", 101, 110), 11, 20,
                TextCigarCodec.decode("10S10M15S"),
                true, 60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(1), new AlignmentInterval(
                new SimpleInterval("1", 111, 120), 26, 35,
                TextCigarCodec.decode("25S10M"),
                true, 60, NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));

        // ending with 'I'
        cigar = TextCigarCodec.decode("10M5I10M10I");
        alignmentInterval = new AlignmentInterval(
                new SimpleInterval("1", 101, 120), 1, 25, cigar,
                true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);

        generatedARList = Utils.stream(ContigAlignmentsModifier.splitGappedAlignment(alignmentInterval, 1, cigar.getReadLength())).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);
        Assert.assertEquals(generatedARList.get(0), new AlignmentInterval(
                new SimpleInterval("1", 101, 110), 1, 10,
                TextCigarCodec.decode("10M25S"), true, 60,
                NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(1), new AlignmentInterval(
                new SimpleInterval("1", 111, 120), 16, 25,
                TextCigarCodec.decode("15S10M10S"), true, 60,
                NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_TerminalInsertionNeighboringClippings(){

        Cigar cigar = TextCigarCodec.decode("10H20S30I40M50I60S70H");
        AlignmentInterval alignmentInterval = new AlignmentInterval(
                new SimpleInterval("1", 101, 140), 61, 100, cigar,
                true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        List<AlignmentInterval> generatedARList = Utils.stream(ContigAlignmentsModifier.splitGappedAlignment(alignmentInterval, 1, cigar.getReadLength()+10+70))
                .collect(Collectors.toList());
        // no internal gap, so nothing should change
        Assert.assertEquals(generatedARList.size(), 1);
        Assert.assertEquals(generatedARList.get(0), alignmentInterval);

        cigar = TextCigarCodec.decode("10H20S30I40M5D15M50I60S70H");
        alignmentInterval = new AlignmentInterval(
                new SimpleInterval("1", 101, 160), 61, 115, cigar,
                true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        generatedARList = Utils.stream(ContigAlignmentsModifier.splitGappedAlignment(alignmentInterval, 1, cigar.getReadLength()+10+70)).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.size(), 2);

        Assert.assertEquals(generatedARList.get(0), new AlignmentInterval(
                new SimpleInterval("1", 101, 140), 61, 100,
                TextCigarCodec.decode("10H50S40M125S70H"), true, 60,
                NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(1), new AlignmentInterval(
                new SimpleInterval("1", 146, 160), 101, 115,
                TextCigarCodec.decode("10H90S15M110S70H"), true, 60,
                NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_NegativeStrand() {
        // read data with AlignedAssembly.AlignmentInterval.toString():
        // 19149	contig-8	chrUn_JTFH01000557v1_decoy	21	1459	-	10S1044M122I395M75I	60	11	1646	200
        final Cigar cigar = TextCigarCodec.decode("10S1044M122I395M75I");
        final AlignmentInterval alignmentInterval = new AlignmentInterval(
                new SimpleInterval("chrUn_JTFH01000557v1_decoy", 21, 1459), 11,
                1646, cigar, false, 60, 200, 100, ContigAlignmentsModifier.AlnModType.NONE);
        final List<AlignmentInterval> generatedARList = Utils.stream(ContigAlignmentsModifier.splitGappedAlignment(alignmentInterval,
                1, cigar.getReadLength())).collect(Collectors.toList());
        Assert.assertEquals(generatedARList.get(0), new AlignmentInterval(
                new SimpleInterval("chrUn_JTFH01000557v1_decoy", 416, 1459), 11,
                1054, TextCigarCodec.decode("10S1044M592S"), false, 60,
                NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
        Assert.assertEquals(generatedARList.get(1), new AlignmentInterval(
                new SimpleInterval("chrUn_JTFH01000557v1_decoy", 21, 415), 1177,
                1571, TextCigarCodec.decode("1176S395M75S"), false, 60,
                NO_NM, 100, ContigAlignmentsModifier.AlnModType.FROM_SPLIT_GAPPED_ALIGNMENT));
    }

    @Test(groups = "sv")
    public void testGappedAlignmentBreaker_ExpectException() {
        int exceptionCount = 0;

        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10H10D10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10S10D10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10D10S"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10D10H"));} catch (final Exception ex) {++exceptionCount;}

        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10I10S"));} catch (final Exception ex) {++exceptionCount;}
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10D10S"));} catch (final Exception ex) {++exceptionCount;}

        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10I10M10D10S"));} catch (final Exception ex) {++exceptionCount;}
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10D10M10I10S"));} catch (final Exception ex) {++exceptionCount;}

        // these 4 are fine now
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10H10I10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10S10I10M"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10I10S"));} catch (final Exception ex) {++exceptionCount;}
        try {willThrowOnInvalidCigar(TextCigarCodec.decode("10M10I10H"));} catch (final Exception ex) {++exceptionCount;}

        // last two are valid
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10M10I10M10S"));} catch (final Exception ex) {++exceptionCount;}
        try{willThrowOnInvalidCigar(TextCigarCodec.decode("10H10M10D10M10S"));} catch (final Exception ex) {++exceptionCount;}

        Assert.assertEquals(exceptionCount, 8);
    }

    private static Iterable<AlignmentInterval> willThrowOnInvalidCigar(final Cigar cigar) throws GATKException {
        final AlignmentInterval detailsDoesnotMatter = new AlignmentInterval(
                new SimpleInterval("1", 1, 110), 21, 30, cigar,
                true, 60, 0, 100, ContigAlignmentsModifier.AlnModType.NONE);
        return ContigAlignmentsModifier.splitGappedAlignment(detailsDoesnotMatter, 1, cigar.getReadLength() + SvCigarUtils.getTotalHardClipping(cigar));
    }

    //==================================================================================================================

    @DataProvider(name = "forComputeNewRefSpanAndCigar")
    private Object[][] createTestDataForComputeNewRefSpanAndCigar() {

        final List<Object[]> data = new ArrayList<>(20);

        AlignmentInterval alignment = new AlignmentInterval(new SimpleInterval("chr1", 175417007, 175417074),
                14, 81, TextCigarCodec.decode("13H68M394H"),
                true, 60, 0, 68, ContigAlignmentsModifier.AlnModType.NONE);
        SimpleInterval refSpan = new SimpleInterval("chr1", 175417007, 175417043);
        data.add(new Object[]{alignment, 31, true, refSpan, TextCigarCodec.decode("13H37M31S394H")});

        alignment = new AlignmentInterval(new SimpleInterval("chr2", 122612655, 122612751),
                9, 105, TextCigarCodec.decode("8H97M138H"),
                false, 60, 0, 97, ContigAlignmentsModifier.AlnModType.NONE);
        refSpan = new SimpleInterval("chr2", 122612659, 122612751);
        data.add(new Object[]{alignment, 4, true, refSpan, TextCigarCodec.decode("8H93M4S138H")});

        alignment = new AlignmentInterval(new SimpleInterval("chr6", 66782514, 66782679),
                32, 197, TextCigarCodec.decode("31S166M"),
                false, 60, 3, 151, ContigAlignmentsModifier.AlnModType.NONE);
        refSpan = new SimpleInterval("chr6", 66782514, 66782675);
        data.add(new Object[]{alignment, 4, false, refSpan, TextCigarCodec.decode("35S162M")});

        alignment = new AlignmentInterval(new SimpleInterval("chr2", 91421528, 91421734),
                271, 477, TextCigarCodec.decode("270H207M"),
                true, 40, 12, 147, ContigAlignmentsModifier.AlnModType.NONE);
        refSpan = new SimpleInterval("chr2", 91421560, 91421734);
        data.add(new Object[]{alignment, 32, false, refSpan, TextCigarCodec.decode("270H32S175M")});

        final SimpleInterval originalRefSpan = new SimpleInterval("chr2", 128791173, 128792506);
        alignment = new AlignmentInterval(originalRefSpan,
                1, 1332, TextCigarCodec.decode("1190M4D53M2I26M2I31M2D28M1422S"),
                true, 60, 13, 1239, ContigAlignmentsModifier.AlnModType.NONE);
        refSpan = new SimpleInterval("chr2", 128791173, 128792476);
        data.add(new Object[]{alignment, 28, true, refSpan, TextCigarCodec.decode("1190M4D53M2I26M2I31M1450S")});

        alignment = new AlignmentInterval(originalRefSpan,
                1, 1334, TextCigarCodec.decode("1190M4D53M2I26M2I31M2I28M1422S"),
                true, 60, 13, 1239, ContigAlignmentsModifier.AlnModType.NONE);
        refSpan = new SimpleInterval("chr2", 128791173, 128792476);
        data.add(new Object[]{alignment, 28, true, refSpan, TextCigarCodec.decode("1190M4D53M2I26M2I31M1452S")});

        alignment = new AlignmentInterval(originalRefSpan,
                1, 1332, TextCigarCodec.decode("1190M4D53M2I26M2I31M2D28M1422S"),
                true, 60, 13, 1239, ContigAlignmentsModifier.AlnModType.NONE);
        refSpan = new SimpleInterval("chr2", 128792367, 128792506);
        data.add(new Object[]{alignment, 1190, false, refSpan, TextCigarCodec.decode("1190S53M2I26M2I31M2D28M1422S")});

        alignment = new AlignmentInterval(originalRefSpan,
                1, 1334, TextCigarCodec.decode("1190M4I53M2I26M2I31M2I28M1422S"),
                true, 60, 13, 1239, ContigAlignmentsModifier.AlnModType.NONE);
        refSpan = new SimpleInterval("chr2", 128792367, 128792506);
        data.add(new Object[]{alignment, 1190, false, refSpan, TextCigarCodec.decode("1194S53M2I26M2I31M2I28M1422S")});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "forComputeNewRefSpanAndCigar", groups = "sv")
    public void testComputeNewRefSpanAndCigar(final AlignmentInterval interval, final int clipLength, final boolean clipFrom3PrimeEnd,
                                              final SimpleInterval expectedRefSpan, final Cigar expectedCigar) {

        final Tuple2<SimpleInterval, Cigar> x = ContigAlignmentsModifier.computeNewRefSpanAndCigar(interval, clipLength, clipFrom3PrimeEnd);
        Assert.assertEquals(x._1, expectedRefSpan);
        Assert.assertEquals(x._2, expectedCigar);
    }

    @DataProvider(name = "forCigarExtraction")
    private Object[][] createTestDataForCigarExtraction() {

        final List<Object[]> data = new ArrayList<>(20);

        SimpleInterval refSpan = new SimpleInterval("chr1", 82666357, 82666765);
        AlignmentInterval alignment = new AlignmentInterval(refSpan, 69, 472,
                TextCigarCodec.decode("68S122M5D282M"), true, 60, 11,
                353, ContigAlignmentsModifier.AlnModType.NONE);
        List<CigarElement> left = Arrays.asList(new CigarElement(68, CigarOperator.S));
        List<CigarElement> middle = Arrays.asList(new CigarElement(122, CigarOperator.M), new CigarElement(5, CigarOperator.D), new CigarElement(282, CigarOperator.M));
        List<CigarElement> right = Collections.emptyList();
        data.add(new Object[]{alignment, left, middle, right});

        refSpan = new SimpleInterval("chr3", 61792401, 61792448);
        alignment = new AlignmentInterval(refSpan, 43, 90,
                TextCigarCodec.decode("42H48M382H"), true, 46, 1,
                43, ContigAlignmentsModifier.AlnModType.NONE);
        left = Arrays.asList(new CigarElement(42, CigarOperator.H));
        middle = Arrays.asList(new CigarElement(48, CigarOperator.M));
        right = Arrays.asList(new CigarElement(382, CigarOperator.H));
        data.add(new Object[]{alignment, left, middle, right});

        refSpan = new SimpleInterval("chrY", 26303624, 26303671);
        alignment = new AlignmentInterval(refSpan, 1, 48,
                TextCigarCodec.decode("48M424H"), true, 0, 2,
                38, ContigAlignmentsModifier.AlnModType.NONE);
        left = Collections.emptyList();
        middle = Arrays.asList(new CigarElement(48, CigarOperator.M));
        right = Arrays.asList(new CigarElement(424, CigarOperator.H));
        data.add(new Object[]{alignment, left, middle, right});

        return data.toArray(new Object[data.size()][]);
    }

    @Test(dataProvider = "forCigarExtraction", groups = "sv")
    public void testExtractCigar(final AlignmentInterval interval, final List<CigarElement> expectedLeft,
                                 final List<CigarElement> expectedMiddle, final List<CigarElement> expectedRight) {

        final Tuple3<List<CigarElement>, List<CigarElement>, List<CigarElement>> x =
                ContigAlignmentsModifier.splitCigarByLeftAndRightClipping(interval);
        Assert.assertEquals(x._1(), expectedLeft);
        Assert.assertEquals(x._2(), expectedMiddle);
        Assert.assertEquals(x._3(), expectedRight);
    }
}
