package com.github.eugene.kamenev.tsmp4j.algo.fluss;

import com.github.eugene.kamenev.tsmp4j.BaseTest;
import com.github.eugene.kamenev.tsmp4j.algo.mp.stomp.STOMP;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

class FLUSSTest extends BaseTest {

    private static final double[] FLUSS_DATA = loadData("fluss_data_tilt_abp.csv", rows -> Arrays.stream(rows)
        .mapToDouble(s -> Double.parseDouble(s[0]))
        .toArray(), FLUSSTest.class);

    private static final double[] FLUSS_CAC_CHECK = loadData("fluss_cac.csv", rows -> Arrays.stream(rows)
        .mapToDouble(s -> Double.parseDouble(s[0]))
        .toArray(), FLUSSTest.class);

    @Test
    void testFluss() {
        int windowSize = 10;
        int buffSize = 1000;
        int numSegments = 2;

        var stomp = new STOMP(windowSize, buffSize);

        Arrays.stream(FLUSS_DATA)
            .limit(buffSize)
            .forEach(stomp::update);

        var mp = stomp.get();

        var fluss = new FLUSS(windowSize, numSegments).apply(mp);

        equals(fluss.cac(), FLUSS_CAC_CHECK);
        equals(fluss.changePoints(), new int[]{940, 874});
    }
}
