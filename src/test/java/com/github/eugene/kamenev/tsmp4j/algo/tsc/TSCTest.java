package com.github.eugene.kamenev.tsmp4j.algo.tsc;

import com.github.eugene.kamenev.tsmp4j.BaseTest;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

import static org.junit.jupiter.api.Assertions.assertEquals;

class TSCTest extends BaseTest {

    private static final double[] TS = load1D("data.csv", TSCTest.class);
    private static final int[] VALID_CHAIN = {5122, 5318, 5518, 5917, 6312, 6507, 6700, 7457, 8034, 8613, 9396};

    @Test
    void testRobustTimeSeriesChains() {
        var tsc = new TSC(180, TS.length, 0.25d);
        Arrays.stream(TS)
            .forEach(tsc::update);
        var chains = tsc.get();
        var score = TSC.bestScore(TS, chains, 180);

        assertEquals(9396, score.index());
        assertEquals(4.0d, score.score(), 0d);
        equals(score.predIdxOurs(), VALID_CHAIN);
    }
}
