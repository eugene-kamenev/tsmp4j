package com.github.eugene.kamenev.tsmp4j.algo.cp;

import java.util.Arrays;
import java.util.Random;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

class PlatoNearestNeighborTest {

    @Test
    void rollingScanMatchesNaiveReference() {
        for (int m : new int[]{12, 20}) {
            for (int length : new int[]{60, 101, 150, 240, 400}) {
                for (double offset : new double[]{0.0d, 1e8d}) {
                    var query = series(length, 7L + length, offset, 5.0d);
                    var pattern = Arrays.copyOfRange(query, 37, 37 + m);
                    assertEquals(naiveMinDistance(query, pattern),
                        PlatoNearestNeighbor.minDistance(query, pattern), 1e-9d,
                        "m=" + m + " length=" + length + " offset=" + offset);
                }
            }
        }
    }

    @Test
    void detectsSubsequenceCarriedOnLargeOffset() {
        var pattern = new double[20];
        for (int i = 0; i < pattern.length; i++) {
            pattern[i] = 3.0d * Math.sin(i / 3.0d);
        }
        var query = series(200, 11L, 1e6d, 0.1d);
        System.arraycopy(pattern, 0, query, 90, pattern.length);
        for (int i = 0; i < pattern.length; i++) {
            query[90 + i] += 1e6d;
        }

        assertEquals(0.0d, PlatoNearestNeighbor.minDistance(query, pattern), 1e-6d);
    }

    @Test
    void rejectsQueryShorterThanPlato() {
        assertThrows(IllegalArgumentException.class,
            () -> PlatoNearestNeighbor.minDistance(new double[]{1.0d, 2.0d, 3.0d}, new double[20]));
        assertThrows(IllegalArgumentException.class,
            () -> PlatoNearestNeighbor.classify(new double[]{1.0d, 2.0d},
                new double[][]{new double[20], new double[20]}));
    }

    @Test
    void emptyPatternMatchesTrivially() {
        assertEquals(0.0d,
            PlatoNearestNeighbor.minDistance(series(60, 3L, 0.0d, 5.0d), new double[0]), 0.0d);
        assertEquals(0.0d,
            PlatoNearestNeighbor.minDistance(new double[0], new double[0]), 0.0d);
    }

    @Test
    void keepsRunningStatisticsFiniteAtExtremeMagnitudes() {
        // Every window sits on a value so large that (low + high) / 2 overflows to infinity: the
        // running residuals must be carried around a midpoint that stays representable.
        var extreme = Double.MAX_VALUE * 0.9d;
        var constantQuery = new double[60];
        Arrays.fill(constantQuery, extreme);
        var constantPlato = new double[20];
        Arrays.fill(constantPlato, extreme);

        assertEquals(0.0d,
            PlatoNearestNeighbor.minDistance(constantQuery, constantPlato), 1e-9d);

        var alternating = new double[20];
        for (int i = 0; i < alternating.length; i++) {
            alternating[i] = i % 2 == 0 ? 1.0d : -1.0d;
        }
        var distance = PlatoNearestNeighbor.minDistance(constantQuery, alternating);
        assertTrue(Double.isFinite(distance));
        assertEquals(Math.sqrt(alternating.length), distance, 1e-9d);
    }

    @Test
    void treatsNearConstantWindowAsConstant() {
        var pattern = new double[20];
        for (int i = 0; i < pattern.length; i++) {
            pattern[i] = i % 2 == 0 ? 1.0d : -1.0d;
        }
        // Jitter far below any real signal: such a window carries no shape information and must
        // not be amplified into arbitrary normalized values.
        var query = new double[60];
        Arrays.fill(query, 1.0d);
        var random = new Random(5L);
        for (int i = 0; i < query.length; i++) {
            query[i] += (random.nextDouble() * 2 - 1) * 1e-13d;
        }

        double distance = PlatoNearestNeighbor.minDistance(query, pattern);
        assertTrue(Double.isFinite(distance));
        assertEquals(Math.sqrt(pattern.length), distance, 1e-9d);
    }

    private static double[] series(int length, long seed, double offset, double amplitude) {
        var random = new Random(seed);
        var ts = new double[length];
        for (int i = 0; i < length; i++) {
            ts[i] = offset + amplitude * Math.sin(i / 7.0d)
                + (random.nextDouble() * 2 - 1) * 0.1d * amplitude;
        }
        return ts;
    }

    /**
     * Reference computation: z-normalize every window independently, without running statistics.
     */
    private static double naiveMinDistance(double[] query, double[] pattern) {
        int m = pattern.length;
        var normalizedPattern = zNormalize(pattern);
        var min = Double.POSITIVE_INFINITY;
        for (int i = 0; i + m <= query.length; i++) {
            var normalizedWindow = zNormalize(Arrays.copyOfRange(query, i, i + m));
            var sum = 0.0d;
            for (int k = 0; k < m; k++) {
                var diff = normalizedWindow[k] - normalizedPattern[k];
                sum += diff * diff;
            }
            min = Math.min(min, Math.sqrt(sum));
        }
        return min;
    }

    private static double[] zNormalize(double[] data) {
        var mean = Arrays.stream(data).average().orElse(0);
        var variance = Arrays.stream(data).map(x -> (x - mean) * (x - mean)).average().orElse(0);
        var std = Math.sqrt(variance);
        if (std == 0) {
            return new double[data.length];
        }
        return Arrays.stream(data).map(x -> (x - mean) / std).toArray();
    }
}
