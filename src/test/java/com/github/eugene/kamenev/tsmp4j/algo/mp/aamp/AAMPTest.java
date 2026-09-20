package com.github.eugene.kamenev.tsmp4j.algo.mp.aamp;

import com.github.eugene.kamenev.tsmp4j.BaseTest;
import com.github.eugene.kamenev.tsmp4j.stats.RollingWindowWithoutStatistics;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;

class AAMPTest extends BaseTest {

    @Test
    void testAampSelfJoin() {
        var limit = 200;
        var windowSize = 30;
        var check = loadMP("aamp_self_join.csv", AAMPTest.class);

        var aamp = new AAMP(new RollingWindowWithoutStatistics(windowSize, limit), 2.0);

        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(aamp::update);

        var mp = aamp.get();

        equals(check.mp(), mp.profile());
        equals(check.pi(), mp.indexes());
    }

    @Test
    void testAampAbJoin() {
        var limit = 200;
        var windowSize = 30;
        var queryLimit = 60;
        var check = loadMP("aamp_ab_join.csv", AAMPTest.class);
        var checkQuery = loadMP("aamp_ba_join.csv", AAMPTest.class);

        var query = new RollingWindowWithoutStatistics(windowSize, queryLimit);
        data.stream()
            .mapToDouble(t -> t.y())
            .limit(queryLimit)
            .forEach(query::apply);

        var aamp = new AAMP(new RollingWindowWithoutStatistics(windowSize, limit), 2.0);

        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(aamp::update);

        var mp = aamp.get(query);

        equals(check.mp(), mp.profile());
        equals(check.pi(), mp.indexes());
        equals(checkQuery.mp(), mp.leftProfile());
        equals(checkQuery.pi(), mp.leftIndexes());
    }

    @Test
    void testAampAbJoinManhattanDistance() {
        var limit = 200;
        var windowSize = 30;
        var queryLimit = 60;
        var check = loadMP("aamp_ab_join_p1.csv", AAMPTest.class);
        var checkQuery = loadMP("aamp_ba_join_p1.csv", AAMPTest.class);

        var query = new RollingWindowWithoutStatistics(windowSize, queryLimit);
        data.stream()
            .mapToDouble(t -> t.y())
            .limit(queryLimit)
            .forEach(query::apply);

        var aamp = new AAMP(new RollingWindowWithoutStatistics(windowSize, limit), 1.0);

        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(aamp::update);

        var mp = aamp.get(query);

        equals(check.mp(), mp.profile());
        equals(check.pi(), mp.indexes());
        equals(checkQuery.mp(), mp.leftProfile());
        equals(checkQuery.pi(), mp.leftIndexes());
    }

    @Test
    void testAampAbJoinWithQueryAsArray() {
        var limit = 200;
        var windowSize = 30;
        var queryLimit = 60;
        var check = loadMP("aamp_ab_join.csv", AAMPTest.class);

        var query = data.stream()
            .mapToDouble(t -> t.y())
            .limit(queryLimit)
            .toArray();

        var aamp = new AAMP(new RollingWindowWithoutStatistics(windowSize, limit), 2.0);

        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(aamp::update);

        var mp = aamp.get(query);

        equals(check.mp(), mp.profile());
        equals(check.pi(), mp.indexes());
    }

    @Test
    void testAampAbJoinMatchesNaiveComputation() {
        var windowSize = 10;
        var limit = 50;
        var queryLimit = 25;

        var ts = data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .toArray();
        var query = data.stream()
            .mapToDouble(t -> t.y())
            .limit(queryLimit)
            .toArray();

        var aamp = new AAMP(RollingWindowWithoutStatistics.of(ts, windowSize), 2.0);
        var mp = aamp.get(query);

        var check = naiveProfile(ts, query, windowSize);
        equals(check.mp(), mp.profile());
        equals(check.pi(), mp.indexes());
    }

    @Test
    void testAampAbJoinNotReadyReturnsNull() {
        var windowSize = 30;

        var aamp = new AAMP(new RollingWindowWithoutStatistics(windowSize, 100), 2.0);
        data.stream()
            .mapToDouble(t -> t.x())
            .limit(50)
            .forEach(aamp::update);

        assertNull(aamp.get(new RollingWindowWithoutStatistics(windowSize, 20)));
    }

    @Test
    void testAampAbJoinQueryShorterThanWindow() {
        var windowSize = 30;

        var aamp = new AAMP(new RollingWindowWithoutStatistics(windowSize, 200), 2.0);
        data.stream()
            .mapToDouble(t -> t.x())
            .limit(200)
            .forEach(aamp::update);

        assertThrows(IllegalArgumentException.class,
            () -> aamp.get(RollingWindowWithoutStatistics.of(new double[]{1, 2, 3}, windowSize)));
    }

    /**
     * Reference computation of the AB join: distance between every pair of subsequences,
     * without any optimization.
     */
    private static NaiveProfile naiveProfile(double[] ts, double[] query, int windowSize) {
        var n = ts.length - windowSize + 1;
        var qn = query.length - windowSize + 1;
        var mp = new double[n];
        var pi = new int[n];
        for (int i = 0; i < n; i++) {
            var min = Double.POSITIVE_INFINITY;
            var minIndex = -1;
            for (int j = 0; j < qn; j++) {
                var d = 0d;
                for (int k = 0; k < windowSize; k++) {
                    d += Math.pow(ts[i + k] - query[j + k], 2);
                }
                if (d < min) {
                    min = d;
                    minIndex = j;
                }
            }
            mp[i] = Math.sqrt(min);
            pi[i] = minIndex;
        }
        return new NaiveProfile(mp, pi);
    }

    private record NaiveProfile(double[] mp, int[] pi) { }
}
