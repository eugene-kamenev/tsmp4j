package com.github.eugene.kamenev.tsmp4j.algo.mp.mpx;

import com.github.eugene.kamenev.tsmp4j.BaseTest;
import com.github.eugene.kamenev.tsmp4j.algo.mp.DistanceProfileFunction;
import org.junit.jupiter.api.Test;

import java.util.stream.Stream;

import static org.junit.jupiter.api.Assertions.assertEquals;

class MPXTest extends BaseTest {

    public static final double ERROR = Math.pow(10, -12);

    @Override
    protected double error() {
        return ERROR;
    }

    @Test
    void testMovingStatsMethod() {
        var windowSize = 30;
        var limit = 200;
        var statsCheck = loadData("stats.csv", rows -> {
            double[] mean = new double[rows.length];
            double[] std = new double[rows.length];
            for (int i = 0; i < rows.length; i++) {
                mean[i] = Double.parseDouble(rows[i][0]);
                std[i] = Double.parseDouble(rows[i][1]);
            }
            return new Stats(mean, std);
        }, MPXTest.class);

        var stats = new MPXRollingWindowStatistics(windowSize, 1);
        var ts = data.stream()
            .limit(limit)
            .mapToDouble(t -> t.x())
            .toArray();
        var result = new double[2][ts.length - windowSize + 1];
        for (int i = 0, k = 0; i < ts.length; i++) {
            var s = stats.apply(ts[i]);
            if (stats.isReady()) {
                result[0][k] = s.mean();
                result[1][k++] = s.stdDev();
            }
        }

        equals(statsCheck.mean(), result[0]);
        equals(statsCheck.std(), result[1]);
    }

    @Test
    void testMpxNoQueryCrossCorrelationFalse() {
        var windowSize = 30;
        var limit = 200;
        var checkMp = loadMP("mpx_toy_euclidean.csv", MPXTest.class);

        var mpx = new MPX(windowSize, limit, false);
        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(mpx::update);
        var mp = mpx.get();

        equals(checkMp.mp(), mp.profile());
        equals(checkMp.pi(), mp.indexes());
    }

    @Test
    void testMpxNoQueryCrossCorrelationTrue() {
        var windowSize = 30;
        var limit = 200;
        var checkMp = loadMP("mpx_toy_pearson.csv", MPXTest.class);

        var mpx = new MPX(windowSize, limit, true);
        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(mpx::update);
        var mp = mpx.get();

        equals(checkMp.mp(), mp.profile());
        equals(checkMp.pi(), mp.indexes());
    }

    @Test
    void testMpxWithQueryCrossCorrelationFalse() {
        var limit = 200;
        var windowSize = 30;
        var check = MPQuery.load("mpx_toy_query_euclidean.csv");
        var query = data.stream()
            .limit(limit)
            .mapToDouble(c -> c.y())
            .toArray();

        var algo = new MPX(windowSize, limit, false);
        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(algo::update);
        var mpx = algo.get(query);

        equals(mpx.profile(), check.mp());
        equals(mpx.indexes(), check.pi());
        equals(mpx.leftIndexes(), check.pib());
        equals(mpx.leftProfile(), check.mpb());
    }

    @Test
    void testMpxWithQueryCrossCorrelationTrue() {
        var limit = 200;
        var windowSize = 30;
        var check = MPQuery.load("mpx_toy_query_pearson.csv");
        var query = data.stream()
            .limit(limit)
            .mapToDouble(c -> c.y())
            .toArray();

        var algo = new MPX(windowSize, limit, true);
        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(algo::update);
        var mpx = algo.get(query);

        equals(mpx.profile(), check.mp());
        equals(mpx.indexes(), check.pi());
        equals(mpx.leftIndexes(), check.pib());
        equals(mpx.leftProfile(), check.mpb());
    }

    @Test
    void testMpxStreamingProducesSameMatrixProfile() {
        var windowSize = 30;
        var limit = 200;

        var mpx = new MPX(windowSize, limit, false);
        var mpx2 = new MPX(windowSize, limit, false);
        var ystream = data.stream()
            .mapToDouble(t -> t.y())
            .limit(10); // additional 10 points at the beginning
        var xstream = data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit);

        Stream.concat(ystream.boxed(), xstream.boxed())
            .forEach(mpx2::update);
        data.stream()
            .limit(limit)
            .mapToDouble(t -> t.x())
            .forEach(mpx::update);
        var mp = mpx.get();
        var mp2 = mpx2.get();

        equals(mp2.profile(), mp.profile());
        equals(mp2.indexes(), mp.indexes());
    }

    @Test
    void testMpdistSameSize() {
        var windowSize = 30;
        var ts = new MPXRollingWindowStatistics(windowSize, data.size());
        var qts = new MPXRollingWindowStatistics(windowSize, data.size());

        data.stream()
            .mapToDouble(t -> t.x())
            .forEach(ts::apply);
        data.stream()
            .mapToDouble(t -> t.y())
            .forEach(qts::apply);

        var query = new DistanceProfileFunction.DistanceProfileQuery<>(ts, qts, windowSize);
        var dist = new MPX().apply(query);

        assertEquals(MPDist.load("mpdist_same_size.csv", MPXTest.class).x()[0], dist.profile()[0], ERROR);
    }

    @Test
    void testMpdistDiffSize() {
        var windowSize = 30;
        var skip = 150;
        var limit = 50;
        var ts = new MPXRollingWindowStatistics(windowSize, data.size());
        var qts = new MPXRollingWindowStatistics(windowSize, limit);

        data.stream()
            .mapToDouble(t -> t.y())
            .forEach(ts::apply);
        data.stream()
            .skip(skip)
            .limit(limit)
            .mapToDouble(t -> t.x())
            .forEach(qts::apply);

        var query = new DistanceProfileFunction.DistanceProfileQuery<>(ts, qts, windowSize);
        var dist = new MPX().apply(query);

        assertEquals(MPDist.load("mpdist_diff_size.csv", MPXTest.class).x()[0], dist.profile()[0], ERROR);
    }

    public record MPQuery(double[] mp, int[] pi, double[] mpb, int[] pib) {

        public static MPQuery load(String file) {
            return loadData(file, rows -> {
                double[] mp = new double[rows.length];
                int[] pi = new int[rows.length];
                double[] mpb = new double[rows.length];
                int[] mpib = new int[rows.length];
                for (int i = 0; i < rows.length; i++) {
                    mp[i] = Double.parseDouble(rows[i][0]);
                    pi[i] = Integer.parseInt(rows[i][1]) - 1; // in R indexing starts from 1
                    mpb[i] = Double.parseDouble(rows[i][2]);
                    mpib[i] = Integer.parseInt(rows[i][3]) - 1; // in R indexing starts from 1
                }
                return new MPQuery(mp, pi, mpb, mpib);
            }, MPXTest.class);
        }
    }

    public record Stats(double[] mean, double[] std) { }
}
