package com.github.eugene.kamenev.tsmp4j.algo.mp.mpx

import com.github.eugene.kamenev.tsmp4j.BaseSpec
import com.github.eugene.kamenev.tsmp4j.algo.mp.DistanceProfileFunction

import java.util.stream.Stream

import static org.hamcrest.Matchers.closeTo


class MPXSpec extends BaseSpec {

    public static double ERROR = Math.pow(10, -13);

    def 'test moving stats method'() {
        given:
        var windowSize = 30
        var limit = 200
        var statsCheck = loadData('stats.csv', (rows) -> {
            double[] mean = new double[rows.length]
            double[] std = new double[rows.length]
            int counter = 0
            Arrays.stream(rows).forEach(row -> {
                mean[counter] = row[0] as double
                std[counter++] = row[1] as double
            })
            return new Stats(mean, std)
        }, MPXSpec)

        when:
        var stats = new MPXRollingWindowStatistics(windowSize, 1)
        var ts = data
                .stream()
                .limit(limit)
                .mapToDouble(t -> t.x())
                .toArray()
        var result = new double[2][ts.length - windowSize + 1]
        for (var i = 0, k = 0; i < ts.length; i++) {
            var s = stats.apply(ts[i])
            if (stats.isReady()) {
                result[0][k] = s.mean()
                result[1][k++] = s.stdDev()
            }
        }

        then:
        equals(statsCheck.mean(), result[0])
        equals(statsCheck.std(), result[1])
    }

    def 'test mpx no query cross correlation false'() {
        given:
        var windowSize = 30
        var limit = 200
        var checkMp = loadMP('mpx_toy_euclidean.csv', MPXSpec)

        when:
        var mpx = new MPX(windowSize, limit, false)
        data.stream()
                .mapToDouble(t -> t.x())
                .limit(limit)
                .forEach(mpx::update)
        var mp = mpx.get()

        then:
        equals(checkMp.mp(), mp.profile())
        equals(checkMp.pi(), mp.indexes())
    }

    def 'test mpx no query cross correlation true'() {
        given:
        var windowSize = 30
        var limit = 200
        var checkMp = loadMP('mpx_toy_pearson.csv', MPXSpec)

        when:
        var mpx = new MPX(windowSize, limit, true)
        data.stream()
                .mapToDouble(t -> t.x())
                .limit(limit)
                .forEach(mpx::update)
        var mp = mpx.get()

        then:
        equals(checkMp.mp(), mp.profile())
        equals(checkMp.pi(), mp.indexes())
    }

    def 'test mpx with query cross correlation false'() {
        given:
        var limit = 200
        var windowSize = 30
        var check = MPQuery.load('mpx_toy_query_euclidean.csv')
        var query = data.stream()
                .limit(limit)
                .mapToDouble(c -> c.y())
                .toArray()

        when:
        var algo = new MPX(windowSize, limit, false)
        data.stream()
                .mapToDouble(t -> t.x())
                .limit(limit)
                .forEach(algo::update)
        var mpx = algo.get(query)

        then:
        equals(mpx.profile(), check.mp())
        equals(mpx.indexes(), check.pi())
        equals(mpx.leftIndexes(), check.pib())
        equals(mpx.leftProfile(), check.mpb())
    }

    def 'test mpx with query cross correlation true'() {
        given:
        var limit = 200
        var windowSize = 30
        var check = MPQuery.load('mpx_toy_query_pearson.csv')
        var query = data.stream()
                .limit(limit)
                .mapToDouble(c -> c.y())
                .toArray()

        when:
        var algo = new MPX(windowSize, limit, true)
        data.stream()
                .mapToDouble(t -> t.x())
                .limit(limit)
                .forEach(algo::update)
        var mpx = algo.get(query)

        then:
        equals(mpx.profile(), check.mp())
        equals(mpx.indexes(), check.pi())
        equals(mpx.leftIndexes(), check.pib())
        equals(mpx.leftProfile(), check.mpb())

    }

    def 'test mpx streaming produces same matrix profile'() {
        given:
        var windowSize = 30
        var limit = 200

        when:
        var mpx = new MPX(windowSize, limit, false)
        var mpx2 = new MPX(windowSize, limit, false)
        var ystream = data.stream()
                .mapToDouble(t -> t.y())
                .limit(10) // additional 10 points at the beginning
        var xstream = data.stream()
                .mapToDouble(t -> t.x())
                .limit(limit)

        Stream.concat(ystream.boxed(), xstream.boxed())
                .forEach(mpx2::update)
        data.stream()
                .limit(limit)
                .mapToDouble(t -> t.x())
                .forEach(mpx::update)
        var mp = mpx.get()
        var mp2 = mpx2.get()

        then:
        equals(mp2.profile(), mp.profile())
        equals(mp2.indexes(), mp.indexes())

    }

    def 'test mpdist same size'() {
        given:
        var windowSize = 30
        var ts = new MPXRollingWindowStatistics(windowSize, data.size())
        var qts = new MPXRollingWindowStatistics(windowSize, data.size())

        when:
        data.stream()
                .mapToDouble(t -> t.x())
                .forEach(ts::apply)
        data.stream()
                .mapToDouble(t -> t.y())
                .forEach(qts::apply)

        var query = new DistanceProfileFunction.DistanceProfileQuery(ts, qts, windowSize)
        var dist = new MPX().apply(query)

        then:
        closeTo(dist.profile()[0], ERROR).matches(MPDist.load('mpdist_same_size.csv', MPXSpec).x()[0])
    }

    def 'test mpdist diff size'() {
        given:
        var windowSize = 30
        var skip = 150
        var limit = 50
        var ts = new MPXRollingWindowStatistics(windowSize, data.size())
        var qts = new MPXRollingWindowStatistics(windowSize, limit)

        when:
        data.stream()
                .mapToDouble(t -> t.y())
                .forEach(ts::apply)
        data.stream()
                .skip(skip)
                .limit(limit)
                .mapToDouble(t -> t.x())
                .forEach(qts::apply)

        var query = new DistanceProfileFunction.DistanceProfileQuery(ts, qts, windowSize)
        var dist = new MPX().apply(query)

        then:
        closeTo(dist.profile()[0], ERROR).matches(MPDist.load('mpdist_diff_size.csv', MPXSpec).x()[0])
    }

    @Override
    boolean equals(double[] a, double[] b) {
        return super.equals(a, b, ERROR)
    }

    static record MPQuery(double[] mp, int[] pi, double[] mpb, int[] pib) {
        static MPQuery load(String file) {
            return loadData(file, (rows) -> {
                double[] mp = new double[rows.length]
                int[] pi = new int[rows.length]
                double[] mpb = new double[rows.length]
                int[] mpib = new int[rows.length]
                int counter = 0
                Arrays.stream(rows).forEach(row -> {
                    mp[counter] = row[0] as double
                    pi[counter] = (row[1] as int) - 1 // in R indexing starts from 1
                    mpb[counter] = row[2] as double
                    mpib[counter++] = (row[3] as int) - 1 // in R indexing starts from 1
                })
                return new MPQuery(mp, pi, mpb, mpib)
            }, MPXSpec)
        }
    }

    static record Stats(double[] mean, double[] std) {}
}
