package com.github.eugene.kamenev.tsmp4j.algo.cp

import com.github.eugene.kamenev.tsmp4j.BaseSpec
import com.github.eugene.kamenev.tsmp4j.algo.mp.mpx.MPXRollingWindowStatistics
import com.github.eugene.kamenev.tsmp4j.stats.BaseRollingWindowStatistics

class ContrastProfileSpec extends BaseSpec {

    static double ERROR = Math.pow(10, -11);

    def 'test contrast profile'() {
        given:
        var error = Math.pow(10, -14);
        var cmpCheck = MPDist.load("cmp.csv", ContrastProfileSpec)
        var platoCheck = MPDist.load("cmp_plato.csv", ContrastProfileSpec)
        var platoTwinCheck = MPDist.load("cmp_plato_twin.csv", ContrastProfileSpec)
        var windowSize = 30
        var positiveTs = new MPXRollingWindowStatistics(windowSize, data.size())
        var negativeTs = new MPXRollingWindowStatistics(windowSize, data.size())

        when:
        data.stream()
            .mapToDouble(t -> t.y())
            .forEach(positiveTs::apply)
        data.stream()
            .mapToDouble(t -> t.x())
            .forEach(negativeTs::apply)

        var profile = new ContrastProfileAlgorithm().apply(positiveTs, negativeTs)

        then:
        equals(profile.profile(), cmpCheck.x(), error)
        equals(profile.plato(), platoCheck.x(), error)
        equals(profile.platoTwin(), platoTwinCheck.x(), error)

    }

    def 'test pan contrast profile'() {
        given:
        var positiveTs = data.stream()
            .mapToDouble(t -> t.y())
            .toArray()
        var negativeTs = data.stream()
            .mapToDouble(t -> t.x())
            .toArray()

        when:
        var profile = new PanContrastProfileAlgorithm(10, 30, 10).apply(positiveTs, negativeTs)

        then:
        profile.profile().length == 11

    }

    def 'test relative frequency matrix profile'() {
        given:
        var rfmpCheck = loadData("rfmp_profile.csv", (rows) -> {
            var profile = new double[rows[0].length][rows.length];
            for (int i = 0; i < rows.length; i++) {
                for (int j = 0; j < rows[0].length; j++) {
                    profile[j][i] = parseDouble(rows[i][j]);
                }
            }
            return profile;
        }, ContrastProfileSpec)
        var rfmpIndexesCheck = loadData("rfmp_indexes.csv", (rows) -> {
            var indexes = new int[rows[0].length][rows.length];
            for (int i = 0; i < rows.length; i++) {
                for (int j = 0; j < rows[0].length; j++) {
                    indexes[j][i] = parseInt(rows[i][j]) - 1; // in R and Matlab indexing starts from 1
                }
            }
            return indexes;
        }, ContrastProfileSpec)
        var positiveTs = data.stream()
            .mapToDouble(t -> t.y())
            .toArray()
        var negativeTs = data.stream()
            .mapToDouble(t -> t.x())
            .toArray()

        when:
        var profile = new RelativeFrequencyMatrixProfileAlgorithm(3, 30, false, false).apply(positiveTs, negativeTs)

        then:
        for (int i = 0; i < profile.profile().length; i++) {
            equals(profile.indexes()[i], rfmpIndexesCheck[i])
            equals(profile.profile()[i], rfmpCheck[i])
        }
    }

    def 'test relative frequency contrast profile'() {
        given:
        var rfcpPlatoCheck = loadData("rfcp_plato.csv", (rows) -> {
            var plato = new double[rows.length];
            for (int i = 0; i < rows.length; i++) {
                plato[i] = parseDouble(rows[i][0]);
            }
            return plato;
        }, ContrastProfileSpec)
        var rfcpCpCheck = loadData("rfcp_cp.csv", (rows) -> {
            var cp = new double[rows[0].length][rows.length];
            for (int i = 0; i < rows.length; i++) {
                for (int j = 0; j < rows[0].length; j++) {
                    cp[j][i] = parseDouble(rows[i][j]);
                }
            }
            return cp;
        }, ContrastProfileSpec)
        var rfcpIndexCheck = loadData("rfcp_indexes.csv", (rows) -> {
            var index = new int[rows[0].length][rows.length];
            for (int i = 0; i < rows.length; i++) {
                for (int j = 0; j < rows[0].length; j++) {
                    index[j][i] = parseInt(rows[i][j]) - 1; // in R and Matlab indexing starts from 1
                }
            }
            return index;
        }, ContrastProfileSpec)
        var positiveTs = data.stream()
                .mapToDouble(t -> t.y())
                .toArray()
        var negativeTs = data.stream()
                .mapToDouble(t -> t.x())
                .toArray()
        when:
        var profile = new RelativeFrequencyContrastProfileAlgorithm(30, 3, false).apply(positiveTs, negativeTs)

        then:
        for (int i = 0; i < profile.indexes().length; i++) {
            equals(profile.indexes()[i], rfcpIndexCheck[i])
            equals(profile.profile()[i], rfcpCpCheck[i])
        }
        equals(profile.plato(), rfcpPlatoCheck)
    }

    @Override
    boolean equals(double[] a, double[] b) {
        return super.equals(a, b, ERROR)
    }
}
