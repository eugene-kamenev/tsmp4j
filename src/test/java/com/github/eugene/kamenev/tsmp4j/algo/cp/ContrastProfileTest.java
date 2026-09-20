package com.github.eugene.kamenev.tsmp4j.algo.cp;

import com.github.eugene.kamenev.tsmp4j.BaseTest;
import com.github.eugene.kamenev.tsmp4j.algo.mp.mpx.MPXRollingWindowStatistics;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class ContrastProfileTest extends BaseTest {

    private static final double ERROR = Math.pow(10, -11);

    @Override
    protected double error() {
        return ERROR;
    }

    @Test
    void testContrastProfile() {
        var cmpCheck = MPDist.load("cmp.csv", ContrastProfileTest.class);
        var platoCheck = MPDist.load("cmp_plato.csv", ContrastProfileTest.class);
        var platoTwinCheck = MPDist.load("cmp_plato_twin.csv", ContrastProfileTest.class);
        var windowSize = 30;
        var positiveTs = new MPXRollingWindowStatistics(windowSize, data.size());
        var negativeTs = new MPXRollingWindowStatistics(windowSize, data.size());

        data.stream()
            .mapToDouble(t -> t.y())
            .forEach(positiveTs::apply);
        data.stream()
            .mapToDouble(t -> t.x())
            .forEach(negativeTs::apply);

        var profile = new ContrastProfileAlgorithm().apply(positiveTs, negativeTs);

        equals(profile.profile(), cmpCheck.x());
        equals(profile.plato(), platoCheck.x());
        equals(profile.platoTwin(), platoTwinCheck.x());
    }

    @Test
    void testPanContrastProfile() {
        var positiveTs = data.stream()
            .mapToDouble(t -> t.y())
            .toArray();
        var negativeTs = data.stream()
            .mapToDouble(t -> t.x())
            .toArray();

        var profile = new PanContrastProfileAlgorithm(10, 30, 10).apply(positiveTs, negativeTs);

        assertEquals(11, profile.profile().length);
    }

    @Test
    void testRelativeFrequencyMatrixProfile() {
        var rfmpCheck = loadData("rfmp_profile.csv", rows -> {
            var profile = new double[rows[0].length][rows.length];
            for (int i = 0; i < rows.length; i++) {
                for (int j = 0; j < rows[0].length; j++) {
                    profile[j][i] = parseDouble(rows[i][j]);
                }
            }
            return profile;
        }, ContrastProfileTest.class);
        var rfmpIndexesCheck = loadData("rfmp_indexes.csv", rows -> {
            var indexes = new int[rows[0].length][rows.length];
            for (int i = 0; i < rows.length; i++) {
                for (int j = 0; j < rows[0].length; j++) {
                    indexes[j][i] = parseInt(rows[i][j]) - 1; // in R and Matlab indexing starts from 1
                }
            }
            return indexes;
        }, ContrastProfileTest.class);
        var positiveTs = data.stream()
            .mapToDouble(t -> t.y())
            .toArray();
        var negativeTs = data.stream()
            .mapToDouble(t -> t.x())
            .toArray();

        var profile = new RelativeFrequencyMatrixProfileAlgorithm(3, 30, false, false).apply(positiveTs, negativeTs);

        for (int i = 0; i < profile.profile().length; i++) {
            equals(profile.indexes()[i], rfmpIndexesCheck[i]);
            equals(profile.profile()[i], rfmpCheck[i]);
        }
    }

    @Test
    void testRelativeFrequencyContrastProfile() {
        var rfcpPlatoCheck = loadData("rfcp_plato.csv", rows -> {
            var plato = new double[rows.length];
            for (int i = 0; i < rows.length; i++) {
                plato[i] = parseDouble(rows[i][0]);
            }
            return plato;
        }, ContrastProfileTest.class);
        var rfcpCpCheck = loadData("rfcp_cp.csv", rows -> {
            var cp = new double[rows[0].length][rows.length];
            for (int i = 0; i < rows.length; i++) {
                for (int j = 0; j < rows[0].length; j++) {
                    cp[j][i] = parseDouble(rows[i][j]);
                }
            }
            return cp;
        }, ContrastProfileTest.class);
        var rfcpIndexCheck = loadData("rfcp_indexes.csv", rows -> {
            var index = new int[rows[0].length][rows.length];
            for (int i = 0; i < rows.length; i++) {
                for (int j = 0; j < rows[0].length; j++) {
                    index[j][i] = parseInt(rows[i][j]) - 1; // in R and Matlab indexing starts from 1
                }
            }
            return index;
        }, ContrastProfileTest.class);
        var positiveTs = data.stream()
            .mapToDouble(t -> t.y())
            .toArray();
        var negativeTs = data.stream()
            .mapToDouble(t -> t.x())
            .toArray();

        var profile = new RelativeFrequencyContrastProfileAlgorithm(30, 3, false).apply(positiveTs, negativeTs);

        for (int i = 0; i < profile.indexes().length; i++) {
            equals(profile.indexes()[i], rfcpIndexCheck[i]);
            equals(profile.profile()[i], rfcpCpCheck[i]);
        }
        equals(profile.plato(), rfcpPlatoCheck);
    }
}
