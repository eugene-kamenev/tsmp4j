package com.github.eugene.kamenev.tsmp4j.algo.mp.stamp;

import com.github.eugene.kamenev.tsmp4j.BaseTest;
import com.github.eugene.kamenev.tsmp4j.algo.mp.BaseMatrixProfile;
import com.github.eugene.kamenev.tsmp4j.stats.BaseRollingWindowStatistics;
import com.github.eugene.kamenev.tsmp4j.stats.BaseWindowStatistic;
import org.junit.jupiter.api.Test;

import java.util.Arrays;

class STAMPTest extends BaseTest {

    @Test
    void testStampSelfJoin() {
        var limit = 200;
        var windowSize = 30;
        var check = loadCheck("stamp_self_join.csv");

        var stamp = new STAMP(windowSize, limit);

        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(stamp::update);

        var mp = stamp.get();

        equals(mp.profile(), check.profile());
        equals(mp.indexes(), check.indexes());
        equals(mp.rightIndexes(), check.rightIndexes());
        equals(mp.rightProfile(), check.rightProfile());
        equals(mp.leftProfile(), check.leftProfile());
        equals(mp.leftIndexes(), check.leftIndexes());
    }

    @Test
    void testStampAbJoin() {
        var limit = 200;
        var windowSize = 30;
        var check = loadCheck("stamp_ab_join.csv", true);

        var stamp = new STAMP(windowSize, limit);
        BaseRollingWindowStatistics<BaseWindowStatistic> query = new BaseRollingWindowStatistics<>(windowSize, 60);

        data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .forEach(stamp::update);

        data.stream()
            .mapToDouble(t -> t.y())
            .limit(60)
            .forEach(query::apply);

        var mp = stamp.get(query);

        equals(mp.profile(), check.profile());
        equals(mp.indexes(), check.indexes());
    }

    @Test
    void testStampWithNaNValuesInData() {
        var limit = 200;
        var windowSize = 30;
        var check = loadCheck("stamp_self_join_nan.csv");
        var ts = data.stream()
            .mapToDouble(t -> t.x())
            .limit(limit)
            .toArray();

        ts[100] = Double.NaN;

        var stamp = new STAMP(windowSize, limit);

        Arrays.stream(ts)
            .forEach(stamp::update);

        var mp = stamp.get();

        equals(mp.profile(), check.profile());
        equals(mp.indexes(), check.indexes());
        equals(mp.rightIndexes(), check.rightIndexes());
        equals(mp.rightProfile(), check.rightProfile());
        equals(mp.leftProfile(), check.leftProfile());
        equals(mp.leftIndexes(), check.leftIndexes());
    }

    static BaseMatrixProfile loadCheck(String fileName) {
        return loadCheck(fileName, false);
    }

    static BaseMatrixProfile loadCheck(String fileName, boolean partial) {
        return loadData(fileName, rows -> {
            var mp = new double[rows.length];
            var lmp = new double[rows.length];
            var rmp = new double[rows.length];
            var pi = new int[rows.length];
            var lpi = new int[rows.length];
            var rpi = new int[rows.length];
            for (int counter = 0; counter < rows.length; counter++) {
                var row = rows[counter];
                // mp,mpi,lmp,lpi,rmp,rpi
                mp[counter] = parseDouble(row[0]);
                pi[counter] = parseInt(row[1]);
                pi[counter] = pi[counter] >= 0 ? pi[counter] - 1 : -1;
                if (!partial) {
                    lmp[counter] = parseDouble(row[2]);
                    lpi[counter] = parseInt(row[3]);
                    lpi[counter] = lpi[counter] >= 0 ? lpi[counter] - 1 : -1;
                    rmp[counter] = parseDouble(row[4]);
                    rpi[counter] = parseInt(row[5]);
                    rpi[counter] = rpi[counter] >= 0 ? rpi[counter] - 1 : -1;
                }
            }
            return new BaseMatrixProfile(0, 0d, mp, pi, rmp, lmp, rpi, lpi);
        }, STAMPTest.class);
    }
}
