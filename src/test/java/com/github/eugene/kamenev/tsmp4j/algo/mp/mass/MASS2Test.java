package com.github.eugene.kamenev.tsmp4j.algo.mp.mass;

import com.github.eugene.kamenev.tsmp4j.BaseTest;
import com.github.eugene.kamenev.tsmp4j.algo.mp.DistanceProfileFunction;
import com.github.eugene.kamenev.tsmp4j.stats.BaseRollingWindowStatistics;
import com.github.eugene.kamenev.tsmp4j.stats.BaseWindowStatistic;
import org.junit.jupiter.api.Test;

class MASS2Test extends BaseTest {

    @Test
    void testMass2DistanceProfile() {
        var windowSize = 51;
        var skip = 149;
        var limit = windowSize;
        BaseRollingWindowStatistics<BaseWindowStatistic> ts =
            new BaseRollingWindowStatistics<>(windowSize, data.size());
        BaseRollingWindowStatistics<BaseWindowStatistic> qts =
            new BaseRollingWindowStatistics<>(windowSize, limit);
        var check = MPDist.load("mass_2_profile.csv", MASS2Test.class);

        data.stream()
            .mapToDouble(t -> t.y())
            .forEach(ts::apply);
        data.stream()
            .skip(skip)
            .limit(limit)
            .mapToDouble(t -> t.x())
            .forEach(qts::apply);

        var query = new DistanceProfileFunction.DistanceProfileQuery<>(ts, qts, windowSize);

        var dist = new MASS2<BaseWindowStatistic>().apply(query);

        equals(dist.profile(), check.x());
    }
}
