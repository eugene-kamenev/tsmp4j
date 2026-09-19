package com.github.eugene.kamenev.tsmp4j.algo.mp.aamp;

import com.github.eugene.kamenev.tsmp4j.BaseTest;
import com.github.eugene.kamenev.tsmp4j.stats.RollingWindowWithoutStatistics;
import org.junit.jupiter.api.Test;

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
}
