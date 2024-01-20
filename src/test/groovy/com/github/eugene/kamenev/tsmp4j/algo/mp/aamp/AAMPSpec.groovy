package com.github.eugene.kamenev.tsmp4j.algo.mp.aamp

import com.github.eugene.kamenev.tsmp4j.BaseSpec
import com.github.eugene.kamenev.tsmp4j.stats.NoStatistic
import com.github.eugene.kamenev.tsmp4j.stats.RollingWindowWithoutStatistics

class AAMPSpec extends BaseSpec {

        def 'test aamp distance profile'() {
            given:
            var windowSize = 30
            var limit = 200
            var checkMp = loadMP('aamp_self_join.csv', AAMPSpec)

            when:
            var aamp = new AAMP(new RollingWindowWithoutStatistics(windowSize, new NoStatistic[limit]), 2)

            data.stream()
                .limit(limit)
                .forEach(t -> {
                    aamp.update(t.x())
                })

            var mp = aamp.get()

            then:
            equals(checkMp.mp(), mp.profile())
            equals(checkMp.pi(), mp.indexes())
        }
}
