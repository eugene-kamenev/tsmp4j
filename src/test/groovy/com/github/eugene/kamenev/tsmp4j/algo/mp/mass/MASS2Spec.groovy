package com.github.eugene.kamenev.tsmp4j.algo.mp.mass

import com.github.eugene.kamenev.tsmp4j.BaseSpec
import com.github.eugene.kamenev.tsmp4j.algo.mp.DistanceProfileFunction
import com.github.eugene.kamenev.tsmp4j.stats.BaseRollingWindowStatistics

class MASS2Spec extends BaseSpec {

    def 'test mass2 distance profile'() {
        given:
        var windowSize = 51
        var skip = 149
        var limit = windowSize
        var ts = new BaseRollingWindowStatistics(windowSize, data.size())
        var qts = new BaseRollingWindowStatistics(windowSize, limit)
        var check = MPDist.load('mass_2_profile.csv', MASS2Spec)

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

        var dist = new MASS2().apply(query)

        then:
        equals(dist, check.x())

    }


}
