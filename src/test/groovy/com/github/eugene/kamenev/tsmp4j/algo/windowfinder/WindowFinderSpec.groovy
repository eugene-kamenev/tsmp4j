package com.github.eugene.kamenev.tsmp4j.algo.windowfinder

import com.github.eugene.kamenev.tsmp4j.BaseSpec

class WindowFinderSpec extends BaseSpec {

    def 'test windowfinder'() {
        when:
        var ts = data.stream()
                .mapToDouble(t -> t.x())
                .toArray()
        var result = MWF.mwf(ts, 10, 500)

        then:
        result == 45
    }
}
