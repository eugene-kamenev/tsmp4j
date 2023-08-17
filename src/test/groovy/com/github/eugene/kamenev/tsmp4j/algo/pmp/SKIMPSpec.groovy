package com.github.eugene.kamenev.tsmp4j.algo.pmp

import com.github.eugene.kamenev.tsmp4j.BaseSpec
import com.github.eugene.kamenev.tsmp4j.utils.Util

class SKIMPSpec extends BaseSpec {

    def 'test skimp no cross correlation'() {
        when:
        var limit = 200
        var windows = Util.createRange(4, 6, 1)
        var skimp = new SKIMP(limit, false, windows)
        data.stream()
                .mapToDouble(t -> t.x())
                .limit(limit)
                .forEach(skimp::update)
        var result = skimp.get()
        var dst = result.profile()
        var ids = result.indexes()

        then:
        dst.length == 3
        ids.length == 3
        dst[0].length == limit - 3
        dst[1].length == limit - 4
        dst[2].length == limit - 5
        ids[0].length == limit - 3
        ids[1].length == limit - 4
        ids[2].length == limit - 5
    }

    def 'test skimp cross correlation'() {
        when:
        var limit = 200
        var windows = Util.createRange(4, 6, 1)
        var skimp = new SKIMP(limit, true, windows)
        data.stream()
                .mapToDouble(t -> t.x())
                .limit(limit)
                .forEach(skimp::update)
        var result = skimp.get()
        var dst = result.profile()
        var ids = result.indexes()

        then:
        dst.length == 3
        ids.length == 3
        dst[0].length == limit - 3
        dst[1].length == limit - 4
        dst[2].length == limit - 5
        ids[0].length == limit - 3
        ids[1].length == limit - 4
        ids[2].length == limit - 5
    }
}
