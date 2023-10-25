package com.github.eugene.kamenev.tsmp4j.algo.fluss

import com.github.eugene.kamenev.tsmp4j.BaseSpec
import com.github.eugene.kamenev.tsmp4j.algo.mp.stomp.STOMP

class FLUSSSpec extends BaseSpec {

    private static final double[] flussData = loadData('fluss_data_tilt_abp.csv', (rows) -> {
        return Arrays.stream(rows)
                .mapToDouble(s -> {
                    return (s[0] as double)
                }).toArray()
    }, FLUSSSpec)

    private static final double[] flussCacCheck = loadData('fluss_cac.csv', (rows) -> {
        return Arrays.stream(rows)
                .mapToDouble(s -> {
                    return (s[0] as double)
                }).toArray()
    }, FLUSSSpec);

    def 'test fluss'() {
        given:
        int windowSize = 10
        int buffSize = 1000
        int numSegments = 2

        and:
        var stomp = new STOMP(windowSize, buffSize)

        Arrays.stream(flussData)
                .limit(buffSize)
                .forEach(stomp::update)

        var mp = stomp.get()

        when:
        var fluss = new FLUSS(windowSize, numSegments).apply(mp)

        then:
        equals(fluss.cac(), flussCacCheck)
        equals(fluss.changePoints(), new int[]{940, 874})
    }
}
