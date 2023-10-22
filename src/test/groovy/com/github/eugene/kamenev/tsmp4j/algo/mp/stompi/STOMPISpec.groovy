package com.github.eugene.kamenev.tsmp4j.algo.mp.stompi

import com.github.eugene.kamenev.tsmp4j.BaseSpec
import com.github.eugene.kamenev.tsmp4j.algo.mp.stomp.STOMP
import com.github.eugene.kamenev.tsmp4j.stats.BaseRollingWindowStatistics

class STOMPISpec extends BaseSpec {

    public static double ERROR = Math.pow(10, -10)

    def 'test incremental stompi'() {
        when:
        var windowSize = 30
        var initialStats = new BaseRollingWindowStatistics(windowSize, 200)

        var initialData = data.stream()
                .mapToDouble(t -> t.x())
                .limit(200)
                .toArray()

        var stomp = new STOMP(windowSize, 300)

        data.stream()
                .mapToDouble(t -> t.x())
                .limit(300)
                .forEach(stomp::update)

        var mp = stomp.get()

        for (int i = 0; i < initialData.length; i++) {
            initialStats.apply(initialData[i])
        }

        var testData = data.stream()
                .mapToDouble(t -> t.x())
                .skip(200)
                .limit(100)
                .toArray()

        var stompi1 = new STOMPI(initialStats, 0)
        var stompi2 = new STOMPI(initialStats, 0)
        var stompi3 = new STOMPI(initialStats, 100)

        var incremental = stompi1.get()
        var withHistory = stompi3.get()
        for (int i = 0; i < testData.length; i++) {
            stompi1.update(testData[i])
            stompi2.update(testData[i])
            stompi3.update(testData[i])
            withHistory = stompi3.get()
            incremental = stompi1.get()
        }

        var incrementalBatch = stompi2.get()

        then:
        equals(incremental.indexes(), mp.indexes())
        equals(incremental.leftIndexes(), mp.leftIndexes())
        equals(incremental.rightIndexes(), mp.rightIndexes())

        equals(incremental.profile(), mp.profile(), ERROR)
        equals(incremental.leftProfile(), mp.leftProfile(), ERROR)
        equals(incremental.rightProfile(), mp.rightProfile(), ERROR)

        equals(incrementalBatch.indexes(), mp.indexes())
        equals(incrementalBatch.leftIndexes(), mp.leftIndexes())
        equals(incrementalBatch.rightIndexes(), mp.rightIndexes())

        equals(incrementalBatch.profile(), mp.profile())
        equals(incrementalBatch.leftProfile(), mp.leftProfile())
        equals(incrementalBatch.rightProfile(), mp.rightProfile())
    }
}
