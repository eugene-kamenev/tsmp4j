package com.github.eugene.kamenev.tsmp4j.algo.mp.stomp

import com.github.eugene.kamenev.tsmp4j.BaseSpec
import com.github.eugene.kamenev.tsmp4j.algo.mp.BaseMatrixProfile
import com.github.eugene.kamenev.tsmp4j.stats.BaseRollingWindowStatistics

class STOMPSpec extends BaseSpec {

    def 'test stomp self join'() {
        given:
        var limit = 200
        var windowSize = 30
        var check = loadCheck('stomp_self_join.csv')

        when:
        var stomp = new STOMP(windowSize, limit)

        data.stream()
                .mapToDouble(t -> t.x())
                .limit(limit)
                .forEach(stomp::update)

        var mp = stomp.get()

        then:
        mp
    }

    def 'test stomp ab join'() {
        given:
        var limit = 200
        var windowSize = 30
        var check = loadCheck('stomp_ab_join.csv', true)

        when:
        var stomp = new STOMP(windowSize, limit)
        var query = new BaseRollingWindowStatistics(windowSize, 60)

        data.stream()
                .mapToDouble(t -> t.x())
                .limit(limit)
                .forEach(stomp::update)

        data.stream()
                .mapToDouble(t -> t.y())
                .limit(60)
                .forEach(query::apply)

        var mp = stomp.get(query)

        then:
        equals(mp.profile(), check.profile())
        equals(mp.indexes(), check.indexes())
    }

    def 'test stomp with NaN values in data'() {
        var limit = 200
        var windowSize = 30
        var check = loadCheck('stomp_self_join_nan.csv')
        var ts = data.stream()
                .mapToDouble(t -> t.x())
                .limit(limit)
                .toArray()

        ts[100] = Double.NaN

        when:
        var stomp = new STOMP(windowSize, limit)

        Arrays.stream(ts)
                .forEach(stomp::update)

        var mp = stomp.get()

        then:
        equals(mp.profile(), check.profile())
        equals(mp.indexes(), check.indexes())
        equals(mp.rightIndexes(), check.rightIndexes())
        equals(mp.rightProfile(), check.rightProfile())
        equals(mp.leftProfile(), check.leftProfile())
        equals(mp.leftIndexes(), check.leftIndexes())
    }

    static BaseMatrixProfile loadCheck(String fileName, boolean partial = false) {
        loadData(fileName, (rows) -> {
            var mp = new double[rows.length]
            var lmp = new double[rows.length]
            var rmp = new double[rows.length]
            var pi = new int[rows.length]
            var lpi = new int[rows.length]
            var rpi = new int[rows.length]
            var counter = 0
            Arrays.stream(rows).forEach(row -> {
                // mp,mpi,lmp,lpi,rmp,rpi
                mp[counter] = parseDouble(row[0])
                pi[counter] = parseInt(row[1])
                pi[counter] = pi[counter] >= 0 ? pi[counter] - 1 : -1
                if (!partial) {
                    lmp[counter] = parseDouble(row[2])
                    lpi[counter] = parseInt(row[3])
                    lpi[counter] = lpi[counter] >= 0 ? lpi[counter] - 1 : -1
                    rmp[counter] = parseDouble(row[4])
                    rpi[counter] = parseInt(row[5])
                    rpi[counter] = rpi[counter] >= 0 ? rpi[counter] - 1 : -1
                }
                counter++
            })
            return new BaseMatrixProfile(mp, pi, rmp, lmp, rpi, lpi)
        }, STOMPSpec)
    }
}
