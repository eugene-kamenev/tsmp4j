package com.github.eugene.kamenev.tsmp4j.algo.tsc;

import com.github.eugene.kamenev.tsmp4j.BaseSpec

class TSCSpec extends BaseSpec {

    public static final double[] data = load1D("data.csv", TSCSpec)
    public static final int[] validChain = new int[] {5122,5318,5518,5917,6312,6507,6700,7457,8034,8613,9396}

    def 'test robust time series chains'() {
        when:
        def tsc = new TSC(180, data.length, 0.25d)
        Arrays.stream(data)
            .forEach(tsc::update)
        def chains = tsc.get()
        def score = TSC.bestScore(data, chains, 180);

        then:
        score.index() == 9396
        score.score() == 4.0d
        equals(score.predIdxOurs(), validChain)
    }
}