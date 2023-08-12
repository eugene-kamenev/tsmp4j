/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package com.github.eugene.kamenev.tsmp4j.algo.mp;

import com.github.eugene.kamenev.tsmp4j.algo.mp.mass.MASS2;
import com.github.eugene.kamenev.tsmp4j.stats.BaseRollingWindowStatistics;
import com.github.eugene.kamenev.tsmp4j.stats.BaseWindowStatistic;
import com.github.eugene.kamenev.tsmp4j.stats.RollingWindowStatistics;
import com.github.eugene.kamenev.tsmp4j.stats.WindowStatistic;
import java.util.Arrays;

/**
 * STAMP: Streaming Time series Anytime Matrix Profile Reference: Yeh CCM, Zhu Y, Ulanova L, Begum
 * N, Ding Y, Dau HA, et al. Matrix profile I: All pairs similarity joins for time series: A
 * unifying view that includes motifs, discords and shapelets. Proc - IEEE Int Conf Data Mining,
 * ICDM. 2017;1317-22. Reference: Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time
 * Series Chains: A New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1-27.
 * Reference Website: <a
 * href="http://www.cs.ucr.edu/~eamonn/MatrixProfile.html">MatrixProfile.html</a>
 */
public class STAMP extends BaseMatrixProfileAlgorithm<BaseWindowStatistic> {

    private final DistanceProfileFunction<BaseWindowStatistic> distanceProfileFunction;

    public STAMP(int minInstances, int windowSize) {
        super(new BaseRollingWindowStatistics<>(windowSize, minInstances));
        this.distanceProfileFunction = new MASS2<>();
    }

    @Override
    public MatrixProfile get(double[] query) {
        if (this.isReady()) {
            var qs = new BaseRollingWindowStatistics<>(query.length,
                new BaseWindowStatistic[query.length]);
            for (double v : query) {
                qs.apply(v);
            }
            return stamp(qs, this.rollingStatistics, this.rollingStatistics.windowSize(),
                this.distanceProfileFunction);
        }
        return null;
    }

    @Override
    public MatrixProfile get() {
        if (this.isReady()) {
            return stamp(this.rollingStatistics, this.rollingStatistics.windowSize(),
                this.distanceProfileFunction);
        }
        return null;
    }

    public static <S extends WindowStatistic> MatrixProfile stamp(RollingWindowStatistics<S> tsA,
        int window, int w,
        RollingWindowStatistics<S> tsB, boolean trivialMatch, DistanceProfileFunction<S> ds) {
        var n = tsA.getStatsBuffer().getLength();
        var matrixProfile = new double[n - window + 1];
        var matrixProfileIndex = new int[n - window + 1];
        Arrays.fill(matrixProfile, Double.POSITIVE_INFINITY);
        Arrays.fill(matrixProfileIndex, -1);

        double[] distanceProfile;
        int[] distanceProfileIndex;

        var index = 0;
        while (index < w) {
            var dsq = new DistanceProfileFunction.DistanceProfileQuery<>(tsA, tsB, index, window);
            distanceProfile = ds.apply(dsq);
            distanceProfileIndex = getDistanceProfileIndex(n, index, window);

            if (trivialMatch) {
                int startIndex = Math.max(0, index - window / 2);
                int endIndex = Math.min(index + window / 2 + 1, distanceProfile.length);
                for (int i = startIndex; i < endIndex; i++) {
                    distanceProfile[i] = Double.POSITIVE_INFINITY;
                }
            }

            for (int i = 0; i < distanceProfile.length; i++) {
                if (distanceProfile[i] < matrixProfile[i]) {
                    matrixProfile[i] = distanceProfile[i];
                    matrixProfileIndex[i] = distanceProfileIndex[i];
                }
            }
            index++;
        }

        return new MatrixProfile(matrixProfile, matrixProfileIndex);
    }

    public static MatrixProfile of(double[] ts, double[] query) {
        var tsStats = new BaseRollingWindowStatistics<>(ts.length, ts.length);
        for (double v : ts) {
            tsStats.apply(v);
        }
        var queryStats = new BaseRollingWindowStatistics<>(query.length, query.length);
        for (double v : query) {
            queryStats.apply(v);
        }
        return stamp(tsStats, queryStats, query.length, new MASS2<>());
    }

    public static MatrixProfile of(double[] ts, int window) {
        var tsStats = new BaseRollingWindowStatistics<>(window, ts.length);
        for (double v : ts) {
            tsStats.apply(v);
        }
        return stamp(tsStats, window, new MASS2<>());
    }

    public static <S extends WindowStatistic> MatrixProfile stamp(RollingWindowStatistics<S> ts,
        int window, DistanceProfileFunction<S> ds) {
        return stamp(ts, window, ts.getStatsBuffer().getLength() - window + 1, ts, true, ds);
    }

    public static <S extends WindowStatistic> MatrixProfile stamp(RollingWindowStatistics<S> ts,
        RollingWindowStatistics<S> query, int window, DistanceProfileFunction<S> ds) {
        return stamp(query, window, ts.getStatsBuffer().getLength() - window + 1, ts, false, ds);
    }

    private static int[] getDistanceProfileIndex(int n, int index, int w) {
        var distanceProfileIndex = new int[n - w + 1];
        Arrays.fill(distanceProfileIndex, index);
        return distanceProfileIndex;
    }
}
