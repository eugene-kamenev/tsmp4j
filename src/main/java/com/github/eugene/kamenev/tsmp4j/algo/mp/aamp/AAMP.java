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

package com.github.eugene.kamenev.tsmp4j.algo.mp.aamp;

import com.github.eugene.kamenev.tsmp4j.algo.mp.BaseMatrixProfile;
import com.github.eugene.kamenev.tsmp4j.algo.mp.BaseMatrixProfileAlgorithm;
import com.github.eugene.kamenev.tsmp4j.algo.mp.MatrixProfile;
import com.github.eugene.kamenev.tsmp4j.stats.NoStatistic;
import com.github.eugene.kamenev.tsmp4j.stats.RollingWindowStatistics;
import com.github.eugene.kamenev.tsmp4j.stats.WindowStatistic;
import java.util.Arrays;

/**
 * Efficient Algorithms for Knowledge Discovery from Time Series
 * Tanmoy Mondal, Reza Akbarinia and Florent Masseglia
 * https://mondal-tanmoy.github.io/files/pdf/journal/AAMP_Journal.pdf
 * This algorithm computes matrix profile with the "pure" (non-normalized) Euclidean distance.
 * https://sites.google.com/view/aamp-and-acamp/home
 * https://github.com/anoynymREVIEW/ICDM_AAMP_ACAMP
 */
public class AAMP extends BaseMatrixProfileAlgorithm<NoStatistic, MatrixProfile> {

    private final double p;

    public AAMP(RollingWindowStatistics<NoStatistic> rollingWindowStatistics,
        double exclusionZone, double p) {
        super(rollingWindowStatistics, exclusionZone);
        this.p = p;
    }

    public AAMP(RollingWindowStatistics<NoStatistic> rollingWindowStatistics, double p) {
        this(rollingWindowStatistics, 0.5d, p);
    }

    @Override
    public MatrixProfile get(RollingWindowStatistics<NoStatistic> query) {
        if (this.isReady()) {
            return aamp(this.rollingStatistics(), query, this.p);
        }
        return null;
    }

    /**
     * Computes matrix profile join between two time series with the "pure" (non-normalized)
     * distance. As in {@code STAMP} and {@code STOMP}, the exclusion zone is not applicable to a
     * join between two different series, hence it is not used here.
     *
     * @param ts reference time series
     * @param query query time series
     * @param p the power of the distance
     * @return matrix profile of {@code ts} against {@code query} in
     *     {@code profile}/{@code indexes} and matrix profile of {@code query} against {@code ts}
     *     in {@code leftProfile}/{@code leftIndexes}
     */
    public static <S extends WindowStatistic> MatrixProfile aamp(RollingWindowStatistics<S> ts,
        RollingWindowStatistics<S> query, double p) {
        int m = ts.windowSize();
        if (query.dataSize() < m) {
            throw new IllegalArgumentException(
                "Query must be at least as long as the window size.");
        }

        double[] mp = new double[ts.dataSize() - m + 1];
        int[] mpi = new int[mp.length];
        double[] mpb = new double[query.dataSize() - m + 1];
        int[] mpib = new int[mpb.length];
        Arrays.fill(mp, Double.POSITIVE_INFINITY);
        Arrays.fill(mpb, Double.POSITIVE_INFINITY);

        // AB Join
        join(ts, query, mp, mpi, mpb, mpib, m, p);
        // BA Join
        join(query, ts, mpb, mpib, mp, mpi, m, p);

        return new BaseMatrixProfile(m, 0d, normalize(mp, p), mpi, null, normalize(mpb, p), null, mpib);
    }

    private static <S extends WindowStatistic> void join(RollingWindowStatistics<S> a,
        RollingWindowStatistics<S> b, double[] mp, int[] mpi, double[] mpb, int[] mpib, int w, double p) {
        int amx = a.dataSize() - w + 1;
        int bmx = b.dataSize() - w + 1;
        for (int ia = 0; ia < amx; ia++) {
            int mx = Math.min(amx - ia, bmx);
            double d = 0;
            for (int k = 0; k < w; k++) {
                d += Math.pow(Math.abs(a.x(ia + k) - b.x(k)), p);
            }
            for (int ib = 0; ib < mx; ib++) {
                if (ib > 0) {
                    d = d - Math.pow(Math.abs(a.x(ia + ib - 1) - b.x(ib - 1)), p) +
                        Math.pow(Math.abs(a.x(ia + ib - 1 + w) - b.x(ib - 1 + w)), p);
                }
                if (d < mp[ia + ib]) {
                    mp[ia + ib] = d;
                    mpi[ia + ib] = ib;
                }
                if (d < mpb[ib]) {
                    mpb[ib] = d;
                    mpib[ib] = ia + ib;
                }
            }
        }
    }

    private static double[] normalize(double[] dmin, double p) {
        for (int i = 0; i < dmin.length; i++) {
            dmin[i] = Math.max(dmin[i], 0);
            dmin[i] = Math.pow(dmin[i], 1.0 / p);
        }
        return dmin;
    }

    @Override
    public MatrixProfile get() {
        var X = this.rollingStatistics();
        int m = X.windowSize();
        int Nb = this.rollingStatistics().dataSize();
        int s = Nb - m;
        int excZone = this.exclusionZoneSize;
        double[] Dmin = new double[s + 1];
        int[] minind = new int[s + 1];
        Arrays.fill(Dmin, Double.POSITIVE_INFINITY);
        boolean matchFlag = false;
        for (int k = 0; k < s; k++) {
            double D = 0;
            for (int j = 0; j < m; j++) {
                D += Math.pow(Math.abs(X.x(j) - X.x(k + j + 1)), p);
            }

            if (k >= excZone) {
                matchFlag = true;
            }

            if (D < Dmin[0] && matchFlag) {
                Dmin[0] = D;
                minind[0] = k + 1;
            }

            if (D < Dmin[k + 1] && matchFlag) {
                Dmin[k + 1] = D;
                minind[k + 1] = 0;
            }

            for (int i = 1; i < s - k; i++) {
                int kplusi = k + i + 1;
                D = D - Math.pow(Math.abs(X.x(i - 1) - X.x(kplusi - 1)), p) +
                    Math.pow(Math.abs(X.x(m + i - 1) - X.x(m + kplusi - 1)), p);

                if (Dmin[i] > D && matchFlag) {
                    minind[i] = kplusi;
                    Dmin[i] = D;
                }

                if (Dmin[kplusi] > D && matchFlag) {
                    minind[kplusi] = i;
                    Dmin[kplusi] = D;
                }
            }
        }

        return new BaseMatrixProfile(m, exclusionZone, normalize(Dmin, p), minind);
    }
}
