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

package com.github.eugene.kamenev.tsmp4j.algo.mp.mpx;

import com.github.eugene.kamenev.tsmp4j.algo.mp.BaseMatrixProfileAlgorithm;
import com.github.eugene.kamenev.tsmp4j.algo.mp.DistanceProfileFunction;
import com.github.eugene.kamenev.tsmp4j.algo.mp.MatrixProfile;
import com.github.eugene.kamenev.tsmp4j.stats.RollingWindowStatistics;
import com.github.eugene.kamenev.tsmp4j.stats.WindowStatistic;
import java.util.Arrays;

/**
 * Fast implementation of MatrixProfile and MatrixProfileIndex for internal purposes, without FFT
 * Reference: <a
 * href="https://github.com/matrix-profile-foundation/tsmp/blob/master/R/mpx.R">mpx.R</a>
 */
public class MPX extends BaseMatrixProfileAlgorithm<MPXStatistics> implements
    DistanceProfileFunction<MPXStatistics> {

    private final boolean crossCorrelation;

    private final int minlag;

    public MPX(int windowSize, int bufferSize, boolean crossCorrelation) {
        this(new MPXRollingWindowStatistics(windowSize, bufferSize),
            (int) Math.ceil(windowSize / 4.0), crossCorrelation);
    }

    public MPX(RollingWindowStatistics<MPXStatistics> rollingWindowStatistics,
        int minlag, boolean crossCorrelation) {
        super(rollingWindowStatistics);
        this.minlag = minlag;
        this.crossCorrelation = crossCorrelation;
    }

    @Override
    public MatrixProfile get(double[] query) {
        if (this.isReady()) {
            var qs = new MPXRollingWindowStatistics(query.length, query.length);
            for (double v : query) {
                qs.apply(v);
            }
            return compute(this.rollingStatistics, qs, crossCorrelation);
        }
        return null;
    }

    @Override
    public MatrixProfile get() {
        if (this.isReady()) {
            var sb = ((MPXRollingWindowStatistics) this.rollingStatistics);

            int n = sb.dataSize();
            int w = sb.windowSize();
            int profile_len = n - w + 1;
            int diag, offset, col;
            double c, c_cmp;

            double[] mp = new double[profile_len];
            int[] mpi = new int[profile_len];

            for (diag = this.minlag + 1; diag < profile_len; diag++) {
                c = 0;
                for (int i = diag, j = 0; j < w; i++, j++) {
                    c += (sb.x(i) - sb.mean(diag)) * (sb.x(i - diag) - sb.mean(0));
                }
                for (offset = 0; offset < n - w - diag + 1; offset++) {
                    col = offset + diag;
                    c += sb.df(offset) * sb.dg(col) + sb.df(col) * sb.dg(offset);
                    c_cmp = c * sb.stdDev(offset) * sb.stdDev(col);

                    if (c_cmp > mp[offset]) {
                        mp[offset] = c_cmp;
                        mpi[offset] = col;
                    }

                    if (c_cmp > mp[col]) {
                        if (c_cmp > 1.0) {
                            c_cmp = 1.0;
                        }
                        mp[col] = c_cmp;
                        mpi[col] = offset;
                    }
                }
            }

            if (!this.crossCorrelation) {
                for (int i = 0; i < profile_len; i++) {
                    mp[i] = Math.sqrt(2.0 * w * (1.0 - mp[i]));
                }
            }
            return new MatrixProfile(mp, mpi);
        }
        return null;
    }

    @Override
    public double[] apply(DistanceProfileQuery<MPXStatistics> mpxStatisticsDistanceProfileQuery) {
        return compute(mpxStatisticsDistanceProfileQuery.ts(),
            mpxStatisticsDistanceProfileQuery.query(), this.crossCorrelation).profile();
    }

    public static MatrixProfile of(double[] ts, int windowSize, boolean crossCorrelation) {
        var mpx = new MPX(windowSize, ts.length, crossCorrelation);
        Arrays.stream(ts)
            .forEach(mpx::update);
        return mpx.get();
    }

    public static MatrixProfile of(double[] ts, int windowSize) {
        var mpx = new MPX(windowSize, ts.length, false);
        Arrays.stream(ts)
            .forEach(mpx::update);
        return mpx.get();
    }

    private static MatrixProfile compute(RollingWindowStatistics<MPXStatistics> ts,
        RollingWindowStatistics<MPXStatistics> qs, boolean crossCorrelation) {

        int n = ts.dataSize();
        int qn = qs.dataSize();
        int w = ts.windowSize();

        int profile_len = n - w + 1;
        int profile_lenb = qn - w + 1;

        double[] mp = new double[profile_len];
        int[] mpi = new int[profile_len];
        double[] mpb = new double[profile_lenb];
        int[] mpib = new int[profile_lenb];

        Arrays.fill(mp, -1.0);
        Arrays.fill(mpb, -1.0);

        // AB Join
        computeJoin(ts, qs, w, mp, mpi, mpb, mpib);
        // BA Join
        computeJoin(qs, ts, w, mpb, mpib, mp, mpi);

        postProcess(mp, w, crossCorrelation);
        postProcess(mpb, w, crossCorrelation);

        return new MatrixProfile(mp, mpi);
    }

    private static <S extends WindowStatistic> void computeJoin(RollingWindowStatistics<S> ts,
        RollingWindowStatistics<S> query,
        int w, double[] mp, int[] mpi, double[] mpb, int[] mpib) {
        int profile_len = ts.getStatsBuffer().getLength() - w + 1;
        int profile_lenb = query.getStatsBuffer().getLength() - w + 1;
        double cov_, corr_;

        var sa = (MPXRollingWindowStatistics) ts;
        var sb = (MPXRollingWindowStatistics) query;

        for (int i = 0; i < profile_len; i++) {
            int mx = Math.min((profile_len - i), profile_lenb);

            cov_ = 0;
            for (int j = i; j < i + w; j++) {
                cov_ += ((ts.x(j) - ts.mean(i)) * (query.x(j - i) - query.mean(0)));
            }

            for (int j = 0; j < mx; j++) {
                int k = j + i;
                cov_ += sa.df(k) * sb.dg(j) + sa.dg(k) * sb.df(j);
                corr_ = cov_ * ts.stdDev(k) * query.stdDev(j);

                if (corr_ > mp[k]) {
                    mp[k] = corr_;
                    mpi[k] = j;
                }

                if (corr_ > mpb[j]) {
                    mpb[j] = corr_;
                    mpib[j] = k;
                }
            }
        }
    }

    private static void postProcess(double[] mp, int w, boolean crossCorrelation) {
        int len = mp.length;
        for (int i = 0; i < len; i++) {
            if (!crossCorrelation) {
                if (mp[i] == -1.0) {
                    mp[i] = Double.POSITIVE_INFINITY;
                } else {
                    var value = Math.sqrt(2.0 * w * (1.0 - mp[i]));
                    if (Double.isNaN(value)) {
                        mp[i] = 0;
                    } else {
                        mp[i] = value;
                    }
                }
            } else {
                if (mp[i] > 1.0) {
                    mp[i] = 1.0;
                }
            }
        }
    }
}
