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

import com.github.eugene.kamenev.tsmp4j.stats.RollingWindowStatistics;
import java.util.Arrays;

/**
 * Fast implementation of MatrixProfile and MatrixProfileIndex for internal purposes, without FFT
 * Reference: <a
 * href="https://github.com/matrix-profile-foundation/tsmp/blob/master/R/mpx.R">mpx.R</a>
 */
public class MPX extends BaseMatrixProfileAlgorithm {

    private final boolean crossCorrelation;
    private final int minlag;

    public MPX(int minInstances, int windowSize, boolean crossCorrelation) {
        this(minInstances, windowSize,
            new RollingWindowStatistics(windowSize, minInstances, true), crossCorrelation);
    }

    public MPX(int minInstances, int windowSize, RollingWindowStatistics rollingStatistics,
        boolean crossCorrelation) {
        super(minInstances, windowSize, rollingStatistics);
        this.crossCorrelation = crossCorrelation;
        this.minlag = (int) Math.ceil(windowSize / 4.0);
    }

    public static MatrixProfile of(double[] ts, int windowSize, boolean crossCorrelation) {
        var mpx = new MPX(ts.length, windowSize, crossCorrelation);
        Arrays.stream(ts)
            .forEach(mpx::update);
        return mpx.get();
    }

    public static MatrixProfile of(double[] ts, int windowSize) {
        var mpx = new MPX(ts.length, windowSize, false);
        Arrays.stream(ts)
            .forEach(mpx::update);
        return mpx.get();
    }

    private static void computeJoin(RollingWindowStatistics ts, RollingWindowStatistics query,
        int w,
        double[] diff_fa, double[] diff_ga, double[] diff_fb, double[] diff_gb, double[] mp,
        int[] mpi, double[] mpb, int[] mpib) {
        int profile_len = ts.getStatsBuffer().getLength() - w + 1;
        int profile_lenb = query.getStatsBuffer().getLength() - w + 1;
        double cov_, corr_;

        for (int i = 0; i < profile_len; i++) {
            int mx = Math.min((profile_len - i), profile_lenb);

            cov_ = 0;
            for (int j = i; j < i + w; j++) {
                cov_ += ((ts.getX(j) - ts.getMean(i)) * (query.getX(j - i) - query.getMean(0)));
            }

            for (int j = 0; j < mx; j++) {
                int k = j + i;
                cov_ += diff_fa[k] * diff_gb[j] + diff_ga[k] * diff_fb[j];
                corr_ = cov_ * ts.getStdDev(k) * query.getStdDev(j);

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
                mp[i] =
                    (mp[i] == -1.0) ? Double.POSITIVE_INFINITY : Math.sqrt(2.0 * w * (1.0 - mp[i]));
            } else {
                if (mp[i] > 1.0) {
                    mp[i] = 1.0;
                }
            }
        }
    }

    @Override
    public MatrixProfile get(double[] query) {
        if (this.isReady()) {
            int n = this.minInstances;
            int qn = query.length;
            int w = windowSize;

            int profile_len = n - w + 1;
            int profile_lenb = qn - w + 1;

            var ds = this.rollingStatistics;

            var qs = new RollingWindowStatistics(query.length, query.length, true);
            for (double v : query) {
                qs.apply(v);
            }

            double[] diff_fa = new double[profile_len];
            double[] diff_ga = new double[profile_len];
            double[] diff_fb = new double[profile_lenb];
            double[] diff_gb = new double[profile_lenb];

            double[] mp = new double[profile_len];
            int[] mpi = new int[profile_len];
            double[] mpb = new double[profile_lenb];
            int[] mpib = new int[profile_lenb];

            Arrays.fill(mp, -1.0);
            Arrays.fill(mpb, -1.0);

            for (int i = w; i < n; i++) {
                diff_fa[i - w + 1] = 0.5 * (ds.getX(i) - ds.getX(i - w));
                diff_ga[i - w + 1] =
                    (ds.getX(i) - ds.getMean(i - w + 1)) + (ds.getX(i - w) - ds.getMean(i - w));
                if (i < qn) {
                    diff_fb[i - w + 1] = 0.5 * (qs.getX(i) - qs.getX(i - w));
                    diff_gb[i - w + 1] = (query[i] - qs.getMean(i - w + 1) + (qs.getX(i - w)
                        - qs.getMean(
                        i - w)));
                }
            }
            // AB Join
            computeJoin(ds, qs, w, diff_fa, diff_ga, diff_fb, diff_gb, mp, mpi, mpb, mpib);
            // BA Join
            computeJoin(qs, ds, w, diff_fb, diff_gb, diff_fa, diff_ga, mpb, mpib, mp, mpi);

            postProcess(mp, w, crossCorrelation);
            postProcess(mpb, w, crossCorrelation);

            return new MatrixProfile(mp, mpi);
        }
        return null;
    }

    @Override
    public MatrixProfile get() {
        if (this.isReady()) {
            int n = this.minInstances;
            int w = windowSize;
            int profile_len = n - w + 1;
            int diag, offset, col;
            double c, c_cmp;

            double[] df = new double[profile_len];
            double[] dg = new double[profile_len];
            double[] mp = new double[profile_len];
            int[] mpi = new int[profile_len];

            var sb = this.rollingStatistics;

            for (int i = w; i < n; i++) {
                df[i - w + 1] = 0.5 * (sb.getX(i) - sb.getX(i - w));
                dg[i - w + 1] =
                    (sb.getX(i) - sb.getMean(i - w + 1)) + (sb.getX(i - w) - sb.getMean(i - w));
            }

            for (diag = this.minlag + 1; diag < profile_len; diag++) {
                c = 0;
                for (int i = diag, j = 0; j < w; i++, j++) {
                    c += (sb.getX(i) - sb.getMean(diag)) * (sb.getX(i - diag) - sb.getMean(0));
                }

                for (offset = 0; offset < n - w - diag + 1; offset++) {
                    col = offset + diag;
                    c += df[offset] * dg[col] + df[col] * dg[offset];
                    c_cmp = c * sb.getStdDev(offset) * sb.getStdDev(col);

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
}
