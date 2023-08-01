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

package com.github.eugene.kamenev.tsmp4j;

import java.util.Arrays;

/**
 * Fast implementation of MatrixProfile and MatrixProfileIndex for internal purposes, without FFT
 * Reference: <a href="https://github.com/matrix-profile-foundation/tsmp/blob/master/R/mpx.R">mpx.R</a>
 *
 * @param mp Matrix Profile
 * @param mpi Matrix Profile Index
 */
public record MPX(double[] mp, int[] mpi, double[] mpb, int[] mpib) {

    public static MPX mpx(double[] ts, int w, boolean crossCorrelation) {
        var minlag = (int) Math.ceil(w / 4.0);
        return mpx(ts, w, minlag, crossCorrelation);
    }

    public static MPX mpx(double[] ts, int w, int minlag, boolean crossCorrelation) {
        var n = ts.length;
        var profile_len = n - w + 1;
        var stats = muinvn(ts, w);
        var mu = stats[0];
        var sig = stats[1];
        var df = new double[profile_len];
        var dg = new double[profile_len];
        var mp = new double[profile_len];
        var mpi = new int[profile_len];

        for (var i = w; i < n; i++) {
            df[i - w + 1] = 0.5 * (ts[i] - ts[i - w]);
            dg[i - w + 1] = (ts[i] - mu[i - w + 1]) + (ts[i - w] - mu[i - w]);
        }

        double c, c_cmp;
        for (int diag = minlag + 1; diag < profile_len; diag++) {
            c = 0;
            for (int i = diag; i < diag + w; i++) {
                c = c + ((ts[i] - mu[diag]) * (ts[i - diag] - mu[0]));
            }

            for (int offset = 0; offset < n - w - diag + 1; offset++) {
                int col = offset + diag;
                c = c + df[offset] * dg[col] + df[col] * dg[offset];
                c_cmp = c * sig[offset] * sig[col];

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

        if (!crossCorrelation) {
            for (int i = 0; i < profile_len; i++) {
                mp[i] = Math.sqrt(2.0 * w * (1.0 - mp[i]));
            }
        }

        return new MPX(mp, mpi, null, null);
    }

    public static MPX mpx(double[] ts, double[] query, int w, boolean crossCorrelation) {
        int n = ts.length;
        int qn = query.length;
        double cov_, corr_;

        int profile_len = n - w + 1;
        int profile_lenb = qn - w + 1;

        var stats_a = muinvn(ts, w);
        double[] mua = stats_a[0];
        double[] siga = stats_a[1];

        var stats_b = muinvn(query, w);
        double[] mub = stats_b[0];
        double[] sigb = stats_b[1];

        double[] diff_fa = new double[profile_len];
        double[] diff_ga = new double[profile_len];
        double[] diff_fb = new double[profile_lenb];
        double[] diff_gb = new double[profile_lenb];

        double[] mp = new double[profile_len];
        int[] mpi = new int[profile_len];
        double[] mpb = new double[profile_lenb];
        int[] mpib = new int[profile_lenb];

        Arrays.fill(mp, -1.0);
        Arrays.fill(mpi, -1);
        Arrays.fill(mpb, -1.0);
        Arrays.fill(mpib, -1);

        diff_fa[0] = 0;
        diff_ga[0] = 0;
        for (int i = w; i < n; i++) {
            diff_fa[i - w + 1] = (0.5 * (ts[i] - ts[i - w]));
            diff_ga[i - w + 1] = (ts[i] - mua[i - w + 1]) + (ts[i - w] - mua[i - w]);
        }

        diff_fb[0] = 0;
        diff_gb[0] = 0;
        for (int i = w; i < qn; i++) {
            diff_fb[i - w + 1] = (0.5 * (query[i] - query[i - w]));
            diff_gb[i - w + 1] = (query[i] - mub[i - w + 1]) + (query[i - w] - mub[i - w]);
        }

        // AB JOIN
        for (int i = 0; i < profile_len; i++) {
            int mx = Math.min((profile_len - i), profile_lenb);

            cov_ = 0;
            for (int j = i; j < i + w; j++) {
                cov_ += ((ts[j] - mua[i]) * (query[j - i] - mub[0]));
            }

            for (int j = 0; j < mx; j++) {
                int k = j + i;
                cov_ += diff_fa[k] * diff_gb[j] + diff_ga[k] * diff_fb[j];
                corr_ = cov_ * siga[k] * sigb[j];

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

        // BA JOIN
        for (int i = 0; i < profile_lenb; i++) {
            int mx = Math.min((profile_lenb - i), profile_len);

            cov_ = 0;
            for (int j = i; j < i + w; j++) {
                cov_ += ((query[j] - mub[i]) * (ts[j - i] - mua[0]));
            }

            for (int j = 0; j < mx; j++) {
                int k = j + i;
                cov_ += diff_fb[k] * diff_ga[j] + diff_gb[k] * diff_fa[j];
                corr_ = cov_ * sigb[k] * siga[j];

                if (corr_ > mpb[k]) {
                    mpb[k] = corr_;
                    mpib[k] = j;
                }

                if (corr_ > mp[j]) {
                    mp[j] = corr_;
                    mpi[j] = k;
                }
            }
        }

        if (!crossCorrelation) {
            for (int i = 0; i < profile_len; i++) {
                if (mp[i] == -1.0) {
                    mp[i] = Double.POSITIVE_INFINITY;
                } else {
                    mp[i] = Math.sqrt(2.0 * w * (1.0 - mp[i]));
                }
            }

            for (int i = 0; i < profile_lenb; i++) {
                if (mpb[i] == -1.0) {
                    mpb[i] = Double.POSITIVE_INFINITY;
                } else {
                    mpb[i] = Math.sqrt(2.0 * w * (1.0 - mpb[i]));
                }
            }
        } else {
            for (int i = 0; i < profile_len; i++) {
                if (mp[i] > 1.0) {
                    mp[i] = 1.0;
                }
            }

            for (int i = 0; i < profile_lenb; i++) {
                if (mpb[i] > 1.0) {
                    mpb[i] = 1.0;
                }
            }
        }

        return new MPX(mp, mpi, mpb, mpib);
    }

    static double[][] muinvn(double[] a, int w) {
        var n = a.length;
        var h = new double[n];
        var r = new double[n];
        var mu = new double[n - w + 1];
        var sig = new double[n - w + 1];
        var p = a[0];
        double s = 0;
        double x, z, c, a1, a2, a3, mu_a;

        // compute moving mean
        for (var i = 1; i < w; i++) {
            x = p + a[i];
            z = x - p;
            s = s + ((p - (x - z)) + (a[i] - z));
            p = x;
        }

        mu[0] = (p + s) / w;
        for (var i = w; i < n; i++) {
            x = p - a[i - w + 1];
            z = x - p;
            s = s + ((p - (x - z)) - (a[i - w] + z));
            p = x;

            x = p + a[i];
            z = x - p;
            s = s + ((p - (x - z)) + (a[i] - z));
            p = x;

            mu[i - w + 1] = (p + s) / w;
        }

        // compute moving standard deviation
        for (var i = 0; i < n - w + 1; i++) {
            for (var j = i; j < i + w; j++) {
                mu_a = a[j] - mu[i];
                h[j] = mu_a * mu_a;

                c = (Math.pow(2, 27) + 1) * mu_a;
                a1 = (c - (c - mu_a));
                a2 = (mu_a - a1);
                a3 = a1 * a2;
                r[j] = a2 * a2 - (((h[j] - a1 * a1) - a3) - a3);
            }

            p = h[i];
            s = r[i];
            for (var j = i + 1; j < i + w; j++) {
                x = p + h[j];
                z = x - p;
                s = s + (((p - (x - z)) + (h[j] - z)) + r[j]);
                p = x;
            }

            sig[i] = p + s == 0 ? 0 : 1 / Math.sqrt(p + s);
        }

        return new double[][]{mu, sig};
    }
}
