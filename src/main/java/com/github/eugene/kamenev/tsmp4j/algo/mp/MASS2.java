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

import com.github.eugene.kamenev.tsmp4j.stats.Stats;
import java.util.Arrays;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;
import org.apache.commons.math3.transform.TransformType;

/**
 * MASS_V2 Mueen's Algorithm for Similarity Search is The Fastest Similarity Search Algorithm for
 * Time Series Subsequences under Euclidean Distance and Correlation Coefficient. Reference: <a
 * href="https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html">FastestSimilaritySearch.html</a>
 */
public class MASS2 implements DistanceProfileFunction {

    @Override
    public double[] apply(DistanceProfileQuery dsq) {
        var n = dsq.ts().getStatsBuffer().getLength();
        var m = dsq.windowSize();
        var meanB = dsq.query().getMean(dsq.queryIndex());
        var stdDevB = dsq.query().getStdDev(dsq.queryIndex());
        var query = dsq.query().getStatsBuffer()
            .toStreamReversed()
            .skip(dsq.query().getStatsBuffer().getLength() - (dsq.windowSize() + dsq.queryIndex()))
            .limit(dsq.windowSize())
            .mapToDouble(Stats::x)
            .toArray();
        var ts = dsq.ts().getStatsBuffer()
            .toStream()
            .mapToDouble(Stats::x)
            .toArray();
        int padSize = (int) Math.pow(2, Math.ceil(Math.log(n) / Math.log(2)));
        if (padSize > n) {
            ts = Arrays.copyOf(ts, padSize);
        }
        if (n - m > 0 || m < padSize) {
            double[] paddedArray = new double[padSize];
            System.arraycopy(query, 0, paddedArray, 0, query.length);
            query = paddedArray;
        }

        var transformer = new FastFourierTransformer(DftNormalization.STANDARD);
        var tsFft = transformer.transform(ts, TransformType.FORWARD);
        var queryFft = transformer.transform(query, TransformType.FORWARD);
        var prod = new org.apache.commons.math3.complex.Complex[queryFft.length];
        for (int i = 0; i < queryFft.length; i++) {
            prod[i] = queryFft[i].multiply(tsFft[i]);
        }
        var inv = transformer.transform(prod, TransformType.INVERSE);
        double[] z = new double[inv.length];
        for (int i = 0; i < inv.length; i++) {
            z[i] = inv[i].getReal();
        }

        var dist = new double[n - m + 1];
        var shift = m - 1;
        for (var i = shift; i < n; i++) {
            var meanA = dsq.ts().getMean(i - shift);
            var stdDevA = dsq.ts().getStdDev(i - shift);
            var d = 2 * (m - (z[i] - m * meanA * meanB) / (stdDevA * stdDevB));
            dist[i - m + 1] = d > 0 ? Math.sqrt(d) : 0;
        }

        return dist;
    }
}
