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

import org.jtransforms.fft.DoubleFFT_1D;

import java.util.Arrays;

/**
 * MASS_V2
 * Mueen's Algorithm for Similarity Search is The Fastest Similarity Search Algorithm for Time
 * Series Subsequences under Euclidean Distance and Correlation Coefficient.
 * Reference: <a href="https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html">FastestSimilaritySearch.html</a>
 */
public class MASS2 {

    public static double[] mpDist(double[] ts, double[] query) {
        var n = ts.length;
        var m = query.length;
        var meany = Arrays.stream(query).average().orElse(Double.NaN);
        var sigmay = Math.sqrt(Arrays.stream(query).map(x -> Math.pow(x - meany, 2)).sum() / query.length);
        query = Util.reverse(query);
        var stats = movingAverageAndStd(ts, m);

        if(n - m > 0) {
            double[] paddedArray = new double[n];
            System.arraycopy(query, 0, paddedArray, 0, query.length);
            query = paddedArray;
        }

        var complexTs = new Complex[n];
        var complexQuery = new Complex[n];

        for (var i = 0; i < n; i++) {
            complexTs[i] = new Complex(ts[i], 0);
            complexQuery[i] = new Complex(query[i], 0);
        }

        complexTs = fft1D(complexTs);
        complexQuery = fft1D(complexQuery);


        //multiply two fft results
        var complexDot = new Complex[n];
        for (var i = 0; i < n; i++) {
            complexDot[i] = complexTs[i].mult(complexQuery[i]);
        }

        //inverse fft for dot computation
        complexDot = ifft1D(complexDot);
        var realDot = new double[complexDot.length];

        for (var i = 0; i < complexDot.length; i++) {
            realDot[i] = complexDot[i].getReal();
        }

        var dist = new double[n - m + 1];
        var shift = m - 1;
        for (var i = shift; i < n; i++) {
            var meanx = stats[0][i - shift];
            var sigmax = stats[1][i - shift];
            dist[i - m + 1] = Math.sqrt(2 * (m - (realDot[i] - m * meanx * meany) / (sigmax * sigmay)));
        }

        return dist;
    }

    private static Complex[] fft1D(Complex[] signal) {
        var n = signal.length;
        var fourier = new Complex[n];

        var coeff = new double[2 * n];
        var i = 0;
        for (var c : signal) {
            coeff[i++] = c.getReal();
            coeff[i++] = c.getImaginary();
        }

        var fft = new DoubleFFT_1D(n);
        fft.complexForward(coeff);

        for (i = 0; i < 2 * n; i += 2) {
            var c = new Complex(coeff[i], coeff[i + 1]);
            fourier[i / 2] = c;
        }
        return fourier;
    }

    private static Complex[] ifft1D(Complex[] fourier) {
        var n = fourier.length;
        var s = (1.0 / (double) n);

        var signal = new Complex[n];
        var coeff = new double[2 * n];

        var i = 0;
        for (Complex c : fourier) {
            coeff[i++] = c.getReal();
            coeff[i++] = c.getImaginary();
        }

        var fft = new DoubleFFT_1D(n);
        fft.complexInverse(coeff, false);

        for (i = 0; i < 2 * n; i += 2) {
            var c = new Complex(s * coeff[i], s * coeff[i + 1]);
            signal[i / 2] = c;
        }
        return signal;
    }

    static double[][] movingAverageAndStd(double[] data, int windowSize) {
        var result = new double[2][data.length - windowSize + 1];

        for (int i = 0; i < data.length - windowSize + 1; i++) {
            double sum = 0.0, sum2 = 0.0;
            for (var j = 0; j < windowSize; j++) {
                sum += data[i + j];
                sum2 += Math.pow(data[i + j], 2);
            }
            result[0][i] = sum / windowSize;
            var var = (sum2 / windowSize) - Math.pow(result[0][i], 2);
            result[1][i] = Math.sqrt(Math.max(var, 0.0));
        }

        return result;
    }
}
