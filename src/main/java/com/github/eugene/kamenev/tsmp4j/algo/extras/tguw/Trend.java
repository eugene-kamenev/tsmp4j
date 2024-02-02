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

package com.github.eugene.kamenev.tsmp4j.algo.extras.tguw;

import static com.github.eugene.kamenev.tsmp4j.algo.extras.tguw.TGUWDecomposition.inverse;
import static com.github.eugene.kamenev.tsmp4j.algo.extras.tguw.Utils.diff;
import static com.github.eugene.kamenev.tsmp4j.algo.extras.tguw.Utils.longRunSd;
import static com.github.eugene.kamenev.tsmp4j.algo.extras.tguw.Utils.mad;

import java.util.List;

/**
 * trendsegmentR port from R, should support all the functions as original
 * https://cran.r-project.org/web/packages/trendsegmentR/
 * Performs the detection of linear trend changes for univariate time series by implementing
 * the bottom-up tail greedy unbalanced wavelet transformation proposed by H. Maeng and P. Fryzlewicz (2023).
 * https://arxiv.org/abs/1906.01939
 *
 * The estimated number and locations of the change-points are returned with the piecewise-linear estimator
 * for signal.
 *
 * @param changePoints list of change points, note that actual point position in x is point - 1
 * @param estimates estimated trend
 */
public record Trend(List<Integer> changePoints, double[] estimates) {

    /**
     * Defaults
     * @param x timeseries for segmentation
     * @return segmentation result
     */
    public static Trend segment(double[] x) {
        return segment(x, false, false, false, 0.04, -1, 0.0, -1, false);
    }

    public static Trend segment(double[] x, boolean independent, boolean continuous,
        boolean connected, double p, double thConst, double bal, int minsegL, boolean useKurtosis) {
        int n = x.length;
        if (n <= 1) {
            new Trend(List.of(), x);
        }
        if (minsegL <= 0) {
            minsegL = (int) Math.floor(0.9 * Math.log(n));
        }
        double sigma;
        if (independent) {
            sigma = mad(diff(diff(x))) / Math.sqrt(6);
        } else {
            sigma = longRunSd(x, 1.3, true);
        }
        var tguw = new TGUW(x, p);
        tguw.forward();
        var decomposition = tguw.getDecomposition();
        if (thConst <= 0) {
            thConst = decomposition.estimateThreshold(minsegL, -1, bal, false, false, useKurtosis);
        }
        double lambda = sigma * Math.sqrt(2 * Math.log(n)) * thConst;
        decomposition.threshold(lambda, minsegL, bal, connected);
        if (!connected) {
            var cpt = decomposition.findChangePoints(false, -1);
            double[] estimate = continuous ? Utils.estimateContinuous(cpt, x) : Utils.estimate(cpt, x);
            return new Trend(cpt, estimate);
        } else {
            inverse(decomposition.mergeHistory(), decomposition.transformed());
            var cpt = decomposition.findChangePoints(false, -1);
            if (continuous) {
                return new Trend(cpt, Utils.estimateContinuous(cpt, x));
            }
            return new Trend(cpt, decomposition.transformed());
        }
    }
}
