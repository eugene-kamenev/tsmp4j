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

import static com.github.eugene.kamenev.tsmp4j.MASS2.mpDist;

/**
 * STAMP: Streaming Time series Anytime Matrix Profile
 * Reference: Yeh CCM, Zhu Y, Ulanova L, Begum N, Ding Y, Dau HA, et al. Matrix profile I: All
 * pairs similarity joins for time series: A unifying view that includes motifs, discords and
 * shapelets. Proc - IEEE Int Conf Data Mining, ICDM. 2017;1317-22.
 * Reference: Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series Chains: A
 * New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1-27.
 * Reference Website: <a href="http://www.cs.ucr.edu/~eamonn/MatrixProfile.html">MatrixProfile.html</a>
 *
 * @param distances
 * @param indexes
 */
public record STAMP(double[] distances, int[] indexes) {

    public static STAMP stamp(double[] ts, int window) {
        return stamp(ts, window, ts.length - window + 1, ts, true);
    }

    public static STAMP stamp(double[] ts, double[] query, int window) {
        return stamp(query, window, query.length - window + 1, ts, false);
    }

    public static STAMP stamp(double[] tsA, int window, int w, double[] tsB, boolean trivialMatch) {
        var n = tsB.length;
        var matrixProfile = new double[n - window + 1];
        var matrixProfileIndex = new int[n - window + 1];
        Arrays.fill(matrixProfile, Double.POSITIVE_INFINITY);
        Arrays.fill(matrixProfileIndex, -1);

        double[] distanceProfile;
        int[] distanceProfileIndex;

        var index = 0;
        while (index < w) {
            distanceProfile = getDistanceProfile(tsA, tsB, index, window);
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

        return new STAMP(matrixProfile, matrixProfileIndex);
    }

    private static double[] getDistanceProfile(double[] tsA, double[] tsB, int index, int w) {
        var query = Arrays.copyOfRange(tsA, index, index + w);
        return mpDist(tsB, query);
    }

    private static int[] getDistanceProfileIndex(int n, int index, int w) {
        var distanceProfileIndex = new int[n - w + 1];
        Arrays.fill(distanceProfileIndex, index);
        return distanceProfileIndex;
    }
}
