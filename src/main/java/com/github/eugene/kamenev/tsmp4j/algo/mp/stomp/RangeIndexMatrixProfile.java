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

package com.github.eugene.kamenev.tsmp4j.algo.mp.stomp;

import com.github.eugene.kamenev.tsmp4j.algo.mp.MatrixProfile;
import java.util.Arrays;

public record RangeIndexMatrixProfile(
    int windowSize,
    double exclusionZone,
    double[] profile,
    int[] indexes,
    double[] rightProfile,
    double[] leftProfile,
    int[] rightIndexes,
    int[] leftIndexes,
    RangeIndex[] subToLeft,
    RangeIndex[] subToRight,
    RangeIndex[] leftToRight,
    RangeIndex[] rightToLeft
) implements MatrixProfile {

    public static RangeIndex minSoFar(double[] dist, int startInd, int endInd) {
        int step = (endInd < startInd) ? -1 : 1; // Correct direction
        int length = dist.length;
        int chunksize = Math.max(1, (int) Math.ceil(Math.log(length)));

        double[] minSoFarVal = new double[chunksize];
        int[] minSoFarInd = new int[chunksize];

        double curMin = Double.POSITIVE_INFINITY;
        int curLen = 0;
        int curMaxLen = chunksize;

        for (int i = startInd; (step > 0 ? i <= endInd : i >= endInd); i += step) {
            if (dist[i] < curMin) {
                curMin = dist[i];
                curLen++;
                if (curMaxLen < curLen) {
                    curMaxLen += chunksize;
                    minSoFarVal = Arrays.copyOf(minSoFarVal, curMaxLen);
                    minSoFarInd = Arrays.copyOf(minSoFarInd, curMaxLen);
                }
                minSoFarVal[curLen - 1] = curMin;
                minSoFarInd[curLen - 1] = i;
            }
        }

        if (curLen > 0) {
            minSoFarVal = Arrays.copyOf(minSoFarVal, curLen);
            minSoFarInd = Arrays.copyOf(minSoFarInd, curLen);
            return new RangeIndex(minSoFarInd, minSoFarVal);
        } else {
            return null;
        }
    }

    public record RangeIndex(int[] indexes, double[] values) {

    }
}
