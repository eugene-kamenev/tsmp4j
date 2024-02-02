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

import java.util.List;

/**
 * An array of dimension 4 by 3 by n-2 which has the full record of the n-2 merges in the TGUW
 * transformation. Each matrix contains the information of each merge. The first row shows the
 * indices of merged smooth coefficients in increasing order. Second row gives the value of detail
 * filter coefficients which is the weight vector for computing the corresponding detail
 * coefficient. The third row shows the (detail coefficient, first smooth coefficient, second smooth
 * coefficient) obtained by an orthonormal transform. The fourth row gives the balancedness of
 * merging. If it is Type 1 merging (three initial smooth coefficients) then the fourth row is
 * always (1/3, 1/3, 1/3). In Type 2 and Type 3 merging, the values depend on the ratio of the
 * length of the left and right wings to the entire merged region and only first two components of
 * the fourth row are filled with the corresponding ratios (sum to 1) but the third one is left as
 * NA.
 *
 * @param indexes indexes of merged smooth coefficients in increasing order
 * @param details value of detail filter coefficients (high pass coefficients)
 * @param smooths (detail coefficient, first smooth coefficient, second smooth coefficient)
 * @param balanced balancedness of merging
 */
public record MergeHistory(int[][] indexes, double[][] details, double[][] smooths, double[][] balanced) {

    public MergeHistory(int size) {
        this(new int[3][size], new double[3][size], new double[3][size], new double[3][size]);
    }

    public void update(int start, int end,
        List<int[]> edgesSubset, double[][] h, double[][] tc, double[][] balanced, int step) {
        for (int i = start, j = 0; i < end; i+=step) {
            var e = edgesSubset.get(j);
            for (int k = 0; k < 3; k++) {
                update(i, k, e[k], h[j][k], tc[j][k], balanced[k][j]);
            }
            j++;
        }
    }

    private void update(int index, int level, int edgeIndex, double detail, double smooth, double balance) {
        indexes[level][index] = edgeIndex;
        details[level][index] = detail;
        smooths[level][index] = smooth;
        balanced[level][index] = balance;
    }

    public MergeHistory copy() {
        int[][] indexesCopy = new int[3][];
        double[][] detailsCopy = new double[3][];
        double[][] smoothsCopy = new double[3][];
        double[][] balancedCopy = new double[3][];
        for (int i = 0; i < 3; i++) {
            indexesCopy[i] = indexes[i].clone();
            detailsCopy[i] = details[i].clone();
            smoothsCopy[i] = smooths[i].clone();
            balancedCopy[i] = balanced[i].clone();
        }
        return new MergeHistory(indexesCopy, detailsCopy, smoothsCopy, balancedCopy);
    }
}
