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

import java.util.ArrayList;
import java.util.List;

/**
 * Scalable KInetoscopic Matrix Profile (SKIMP) algorithm
 * Reference: <a href="https://sites.google.com/view/pan-matrix-profile/home">SKIMP</a>
 *
 * @param distances distance profiles
 * @param indexes index profiles
 */
public record SKIMP(double[][] distances, int[][] indexes) {

    public static class IndexRange {
        public int lowerBound;
        public int upperBound;

        public IndexRange(int lowerBound, int upperBound) {
            this.lowerBound = lowerBound;
            this.upperBound = upperBound;
        }
    }

    public static List<Integer> binarySplit(int n) {
        List<Integer> index = new ArrayList<>();
        List<IndexRange> intervals = new ArrayList<>();
        // Always begin by exploring the first integer
        index.add(0);

        if (n < 2) {
            return index;
        }

        // After exploring the first integer, split interval 2:n
        intervals.add(new IndexRange(1, n - 1));

        while (!intervals.isEmpty()) {
            IndexRange interval = intervals.remove(0);
            int lowerBound = interval.lowerBound;
            int upperBound = interval.upperBound;
            int middle = (lowerBound + upperBound) / 2;
            index.add(middle);

            if (lowerBound != upperBound) {
                IndexRange L = split(lowerBound, upperBound, middle).L;
                IndexRange R = split(lowerBound, upperBound, middle).R;

                if (L != null) {
                    intervals.add(L);
                }

                if (R != null) {
                    intervals.add(R);
                }
            }
        }

        return index;
    }

    public static class SplitResult {
        public IndexRange L;
        public IndexRange R;

        public SplitResult(IndexRange L, IndexRange R) {
            this.L = L;
            this.R = R;
        }
    }

    public static SplitResult split(int lowerBound, int upperBound, int middle) {
        if (lowerBound == middle) {
            return new SplitResult(null, new IndexRange(middle + 1, upperBound));
        } else if (upperBound == middle) {
            return new SplitResult(new IndexRange(lowerBound, middle - 1), null);
        } else {
            return new SplitResult(new IndexRange(lowerBound, middle - 1), new IndexRange(middle + 1, upperBound));
        }
    }

    public static SKIMP skimp(double[] ts, int[] range, boolean crossCorrelation) {
        var splitIDX = binarySplit(range.length);

        double[][] pmp = new double[range.length][];
        int[][] pmpi = new int[range.length][];
        int[] idx = new int[range.length];
        var lastIndex = splitIDX.size();
        for (int i = 0; i < lastIndex; i++) {
            var splitIDXVal = splitIDX.get(i);
            int subsequenceLength = range[splitIDXVal];
            var pnp = MPX.mpx(ts, subsequenceLength, crossCorrelation);
            var dist = pnp.mp();
            var idxs = pnp.mpi();
            for (int j = 0; j < dist.length; j++) {
                if (pmp[splitIDXVal] == null) {
                    pmp[splitIDXVal] = new double[dist.length];
                    pmpi[splitIDXVal] = new int[idxs.length];
                }
                pmp[splitIDXVal][j] = dist[j];
                pmpi[splitIDXVal][j] = idxs[j];
            }
            int j = splitIDXVal;
            while (j < splitIDX.size() && idx[j] != j) {
                idx[j] = splitIDXVal;
                j++;
            }
        }

        return new SKIMP(pmp, pmpi);
    }
}
