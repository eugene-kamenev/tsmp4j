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

import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import org.apache.commons.math3.stat.regression.SimpleRegression;

/**
 * The transformed x by the Tail Greedy Unbalanced Haar Wavelet transformation
 * with helper methods to operate on it
 */
public record TGUWDecomposition(
    MergeHistory mergeHistory,
    double[] transformed,
    double[] x,
    List<Integer> twotogether) {

    /**
     * Denoise the coefficients by thresholding
     * @param lambda threshold
     * @param minsegL minimum segment length
     * @param bal balance
     * @param connected is connected or not
     */
    public void threshold(double lambda, int minsegL, double bal, boolean connected) {
        var _protected = new HashMap<Integer, Boolean>();
        var indexes = mergeHistory.indexes();
        var coeffs = mergeHistory.smooths();
        var details = coeffs[0].clone();
        var balanced = mergeHistory.balanced();
        for (int i = 0; i < x.length - 2; i++) {
            if (!isProtected(indexes, i, _protected)) {
                double absCoeff = Math.abs(coeffs[0][i]);
                double bal1 = balanced[0][i];
                double bal2 = balanced[1][i];
                double bal3 = balanced[2][i];
                boolean condition = absCoeff > lambda && bal1 > bal && bal2 > bal;
                if (bal3 < 1) {
                    condition = condition && (1 >= minsegL);
                } else {
                    condition = condition && (Math.min(bal1 * bal3, bal2 * bal3) >= minsegL);
                }
                coeffs[0][i] *= condition ? 1 : 0;
            }

            if (Math.abs(coeffs[0][i]) > 0) {
                setProtected(indexes, i, _protected);
            }

        }

        int temp = -1;
        for (int c = 0, j = 0; c < twotogether.size(); c++) {
            var t = twotogether.get(c);
            if (t != 0) {
                j++;
                if (j == 2) {
                    j = 0;
                    var c1 = Math.abs(coeffs[0][temp]);
                    var c2 = Math.abs(coeffs[0][c]);
                    boolean isOverZero = c1 > 0;
                    boolean isOverZero2 = c2 > 0;
                    boolean isZero = c1 == 0;
                    boolean isZero2 = c2 == 0;
                    if ((isOverZero && isOverZero2) || (isZero && isZero2)) {
                        continue;
                    }
                    if (!isOverZero) {
                        coeffs[0][temp] = details[temp];
                    }
                    if (!isOverZero2) {
                        coeffs[0][c] = details[c];
                    }
                }
                temp = c;
            }
        }

        if (!connected) {
            balancedDenoise(minsegL, coeffs, balanced);
        }
    }

    /**
     * Find the change points
     * @param ordered
     * @param maxChangePoints
     * @return
     */
    public List<Integer> findChangePoints(boolean ordered, int maxChangePoints) {
        var twoTogether = getTwoTogether(ordered);
        var indexes = mergeHistory.indexes();
        var smooths = mergeHistory.smooths();
        int indexLen = mergeHistory.indexes()[0].length;
        var changePoints = new LinkedHashSet<Integer>();
        int[] survived = IntStream.range(0, indexLen)
            .filter(i -> Math.abs(smooths[0][i]) > Utils.EPS)
            .toArray();
        if (ordered) {
            int[] ttIndex = new int[survived.length];
            double[] ordDetails = new double[survived.length];
            for (int j = 0; j < survived.length; ) {
                ttIndex[j] = twoTogether.get(survived[j]);
                ordDetails[j] = Math.abs(smooths[0][survived[j]]);
                if (j < survived.length - 1 && ttIndex[j] == twoTogether.get(survived[j + 1])) {
                    ordDetails[j] = ordDetails[j + 1] = Math.max(ordDetails[j], Math.abs(smooths[0][survived[j + 1]]));
                    j += 2;
                } else {
                    j++;
                }
            }
            final var fSurvived = survived;
            survived = IntStream.range(0, ordDetails.length)
                .boxed()
                .sorted(Comparator.comparingDouble(i -> ordDetails[(int) i]).reversed())
                .mapToInt(ele -> fSurvived[ele])
                .toArray();
        }
        if (survived.length > 1) {
            boolean[] matched = new boolean[survived.length];
            for (int i = 0; i < indexLen; i++) {
                int j = 0;
                for (; j < survived.length; j++) {
                    var s = survived[j];
                    if (indexes[0][i] == indexes[1][s] || indexes[0][i] == indexes[2][s]) {
                        matched[j] = true;
                        break;
                    }
                }
            }
            int i = 0;
            while (i < survived.length - 1) {
                var s = survived[i];
                var nextS = survived[i + 1];
                var diff = indexes[1][s] - indexes[0][s];
                var nextDiff = indexes[1][nextS] - indexes[0][nextS];
                var diff2 = indexes[2][s] - indexes[1][s];
                var t = twoTogether.get(s);
                if (t != 0 && diff == 1 && nextDiff == 1) {
                    changePoints.add(indexes[0][s]);
                    changePoints.add(indexes[2][s]);
                    i += 2;
                } else if (t != 0 && diff2 == 1 && nextDiff == 1) {
                    changePoints.add(indexes[0][nextS]);
                    changePoints.add(indexes[2][nextS]);
                    i += 2;
                } else if (t == 0 && diff == 1 && diff2 != 1) {
                    changePoints.add(indexes[0][s]);
                    changePoints.add(indexes[2][s]);
                    i++;
                } else if (t == 0 && diff == 1 && diff2 == 1 && matched[i]) {
                    changePoints.add(indexes[0][s]);
                    changePoints.add(indexes[1][s]);
                    i++;
                } else {
                    changePoints.add(indexes[0][s]);
                    changePoints.add(indexes[1][s]);
                    changePoints.add(indexes[2][s]);
                    i++;
                }
            }
            changePoints.add(this.x.length + 1);
        } else if (survived.length == 1) {
            var s = survived[0];
            var diff = indexes[1][s] - indexes[0][s];
            var diff2 = indexes[2][s] - indexes[1][s];
            var t = twoTogether.get(s);
            if (t == 0 && diff == 1 && diff2 != 1) {
                changePoints.add(indexes[0][s]);
                changePoints.add(indexes[2][s]);
            } else if (t == 0 && diff == 1 && diff2 == 1) {
                changePoints.add(indexes[0][s]);
                changePoints.add(indexes[1][s]);
            }
        }
        Iterator<Integer> it = changePoints.iterator();
        // if the number of observation is three and the detail is survived then cp is the last
        if (this.x.length == 3 && !changePoints.isEmpty() && survived.length == 1) {
            it.remove();
            changePoints.add(this.x.length);
        } else {
            while (it.hasNext()) {
                int cp = it.next();
                if (cp <= 1 || cp >= this.x.length) {
                    it.remove();
                }
            }
        }
        if (maxChangePoints <= 0) {
            maxChangePoints = changePoints.size();
        }
        return changePoints.stream().map(i -> i - 1)
            .limit(maxChangePoints)
            .sorted()
            .collect(Collectors.toList());
    }

    private List<Integer> getTwoTogether(boolean ordered) {
        var twoTogether = twotogether;
        if (ordered) {
            twoTogether = new ArrayList<>(twoTogether.size());
            var ini = 1;
            var j = 0;
            while (j < twotogether.size()) {
                var prev = (int) twotogether.get(j);
                var curr = (int) twotogether.get(j + 1);
                if (prev != 0 && prev == curr) {
                    twoTogether.add(j, ini);
                    twoTogether.add(j + 1, ini);
                    ini++;
                    j += 2;
                } else {
                    twoTogether.add(j, prev);
                    j++;
                }
            }
        }
        return twoTogether;
    }

    /**
     * Estimate the threshold
     * @param minSegLen minimum segment length
     * @param maxChangePoints maximum change points
     * @param bal balance
     * @param continuous is continuous or not
     * @param connected is connected or not
     * @param useKurtosis use kurtosis or not
     * @return threshold
     */

    public double estimateThreshold(int minSegLen, int maxChangePoints, double bal,
        boolean continuous, boolean connected, boolean useKurtosis) {
        int n = x.length;
        if (minSegLen <= 0) {
            minSegLen = (int) Math.floor(0.9 * Math.log(n));
        }
        if (maxChangePoints <= 0) {
            maxChangePoints = (int) Math.ceil(0.1 * n);
        }
        var copy = this.copy();
        copy.threshold(1, minSegLen, bal, connected);
        var cpt = copy.findChangePoints(true, maxChangePoints);
        double[] estimate = !continuous ? Utils.estimate(cpt, x) : Utils.estimateContinuous(cpt, x);
        var diff = IntStream.range(0, estimate.length)
            .mapToDouble(i -> x[i] - estimate[i])
            .toArray();
        var kurt = 1.0;
        if (useKurtosis) {
            kurt = Utils.calculateKurtosis(diff);
        }
        // Fit an AR(1) model using simple linear regression
        SimpleRegression regression = new SimpleRegression();
        for (int i = 1; i < diff.length; i++) {
            regression.addData(diff[i - 1], diff[i]);
        }

        // Get the AR(1) coefficient (slope of the regression)
        double rhohat = regression.getSlope();

        // Calculate the long run standard deviation
        double longsd = Math.sqrt((1 + rhohat) / (1 - rhohat));
        return kurt * 1.3 * 1.3 * longsd;
    }

    /**
     * Inverse the Tail Greedy Unbalanced Haar Wavelet transformation
     * @param mergeHistory
     * @param transformed
     */
    public static void inverse(MergeHistory mergeHistory, double[] transformed) {
        double[] temp = new double[3];
        var smooths = mergeHistory.smooths();
        var indexes = mergeHistory.indexes();
        var details = mergeHistory.details();
        for (int i = transformed.length - 3; i >= 0; i--) {
            int _0 = indexes[0][i] - 1,
                _1 = indexes[1][i] - 1,
                _2 = indexes[2][i] - 1;
            temp[0] = details[0][i];
            temp[1] = details[1][i];
            temp[2] = details[2][i];
            // STEP 1: Obtain the orthonormal matrix of high pass filters
            var inv = Utils.orthmatrix(temp);
            // STEP 2: Choose the corresponding combination of detail and x value
            smooths[1][i] = transformed[_0];
            smooths[2][i] = transformed[_1];
            // STEP 3: Obtain 3 values of x from the combination of detail coefficient and x value
            var tmp1 = smooths[0][i];
            var tmp2 = smooths[1][i];
            var tmp3 = smooths[2][i];
            // Update with transformed values
            transformed[_0] = inv[0][0] * tmp1 + inv[1][0] * tmp2 + inv[2][0] * tmp3;
            transformed[_1] = inv[0][1] * tmp1 + inv[1][1] * tmp2 + inv[2][1] * tmp3;
            transformed[_2] = inv[0][2] * tmp1 + inv[1][2] * tmp2 + inv[2][2] * tmp3;
        }
    }

    private static boolean isProtected(int[][] indexes, int m,
        HashMap<Integer, Boolean> _protected) {
        return _protected.getOrDefault(indexes[0][m] - 1, false) ||
            _protected.getOrDefault(indexes[1][m] - 1, false) ||
            _protected.getOrDefault(indexes[2][m] - 1, false);
    }

    private static void setProtected(int[][] indexes, int m, HashMap<Integer, Boolean> _protected) {
        _protected.put(indexes[0][m] - 1, true);
        _protected.put(indexes[1][m] - 1, true);
    }

    private static void balancedDenoise(int minsegL, double[][] coeffs, double[][] balanced) {
        for (int i = 0; i < coeffs[0].length; i++) {
            if (coeffs[0][i] != 0) {
                if ((balanced[0][i] * balanced[2][i] < minsegL) ||
                    (balanced[1][i] * balanced[2][i]) < minsegL) {
                    coeffs[0][i] = 0;
                }
            }
        }
    }

    private TGUWDecomposition copy() {
        return new TGUWDecomposition(mergeHistory.copy(), transformed.clone(), this.x,
            new ArrayList<>(twotogether));
    }
}
