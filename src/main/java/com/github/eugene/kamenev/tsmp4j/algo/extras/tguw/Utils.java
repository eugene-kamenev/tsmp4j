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
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.stream.IntStream;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.rank.Median;
import org.apache.commons.math3.stat.regression.OLSMultipleLinearRegression;
import org.apache.commons.math3.stat.regression.SimpleRegression;

public final class Utils {

    public static final double EPS = Math.sqrt(Math.ulp(1.0));

    private Utils() {
    }

    /**
     * This method estimates the values of a given data set based on the provided change points.
     *
     * @param changePoints A collection of integers representing the change points in the data.
     * @param y            The data set for which the values are to be estimated.
     * @return A double array representing the estimated values.
     */
    public static double[] estimate(Collection<Integer> changePoints, double[] y) {
        int n = y.length;
        var cp = new int[changePoints.size() * 2 + 2];
        cp[0] = 1;
        int i = 1;
        for (var c : changePoints) {
            cp[i++] = c;
            cp[i++] = c + 1;
        }
        cp[i] = n;
        int temp = -1;
        double[] estimate = new double[n];
        SimpleRegression regression = new SimpleRegression();
        for (i = 0; i < cp.length; i++) {
            if (i % 2 != 0) {
                var start = temp - 1;
                var end = cp[i];
                if (end - start == 1) {
                    estimate[start] = y[start];
                    continue;
                }
                for (int m = start; m < end; m++) {
                    regression.addData(m + 1, y[m]);
                }
                for (int m = start; m < end; m++) {
                    estimate[m] = regression.predict(m + 1);
                }
                regression.clear();
            }
            temp = cp[i];
        }
        return estimate;
    }

    /**
     * This method estimates the continuous change points in a given data set.
     *
     * @param changePoints A collection of integers representing the change points in the data.
     * @param y            The data set for which the change points are to be estimated.
     * @return A double array representing the estimated continuous change points.
     */
    public static double[] estimateContinuous(Collection<Integer> changePoints, double[] y) {
        int n = y.length;
        // Initialize an array to store the change points
        var cp = new double[changePoints.size() + 4];
        cp[0] = 1;
        cp[1] = 1;
        int i = 2;
        // Populate the array with the change points
        for (var c : changePoints) {
            cp[i++] = c;
        }
        cp[i++] = n;
        cp[i] = n;
        // Generate a sequence of integers from 1 to n
        var ts = IntStream.rangeClosed(1, n).mapToDouble(j -> j).toArray();
        // Calculate the B-spline basis for the change points
        var bspline = Utils.bsplinebasis(1, filterDuplicates(cp), ts);
        // Initialize a multiple linear regression model
        var regression = new OLSMultipleLinearRegression(Utils.EPS);
        regression.setNoIntercept(true);
        // Fit the model to the data
        regression.newSampleData(y, bspline);
        // Estimate the residuals of the model
        var residuals = regression.estimateResiduals();
        // Return the estimated continuous change points
        return IntStream.range(0, n).mapToDouble(j -> y[j] - residuals[j]).toArray();
    }

    public static double[][] orthmatrix(double[] d) {
        double[][] M = new double[3][];

        // STEP 1 : Normalisation of high-pass filter row
        M[0] = d;
        double sum = 0;
        for (double v : M[0]) {
            sum += v * v;
        }
        double sqrtSum = Math.sqrt(sum);
        for (int i = 0; i < M[0].length; i++) {
            M[0][i] /= sqrtSum;
        }
        double u = M[0][0];
        double v = M[0][1];
        double w = M[0][2];

        // STEP 2 : Choose the vector which gives zero for the inner product of the first one
        M[1] = new double[]{1 - u * u, -u * v, -u * w};
        M[2] = new double[]{0, -w, v};

        // STEP 3 : Normalisation
        for (int i = 1; i < 3; i++) {
            sum = 0;
            for (double v1 : M[i]) {
                sum += v1 * v1;
            }
            // Reuse the previously calculated sqrtSum
            sqrtSum = Math.sqrt(sum);
            for (int j = 0; j < M[i].length; j++) {
                M[i][j] /= sqrtSum;
            }
        }

        return M;
    }

    /**
     * Calculate the B-spline basis functions, based on code from
     * <a href="https://github.com/eigenmatt/octave-bspline/blob/master/bsplinebasis.m">bsplinebasis.m</a>
     * Evaluate basis for a spline of given degree and knots at given points
     * Each column corresponds to a control point and each row corresponds to one
     * parameter value such that, for a column vector of control points c, the spline can be evaluated.
     * @param degree
     * @param knots
     * @param ts
     * @return B-spline basis functions
     */
    public static double[][] bsplinebasis(int degree, double[] knots, double[] ts) {
        int order = degree + 1;
        int nknots = knots.length;
        int npoints = nknots - order;
        double[][] B = new double[ts.length][npoints];

        for (int i = 0; i < ts.length; i++) {
            double t = ts[i];
            double searcht;
            if (t == knots[nknots - 1]) {
                // special case: if t is last knot, treat it as if it were
                // in the previous span, otherwise it would not be in any
                int j = nknots - 2;
                while (j >= 0 && knots[j] == t) {
                    j--;
                }
                searcht = knots[Math.max(j, 0)];
            } else {
                searcht = t;
            }

            // calculate 1st order basis functions
            // 1 if in knot span, 0 if not
            double[] temp = new double[nknots - 1];
            for (int j = 0; j < nknots - 1; j++) {
                temp[j] = (searcht >= knots[j] && searcht < knots[j + 1]) ? 1.0 : 0.0;
            }

            for (int k = 2; k <= order; k++) {
                // recursively calculate next order basis functions
                // by a linear combination of temp[j] and temp[j+1]
                for (int j = 0; j < nknots - k; j++) {
                    double d = 0.0;
                    double e = 0.0;
                    if (temp[j] != 0.0) {
                        d = ((t - knots[j]) * temp[j]) / (knots[j + k - 1] - knots[j]);
                    }
                    if (temp[j + 1] != 0.0) {
                        e = ((knots[j + k] - t) * temp[j + 1]) / (knots[j + k] - knots[j + 1]);
                    }
                    temp[j] = d + e;
                }
            }

            // Assign temp values to B
            for (int j = 0; j < npoints; j++) {
                B[i][j] = temp[j];
            }
        }

        return B;
    }

    /**
     * Filters duplicates from an array of doubles to avoid
     * singularity in QR decomposition, the same happens in R
     * @param arr
     * @return array where duplicates which counts greater than 2 are removed
     */
    public static double[] filterDuplicates(double[] arr) {
        var list = new ArrayList<Double>(arr.length);
        var countMap = new HashMap<Integer, Integer>();
        for (double num : arr) {
            int n = (int) num;
            var count = countMap.getOrDefault(n, 0);
            if (count < 2) {
                list.add(num);
            }
            countMap.put(n, count + 1);
        }
        return list.stream().mapToDouble(x -> x).toArray();
    }


    public static double[][] balanceNp(Collection<Integer> paired, List<int[]> ee,
        List<Integer> idx,
        int noOfCurrentSteps, int n) {
        var blnc = new double[3][noOfCurrentSteps];

        for (int i = 0; i < ee.size(); i++) {
            var e = ee.get(i);
            var justTwo = paired.contains(e[1]);
            var firstTwo = justTwo && paired.contains(e[0]);
            var lastTwo = justTwo && paired.contains(e[2]);

            if (firstTwo || lastTwo) {
                var index = idx.indexOf(e[2]);
                int prtn = firstTwo ? e[2] - e[0] :
                    (idx.contains(e[2]) && index + 1 < idx.size() ? idx.get(index + 1) : n) - e[1];
                blnc[0][i] = firstTwo ? prtn / (double) (prtn + 1) : 1 / (double) (prtn + 1);
                blnc[1][i] = firstTwo ? 1 / (double) (prtn + 1) : prtn / (double) (prtn + 1);
                blnc[2][i] = prtn + 1;
            } else {
                blnc[0][i] = blnc[1][i] = blnc[2][i] = 1 / 3.0;
            }
        }

        return blnc;
    }

    public static double[][] balanceP(int[][] pr, List<int[]> eeP1, List<Integer> idx,
        List<int[]> eeP2, int n) {
        int m = pr[0].length;
        double[][] blnc = new double[3][m];

        for (int i = 0; i < m; i++) {
            var i_P1 = eeP1.get(i);
            int index = idx.indexOf(eeP2.get(i)[2]);
            int c1 = i_P1[2] - i_P1[0];
            int c2 = (index >= 0 && index + 1 < idx.size()) ? idx.get(index + 1) - i_P1[2]
                : n - i_P1[2] + 1;

            blnc[0][i] = (double) c1 / (c1 + c2);
            blnc[1][i] = (double) c2 / (c1 + c2);
            blnc[2][i] = c1 + c2;
        }

        return blnc;
    }

    /**
     * This method applies a filter to an array of doubles. The filter is applied in a forward
     * direction, but if the result contains NaN values and reverse filtering is allowed, the filter
     * is applied in a reverse direction.
     *
     * @param a            The array of doubles to which the filter is to be applied.
     * @param allowReverse A boolean indicating whether reverse filtering is allowed.
     * @return An array of doubles representing the filtered values.
     */
    public static double[] filter(double[] a, boolean allowReverse) {
        // Calculate the terms used in the filter
        double term1 = a[4] * a[0] - a[3] * a[1];
        double term2 = a[1] * a[5] - a[2] * a[4];
        double term3 = a[3] * a[2] - a[0] * a[5];
        double term1Pow = Math.pow(term1, 2);

        // Apply the filter
        double w = -Math.sqrt(term1Pow / (Math.pow(term2, 2) + Math.pow(term3, 2) + term1Pow));
        double u = w * term2 / term1;
        double v = w * term3 / term1;
        double[] df = {u, v, w};
        // If the result contains NaN values and reverse filtering is allowed, apply the filter in a reverse direction
        if ((Double.isNaN(u) || Double.isNaN(v) || Double.isNaN(w)) && allowReverse) {
            for (int i = 0, j = a.length - 1; i < j; i++, j--) {
                double temp = a[i];
                a[i] = a[j];
                a[j] = temp;
            }
            df = filter(a, false);
            df = new double[]{df[2], df[1], df[0]};
        }
        return df;
    }

    /**
     * This method calculates the long-run standard deviation of a given data set.
     *
     * @param x The data set for which the long-run standard deviation is to be calculated.
     * @param pr The power ratio used in the calculation. The number of segments (m) is determined
     *           by taking the floor of n^(1/pr), where n is the length of the data set.
     * @param robust If true, the median of the absolute differences between consecutive means of
     *               segments is used in the calculation. If false, the mean of these differences is used.
     * @return The long-run standard deviation of the data set.
     */
    public static double longRunSd(double[] x, double pr, boolean robust) {
        int n = x.length;
        int m = (int) Math.floor(Math.pow(n, 1 / pr));
        int kn = Math.floorDiv(n, m);

        double[] A = new double[m];
        for (int i = 0; i < m; i++) {
            double[] segment = new double[kn];
            System.arraycopy(x, kn * i, segment, 0, Math.min(n - kn * i, kn));
            A[i] = new Mean().evaluate(segment);
        }

        double[] diffA = new double[A.length - 1];
        for (int i = 0; i < A.length - 1; i++) {
            diffA[i] = Math.abs(A[i + 1] - A[i]);
        }

        double longRunVar;
        if (robust) {
            longRunVar = (kn / 2.0) * Math.pow(new Median().evaluate(diffA), 2);
        } else {
            longRunVar = (kn / 2.0) * Math.pow(new Mean().evaluate(diffA), 2);
        }

        return Math.sqrt(longRunVar);
    }

    /**
     * This method calculates the kurtosis of a given data set.
     *
     * @param data The data set for which the kurtosis is to be calculated.
     * @return The kurtosis of the data set.
     */
    public static double calculateKurtosis(double[] data) {
        double mean = calculateMean(data);
        double n = data.length;
        double sum1 = 0;
        double sum2 = 0;

        for (double v : data) {
            var diff = v - mean;
            sum1 += Math.pow((diff), 4);
            sum2 += Math.pow((diff), 2);
        }

        return (n * sum1) / Math.pow(sum2, 2);
    }

    /**
     * This method calculates the mean of an array of doubles.
     *
     * @param data The array of doubles for which the mean is to be calculated.
     * @return The mean of the array of doubles.
     */
    public static double calculateMean(double[] data) {
        double sum = 0.0;
        for (double a : data) {
            sum += a;
        }
        return sum / data.length;
    }

    // Function to calculate the differences between consecutive elements of a vector
    public static double[] diff(double[] x) {
        double[] result = new double[x.length - 1];
        for (int i = 0; i < x.length - 1; i++) {
            result[i] = x[i + 1] - x[i];
        }
        return result;
    }

    // Function to calculate the Median Absolute Deviation (MAD) of a numeric vector
    // adds 1.4826 constant same as in R
    public static double mad(double[] x) {
        double median = new Median().evaluate(x);
        double[] absDeviations = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            absDeviations[i] = Math.abs(x[i] - median);
        }
        return 1.4826 * new Median().evaluate(absDeviations);
    }
}
