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

package com.github.eugene.kamenev.tsmp4j.algo.cp;

/**
 * Shared decision rule for the contrast-profile based time series classifiers.
 * <p>
 * Each class is represented by its {@code plato}: the most contrasting subsequence extracted from
 * that class' training series by a contrast profile algorithm. A query series is assigned to the
 * class whose plato has the smallest z-normalized Euclidean distance to any of the query's
 * subsequences. This mirrors the matrix-profile notion of similarity (z-normalized Euclidean
 * distance) that the contrast profiles themselves are built on.
 * </p>
 */
final class PlatoNearestNeighbor {

    /** Tolerance below which a subsequence is treated as constant and normalized to zeros. */
    private static final double STD_EPSILON = 1e-10d;

    private PlatoNearestNeighbor() {
    }

    /**
     * @return index of the closest class plato
     */
    static int classify(double[] query, double[][] platos) {
        double[] distances = distances(query, platos);
        int best = 0;
        for (int i = 1; i < distances.length; i++) {
            if (distances[i] < distances[best]) {
                best = i;
            }
        }
        return best;
    }

    /**
     * @return for every class plato the minimum z-normalized Euclidean distance between the plato
     * and any subsequence of the query
     */
    static double[] distances(double[] query, double[][] platos) {
        double[] distances = new double[platos.length];
        for (int c = 0; c < platos.length; c++) {
            distances[c] = minDistance(query, platos[c]);
        }
        return distances;
    }

    /**
     * @return the smallest z-normalized Euclidean distance between {@code pattern} and any
     * subsequence of {@code query}
     * @throws IllegalArgumentException if {@code query} is shorter than {@code pattern}
     */
    static double minDistance(double[] query, double[] pattern) {
        int m = pattern.length;
        if (query.length < m) {
            throw new IllegalArgumentException(
                "Query of length " + query.length + " is shorter than the plato of length " + m);
        }
        double[] normalizedPattern = zNormalize(pattern);
        // Z-normalization is shift invariant, so the scan carries its running statistics on
        // residuals around the middle of the query's value range. Every accumulated quantity then
        // stays small: neither the running variance nor the deviations lose precision when the
        // series sits on a large constant offset.
        double low = query[0];
        double high = query[0];
        for (double v : query) {
            low = Math.min(low, v);
            high = Math.max(high, v);
        }
        double reference = (low + high) / 2.0d;
        double residualSum = 0.0d;
        for (int i = 0; i < m; i++) {
            residualSum += query[i] - reference;
        }
        double meanResidual = residualSum / m;
        double sumSquaredDeviations = 0.0d;
        for (int i = 0; i < m; i++) {
            double deviation = (query[i] - reference) - meanResidual;
            sumSquaredDeviations += deviation * deviation;
        }
        double minSquared = Double.POSITIVE_INFINITY;
        for (int start = 0; start + m <= query.length; start++) {
            double std = Math.sqrt(Math.max(0.0d, sumSquaredDeviations / m));
            double squaredDistance = 0.0d;
            if (std < STD_EPSILON) {
                for (int k = 0; k < m; k++) {
                    squaredDistance += normalizedPattern[k] * normalizedPattern[k];
                }
            } else {
                double scale = 1.0d / std;
                for (int k = 0; k < m; k++) {
                    double deviation = (query[start + k] - reference) - meanResidual;
                    double diff = deviation * scale - normalizedPattern[k];
                    squaredDistance += diff * diff;
                }
            }
            if (squaredDistance < minSquared) {
                minSquared = squaredDistance;
            }
            if (start + m < query.length) {
                double outgoing = query[start] - reference;
                double incoming = query[start + m] - reference;
                double delta = incoming - outgoing;
                // Sliding update of the sum of squared deviations: ssd' = ssd
                //   + (in - out) * ((in - mean') + (out - mean)).
                double nextMeanResidual = meanResidual + delta / m;
                sumSquaredDeviations +=
                    delta * ((incoming - nextMeanResidual) + (outgoing - meanResidual));
                meanResidual = nextMeanResidual;
            }
        }
        return Math.sqrt(minSquared);
    }

    private static double[] zNormalize(double[] data) {
        if (data.length == 0) {
            return new double[data.length];
        }
        double low = data[0];
        double high = data[0];
        for (double v : data) {
            low = Math.min(low, v);
            high = Math.max(high, v);
        }
        double reference = (low + high) / 2.0d;
        double mean = 0.0d;
        for (double v : data) {
            mean += v - reference;
        }
        mean /= data.length;
        double variance = 0.0d;
        for (double v : data) {
            double deviation = (v - reference) - mean;
            variance += deviation * deviation;
        }
        variance /= data.length;
        double std = Math.sqrt(variance);
        if (std < STD_EPSILON) {
            return new double[data.length];
        }
        double[] normalized = new double[data.length];
        double scale = 1.0d / std;
        for (int i = 0; i < data.length; i++) {
            normalized[i] = ((data[i] - reference) - mean) * scale;
        }
        return normalized;
    }

    /**
     * @return concatenation of every training series except {@code exclude}, used as the negative
     * reference for the excluded class
     */
    static double[] concatOthers(double[][] series, int exclude) {
        int total = 0;
        for (int i = 0; i < series.length; i++) {
            if (i != exclude) {
                total += series[i].length;
            }
        }
        var result = new double[total];
        int pos = 0;
        for (int i = 0; i < series.length; i++) {
            if (i != exclude) {
                System.arraycopy(series[i], 0, result, pos, series[i].length);
                pos += series[i].length;
            }
        }
        return result;
    }
}
