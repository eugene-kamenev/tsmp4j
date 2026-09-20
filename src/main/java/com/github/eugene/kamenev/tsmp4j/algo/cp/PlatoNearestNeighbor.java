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

import java.util.Arrays;

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

    static double minDistance(double[] query, double[] pattern) {
        int m = pattern.length;
        if (query.length < m) {
            return euclidean(zNormalize(query), zNormalize(pattern));
        }
        double[] normalizedPattern = zNormalize(pattern);
        double min = Double.POSITIVE_INFINITY;
        for (int i = 0; i + m <= query.length; i++) {
            double d = euclidean(zNormalize(Arrays.copyOfRange(query, i, i + m)), normalizedPattern);
            if (d < min) {
                min = d;
            }
        }
        return min;
    }

    private static double[] zNormalize(double[] data) {
        double mean = Arrays.stream(data).average().orElse(0);
        double variance = Arrays.stream(data).map(x -> (x - mean) * (x - mean)).average().orElse(0);
        double std = Math.sqrt(variance);
        if (std == 0) {
            return new double[data.length];
        }
        return Arrays.stream(data).map(x -> (x - mean) / std).toArray();
    }

    private static double euclidean(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            double diff = a[i] - b[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
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
