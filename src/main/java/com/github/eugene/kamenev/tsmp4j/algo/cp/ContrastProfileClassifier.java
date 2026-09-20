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

import com.github.eugene.kamenev.tsmp4j.algo.mp.mpx.MPXRollingWindowStatistics;
import java.util.Arrays;

/**
 * Time series classifier built on the Contrast Profile primitive (Matrix Profile XXIII).
 * <p>
 * For every class the classifier runs {@link ContrastProfileAlgorithm} using that class' training
 * series as the positive series and the concatenation of all other classes' training series as the
 * negative series. The resulting {@code plato} is the subsequence that is close to a repeated
 * pattern inside its own class yet far from every other class - a discriminative motif unique to
 * the class. A query is then classified by nearest plato (see {@link PlatoNearestNeighbor}).
 * </p>
 * This corresponds to the "preprocessing unstructured data for classification" downstream use
 * described in the Contrast Profile paper.
 */
public final class ContrastProfileClassifier {

    private final int windowSize;

    private final double[][] platos;

    /**
     * @param windowSize subsequence length used both to extract each class' plato and to compare
     * query subsequences
     * @param classTrainingSeries one training series per class; each must contain at least two
     * instances of the behaviour that defines the class
     */
    public ContrastProfileClassifier(int windowSize, double[][] classTrainingSeries) {
        if (classTrainingSeries.length < 2) {
            throw new IllegalArgumentException("At least two classes are required");
        }
        if (windowSize <= 0) {
            throw new IllegalArgumentException("windowSize must be positive");
        }
        for (double[] series : classTrainingSeries) {
            if (series == null || series.length <= windowSize) {
                throw new IllegalArgumentException(
                    "Each training series must be non-null and longer than windowSize");
            }
        }
        this.windowSize = windowSize;
        this.platos = new double[classTrainingSeries.length][];
        var algorithm = new ContrastProfileAlgorithm();
        for (int c = 0; c < classTrainingSeries.length; c++) {
            var positive = MPXRollingWindowStatistics.of(classTrainingSeries[c], windowSize);
            var negative = MPXRollingWindowStatistics.of(
                PlatoNearestNeighbor.concatOthers(classTrainingSeries, c), windowSize);
            this.platos[c] = algorithm.apply(positive, negative).plato();
        }
    }

    /**
     * @return the index of the class whose plato is closest to the query
     */
    public int classify(double[] query) {
        return PlatoNearestNeighbor.classify(query, platos);
    }

    /**
     * @return minimum distance from the query to every class plato
     */
    public double[] distances(double[] query) {
        return PlatoNearestNeighbor.distances(query, platos);
    }

    public int numClasses() {
        return platos.length;
    }

    public int windowSize() {
        return windowSize;
    }

    /**
     * @return a copy of the discriminative subsequence extracted for the given class
     */
    public double[] plato(int classIndex) {
        return Arrays.copyOf(platos[classIndex], platos[classIndex].length);
    }
}
