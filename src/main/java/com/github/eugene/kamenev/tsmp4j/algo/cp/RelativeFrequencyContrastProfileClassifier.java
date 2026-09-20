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
 * Time series classifier built on the Relative Frequency Contrast Profile primitive.
 * <p>
 * It follows the same nearest-plato decision rule as {@link ContrastProfileClassifier}, but each
 * class' representative subsequence is extracted with
 * {@link RelativeFrequencyContrastProfileAlgorithm}, which accounts for how frequently a pattern
 * recurs (relative frequency) instead of only its closest match. The positive series is the class'
 * own training series and the negative series is the concatenation of all other classes.
 * </p>
 */
public final class RelativeFrequencyContrastProfileClassifier {

    private final int windowSize;

    private final int maxFreq;

    private final double[][] platos;

    /**
     * @param windowSize subsequence length used both to extract each class' plato and to compare
     * query subsequences
     * @param maxFreq number of nearest neighbours used when building the relative frequency matrix
     * profile
     * @param classTrainingSeries one training series per class; each must contain at least two
     * instances of the behaviour that defines the class
     */
    public RelativeFrequencyContrastProfileClassifier(int windowSize, int maxFreq,
        double[][] classTrainingSeries) {
        if (classTrainingSeries.length < 2) {
            throw new IllegalArgumentException("At least two classes are required");
        }
        if (windowSize <= 0) {
            throw new IllegalArgumentException("windowSize must be positive");
        }
        if (maxFreq <= 0) {
            throw new IllegalArgumentException("maxFreq must be positive");
        }
        for (double[] series : classTrainingSeries) {
            if (series == null || series.length <= windowSize) {
                throw new IllegalArgumentException(
                    "Each training series must be non-null and longer than windowSize");
            }
        }
        this.windowSize = windowSize;
        this.maxFreq = maxFreq;
        this.platos = new double[classTrainingSeries.length][];
        var algorithm = new RelativeFrequencyContrastProfileAlgorithm(windowSize, maxFreq, false);
        for (int c = 0; c < classTrainingSeries.length; c++) {
            this.platos[c] = algorithm.apply(classTrainingSeries[c],
                PlatoNearestNeighbor.concatOthers(classTrainingSeries, c)).plato();
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

    public int maxFreq() {
        return maxFreq;
    }

    /**
     * @return a copy of the discriminative subsequence extracted for the given class
     */
    public double[] plato(int classIndex) {
        return Arrays.copyOf(platos[classIndex], platos[classIndex].length);
    }
}
