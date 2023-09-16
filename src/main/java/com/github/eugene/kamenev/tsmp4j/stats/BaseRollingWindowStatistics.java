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

package com.github.eugene.kamenev.tsmp4j.stats;

import com.github.eugene.kamenev.tsmp4j.utils.Buffer;
import com.github.eugene.kamenev.tsmp4j.utils.Buffer.DoubleBuffer;
import com.github.eugene.kamenev.tsmp4j.utils.Buffer.ObjBuffer;

/**
 * Class computes rolling window statistics for a data stream, which is used in Matrix Profile
 * algorithms. Can be used offline/online.
 */
public class BaseRollingWindowStatistics<S extends WindowStatistic>
    implements RollingWindowStatistics<S> {

    private final Buffer.DoubleBuffer dataBuffer;
    private final Buffer.ObjBuffer<S> statsBuffer;
    private final int w;       // Window size
    private double p = 0;      // Accumulated sum
    private double s = 0;      // Compensation for accumulated errors
    private double q = 0;      // Accumulated sum of squares
    private double s_q = 0;    // Compensation for accumulated errors in squared sum
    private long totalDataCount = 0;
    private int toSkip = 0;


    public BaseRollingWindowStatistics(int windowSize, S[] statsBuffer) {
        this.dataBuffer = new DoubleBuffer(windowSize);
        this.statsBuffer = new ObjBuffer<>(statsBuffer);
        this.w = windowSize;
    }

    @SuppressWarnings("unchecked")
    public BaseRollingWindowStatistics(int windowSize, int statsBufferSize) {
        this(windowSize, (S[]) new WindowStatistic[statsBufferSize]);
    }

    @Override
    public S apply(double value) {
        totalDataCount++;
        double mean = 0.0d;
        double populationVariance = 0.0d;

        if (Double.isNaN(value) || Double.isInfinite(value)) {
            toSkip = windowSize();
            value = 0.0d;
        } else {
            toSkip--;
        }

        if (!dataBuffer.isFull()) {
            addValue(value);
            if (dataBuffer.isFull()) {
                mean = getMean();
                populationVariance = getPopulationVariance();
            }
        } else {
            double oldestValue = dataBuffer.head();
            removeValue(oldestValue);
            addValue(value);
            mean = getMean();
            populationVariance = getPopulationVariance();
        }

        var stat = computeStats(value, mean, Math.sqrt(Math.max(0, populationVariance)),
            populationVariance, totalDataCount, toSkip > 0);
        getStatsBuffer().addToEnd(stat);
        return stat;
    }

    private void addValue(double value) {
        addToP(value);
        addToQ(value);
        dataBuffer.addToEnd(value);
    }

    private void removeValue(double value) {
        subtractFromP(value);
        subtractFromQ(value);
    }

    private void addToP(double value) {
        double y = value - s;
        double t = p + y;
        s = (t - p) - y;
        p = t;
    }

    private void subtractFromP(double value) {
        double y = -value - s;
        double t = p + y;
        s = (t - p) - y;
        p = t;
    }

    private void addToQ(double value) {
        double square = value * value;
        double y = square - s_q;
        double t = q + y;
        s_q = (t - q) - y;
        q = t;
    }

    private void subtractFromQ(double value) {
        double square = value * value;
        double y = -square - s_q;
        double t = q + y;
        s_q = (t - q) - y;
        q = t;
    }

    @SuppressWarnings("unchecked")
    protected S computeStats(double x, double mean, double stdDev, double variance, long id,
        boolean skip) {
        return (S) new BaseWindowStatistic(x, mean, stdDev, id, skip);
    }

    @Override
    public DoubleBuffer getDataBuffer() {
        return dataBuffer;
    }

    @Override
    public ObjBuffer<S> getStatsBuffer() {
        return statsBuffer;
    }

    private double getMean() {
        return (p + s) / w;
    }

    private double getPopulationVariance() {
        // not sure, maybe Welford's algorithm is better for variance
        double mean = getMean();
        double sqMean = mean * mean;
        return (q + s_q - w * sqMean) / w;
    }
}
