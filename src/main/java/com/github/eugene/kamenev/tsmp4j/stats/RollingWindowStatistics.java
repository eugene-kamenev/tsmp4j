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
import java.util.function.DoubleFunction;

/**
 * Class computes rolling window statistics for a data stream, which is used in Matrix Profile
 * algorithms. Can be used offline/online.
 */
public class RollingWindowStatistics implements DoubleFunction<Stats> {

    private final Buffer.DoubleBuffer dataBuffer;
    /**
     * Required for MPX algorithm.
     */
    private final boolean useInverseStdDev;
    private final Buffer.ObjBuffer<Stats> statsBuffer;
    private int dataCount = 0;
    private double currentMean = 0;
    private double varianceSum = 0;

    private long totalDataCount = 0;

    /**
     * Constructor for OnlineStatistics.
     *
     * @param windowSize The size of the data window for which statistics are computed.
     */
    public RollingWindowStatistics(int windowSize, int statsBufferSize,
        boolean useInverseStdDev) {
        this.dataBuffer = new Buffer.DoubleBuffer(windowSize);
        this.statsBuffer =
            statsBufferSize > 0 ? new Buffer.ObjBuffer<>(new Stats[statsBufferSize]) : null;
        this.useInverseStdDev = useInverseStdDev;
    }

    /**
     * Applies the online statistics algorithm on the new data point.
     *
     * @param dataPoint The new data point.
     * @return The computed statistics.
     */
    @Override
    public Stats apply(double dataPoint) {
        totalDataCount++;
        double oldestData = Double.NaN;
        if (this.dataBuffer.isFull()) {
            oldestData = this.dataBuffer.head();
            double updatedMean =
                (this.dataCount * this.currentMean - oldestData) / (this.dataCount - 1);
            this.varianceSum -= (oldestData - this.currentMean) * (oldestData - updatedMean);
            this.currentMean = updatedMean;
            this.dataCount--;
        }

        this.dataBuffer.addToEnd(dataPoint);
        this.dataCount++;

        double previousMean = this.currentMean;
        this.currentMean += (dataPoint - this.currentMean) / this.dataCount;
        this.varianceSum += (dataPoint - previousMean) * (dataPoint - this.currentMean);


        double sampleVariance = 0.0d;
        double populationVariance = 0;
        if (this.dataCount > 1) {
            sampleVariance = this.varianceSum / (this.dataCount - 1);
            populationVariance = this.varianceSum / this.dataCount;
        } else {
            sampleVariance = 0.0d;
            populationVariance = 0.0d;
        }
        var inverseStdDev = sanitizeValue(1 / Math.sqrt(varianceSum));
        var stats = new Stats(dataPoint, this.currentMean, this.varianceSum,
            computeStdDev(sampleVariance),
            computeStdDev(populationVariance), inverseStdDev, totalDataCount);
        if (this.statsBuffer != null) {
            this.statsBuffer.addToEnd(stats);
        }
        return stats;
    }

    /**
     * Computes the standard deviation from variance.
     *
     * @param variance The variance value.
     * @return The standard deviation.
     */
    private double computeStdDev(double variance) {
        double stdDev = Math.sqrt(variance);
        return sanitizeValue(stdDev);
    }

    /**
     * Sanitizes the value to ensure it's neither NaN nor Infinite.
     *
     * @param value The value to sanitize.
     * @return The sanitized value.
     */
    private double sanitizeValue(double value) {
        if (Double.isNaN(value) || Double.isInfinite(value)) {
            return 0.0;
        }
        return value;
    }

    /**
     * Checks if the buffer is full.
     *
     * @return True if the buffer is full, otherwise false.
     */
    public boolean isReady() {
        return this.dataBuffer.isFull();
    }

    /**
     * @return if the inverse standard deviation is used.
     */
    public boolean isUseInverseStdDev() {
        return useInverseStdDev;
    }

    public Buffer.ObjBuffer<Stats> getStatsBuffer() {
        return statsBuffer;
    }

    public double getMean(int i) {
        i += this.dataBuffer.getLength() - 1;
        return this.statsBuffer.get(i).meanX();
    }

    public double getStdDev(int i) {
        i += this.dataBuffer.getLength() - 1;
        var stat = this.statsBuffer.get(i);
        return this.useInverseStdDev ? stat.invStdDev() : stat.populationStdDev();
    }

    public double getX(int i) {
        return this.statsBuffer.get(i).x();
    }
}
