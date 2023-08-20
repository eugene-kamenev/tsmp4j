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
import com.github.eugene.kamenev.tsmp4j.utils.Util;

/**
 * Class computes rolling window statistics for a data stream, which is used in Matrix Profile
 * algorithms. Can be used offline/online.
 */
public class BaseRollingWindowStatistics<S extends WindowStatistic>
    implements RollingWindowStatistics<S> {

    private final Buffer.DoubleBuffer dataBuffer;

    private final Buffer.ObjBuffer<S> statsBuffer;
    protected double currentMean = 0;
    protected double varianceSum = 0;
    private int dataCount = 0;
    private long totalDataCount = 0;
    private int toSkip = 0;

    public BaseRollingWindowStatistics(int windowSize, S[] statsBuffer) {
        this.dataBuffer = new DoubleBuffer(windowSize);
        this.statsBuffer = new ObjBuffer<>(statsBuffer);
    }

    @SuppressWarnings("unchecked")
    public BaseRollingWindowStatistics(int windowSize, int statsBufferSize) {
        this(windowSize, (S[]) new WindowStatistic[statsBufferSize]);
    }

    @Override
    public S apply(double dataPoint) {
        this.totalDataCount++;
        if (this.getDataBuffer().isFull()) {
            double oldestData = this.getDataBuffer().head();
            double updatedMean =
                (this.dataCount * this.currentMean - oldestData) / (this.dataCount - 1);
            this.varianceSum -= (oldestData - this.currentMean) * (oldestData - updatedMean);
            this.currentMean = updatedMean;
            this.dataCount--;
        }

        if (Double.isNaN(dataPoint) || Double.isInfinite(dataPoint)) {
            this.toSkip = this.windowSize();
            dataPoint = 0.0d;
        } else {
            this.toSkip--;
        }

        this.getDataBuffer().addToEnd(dataPoint);
        this.dataCount++;

        double previousMean = this.currentMean;
        this.currentMean += (dataPoint - this.currentMean) / this.dataCount;
        this.varianceSum += (dataPoint - previousMean) * (dataPoint - this.currentMean);

        double populationVariance = 0;
        if (this.dataCount > 1) {
            populationVariance = this.varianceSum / this.dataCount;
        }
        var stat = computeStats(dataPoint, this.currentMean, computeStdDev(populationVariance),
            totalDataCount, this.toSkip > 0);
        this.getStatsBuffer().addToEnd(stat);
        return stat;
    }

    @SuppressWarnings("unchecked")
    protected S computeStats(double x, double mean, double stdDev, long id, boolean skip) {
        return (S) new BaseWindowStatistic(x, mean, stdDev, id, skip);
    }

    /**
     * Computes the standard deviation from variance.
     *
     * @param variance The variance value.
     * @return The standard deviation.
     */
    protected double computeStdDev(double variance) {
        double stdDev = Math.sqrt(variance);
        return Util.sanitizeValue(stdDev);
    }

    @Override
    public DoubleBuffer getDataBuffer() {
        return this.dataBuffer;
    }

    @Override
    public ObjBuffer<S> getStatsBuffer() {
        return this.statsBuffer;
    }
}
