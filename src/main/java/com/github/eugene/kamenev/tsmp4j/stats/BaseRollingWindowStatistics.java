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
import java.math.BigDecimal;
import java.math.MathContext;

/**
 * Class computes rolling window statistics for a data stream, which is used in Matrix Profile
 * algorithms. Can be used offline/online.
 */
public class BaseRollingWindowStatistics<S extends WindowStatistic>
    implements RollingWindowStatistics<S> {

    private final Buffer.DoubleBuffer dataBuffer;

    private final Buffer.ObjBuffer<S> statsBuffer;

    /**
     * Usage of {@link MathContext#DECIMAL128} is required to avoid precision loss when computing
     * stats, but still, since we are not using BigDecimal for every computation, some errors are expected.
     */
    private final MathContext mathContext = MathContext.DECIMAL128;
    private long dataCount = 0;
    private long totalDataCount = 0;
    private int toSkip = 0;
    private BigDecimal sum = BigDecimal.ZERO;
    private BigDecimal sumSq = BigDecimal.ZERO;

    protected double varianceSum = 0.0d;

    private double currentMean = 0.0d;

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

        if (Double.isNaN(dataPoint) || Double.isInfinite(dataPoint)) {
            this.toSkip = this.windowSize();
            dataPoint = 0.0d;
        } else {
            this.toSkip--;
        }

        if (this.getDataBuffer().isFull()) {
            var oldestData = this.getDataBuffer().head();
            var x = BigDecimal.valueOf(oldestData);
            this.sum = this.sum.subtract(x, mathContext);
            this.sumSq = this.sumSq.subtract(x.multiply(x, mathContext));
            double updatedMean =
                (this.dataCount * this.currentMean - oldestData) / (this.dataCount - 1);
            this.varianceSum -= (oldestData - this.currentMean) * (oldestData - updatedMean);
            this.currentMean = updatedMean;
            this.dataCount--;
        }
        this.getDataBuffer().addToEnd(dataPoint);
        var x = BigDecimal.valueOf(dataPoint);
        this.dataCount++;
        var count = BigDecimal.valueOf(dataCount);
        this.sum = this.sum.add(x);
        this.sumSq = this.sumSq.add(x.multiply(x, mathContext));
        var mean = this.sum.divide(count, mathContext);
        double previousMean = this.currentMean;
        this.currentMean = mean.doubleValue();
        this.varianceSum += (dataPoint - previousMean) * (dataPoint - this.currentMean);
        var populationVariance = BigDecimal.ZERO;
        var stdDev = BigDecimal.ZERO;
        if (this.dataCount > 1) {
            var meanSq = mean.multiply(mean, mathContext);
            populationVariance = sumSq.divide(count, mathContext)
                .subtract(meanSq, mathContext);
            stdDev = populationVariance.compareTo(BigDecimal.ZERO) >= 0 ? populationVariance.sqrt(
                mathContext) : BigDecimal.ZERO;
        }
        var stat = computeStats(dataPoint, mean.doubleValue(), stdDev.doubleValue(),
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
