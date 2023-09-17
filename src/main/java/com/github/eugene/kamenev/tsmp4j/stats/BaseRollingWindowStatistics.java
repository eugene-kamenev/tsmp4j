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

import static com.mvohm.quadruple.Quadruple.add;
import static com.mvohm.quadruple.Quadruple.divide;
import static com.mvohm.quadruple.Quadruple.multiply;
import static com.mvohm.quadruple.Quadruple.one;
import static com.mvohm.quadruple.Quadruple.subtract;

import com.github.eugene.kamenev.tsmp4j.utils.Buffer;
import com.github.eugene.kamenev.tsmp4j.utils.Buffer.DoubleBuffer;
import com.github.eugene.kamenev.tsmp4j.utils.Buffer.ObjBuffer;
import com.mvohm.quadruple.Quadruple;

/**
 * Class computes rolling window statistics for a data stream, which is used in Matrix Profile
 * algorithms. Can be used offline/online. K-algorithm taken from: <a
 * href="https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Computing_shifted_data">Computing
 * shifted data</a>
 */
public class BaseRollingWindowStatistics<S extends WindowStatistic>
    implements RollingWindowStatistics<S> {

    private final Buffer.DoubleBuffer dataBuffer;
    private final Buffer.ObjBuffer<S> statsBuffer;
    private final Quadruple n = new Quadruple();
    private final Quadruple K = new Quadruple();
    private final Quadruple Ex = new Quadruple();
    private final Quadruple Ex2 = new Quadruple();
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
    public S apply(double value) {
        totalDataCount++;
        if (Double.isNaN(value) || Double.isInfinite(value)) {
            toSkip = windowSize();
            value = 0.0d;
        } else {
            toSkip--;
        }
        if (this.dataBuffer.isFull()) {
            this.removeValue(this.dataBuffer.head());
        }
        this.dataBuffer.addToEnd(value);
        this.addValue(value);
        var mean = getMean();
        var populationVariance = getPopulationVariance();
        var stat = computeStats(value, mean, Math.sqrt(Math.max(0, populationVariance)),
            populationVariance, totalDataCount, toSkip > 0);
        getStatsBuffer().addToEnd(stat);
        return stat;
    }

    private void addValue(double value) {
        var val = new Quadruple(value);
        if (n.doubleValue() == 0.0d) {
            K.assign(val);
        }
        var diff = val.subtract(K);
        n.add(one());
        Ex.add(diff);
        Ex2.add(multiply(diff, diff));
    }

    private void removeValue(double value) {
        n.subtract(one());
        var diff = new Quadruple(value).subtract(K);
        Ex.subtract(diff);
        Ex2.subtract(multiply(diff, diff));
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
        return add(K, divide(Ex, n)).doubleValue();
    }

    private double getPopulationVariance() {
        return subtract(Ex2, multiply(Ex, Ex).divide(n)).divide(n).doubleValue();
    }
}
