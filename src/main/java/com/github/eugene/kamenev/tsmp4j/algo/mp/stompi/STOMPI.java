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

package com.github.eugene.kamenev.tsmp4j.algo.mp.stompi;

import com.github.eugene.kamenev.tsmp4j.algo.mp.BaseMatrixProfile;
import com.github.eugene.kamenev.tsmp4j.algo.mp.DistanceProfileFunction.DistanceProfileQuery;
import com.github.eugene.kamenev.tsmp4j.algo.mp.MatrixProfile;
import com.github.eugene.kamenev.tsmp4j.algo.mp.MatrixProfileAlgorithm;
import com.github.eugene.kamenev.tsmp4j.algo.mp.mass.MASS2;
import com.github.eugene.kamenev.tsmp4j.algo.mp.stomp.STOMP;
import com.github.eugene.kamenev.tsmp4j.stats.BaseRollingWindowStatistics;
import com.github.eugene.kamenev.tsmp4j.stats.BaseWindowStatistic;
import com.github.eugene.kamenev.tsmp4j.stats.RollingWindowStatistics;
import com.github.eugene.kamenev.tsmp4j.utils.Util;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Real-time STOMP algorithm
 */
public class STOMPI implements
    MatrixProfileAlgorithm<BaseWindowStatistic, MatrixProfile> {

    private final RollingWindowStatistics<BaseWindowStatistic> rollingStatistics;

    private MatrixProfile matrixProfile;

    private final List<BaseWindowStatistic> history;

    private int newPoints = 0;

    private final int historySize;

    private final int exclusionZone;

    public STOMPI(BaseRollingWindowStatistics<BaseWindowStatistic> initialStats, int historySize) {
        this.rollingStatistics = new BaseRollingWindowStatistics<>(initialStats, 1);
        this.matrixProfile = new STOMP(initialStats).get();
        this.history = initialStats.getStatsBuffer()
            .toStream()
            .limit(initialStats.dataSize() - 1)
            .collect(Collectors.toList());
        this.historySize = historySize;
        this.exclusionZone = (int) Math.ceil(initialStats.windowSize() * 0.5d);
    }

    @Override
    public MatrixProfile get(RollingWindowStatistics<BaseWindowStatistic> query) {
        throw new UnsupportedOperationException("Not supported for STOMPI");
    }

    @Override
    public RollingWindowStatistics<BaseWindowStatistic> rollingStatistics() {
        return this.rollingStatistics;
    }

    @Override
    public void update(double value) {
        this.newPoints++;
        history.add(this.rollingStatistics().getStatsBuffer().get(0));
        MatrixProfileAlgorithm.super.update(value);
    }

    @Override
    public MatrixProfile get() {
        if (this.newPoints > 0) {
            int winSize = rollingStatistics().windowSize();
            var newBuffer = Stream.concat(
                this.history.stream(),
                rollingStatistics().getStatsBuffer().toStream()
            ).toArray(BaseWindowStatistic[]::new);
            var newStats = new BaseRollingWindowStatistics<>(winSize, newBuffer, true);
            var qIndex = newBuffer.length - rollingStatistics().windowSize() - this.newPoints + 1;
            var fft = Util.forwardFft(newStats, false, 0, Util.padSize(newBuffer.length));
            double[] firstProduct = null;
            double[] lastProduct = null;
            double[] distProfile = null;
            var newMatrixProfile = newMatrixProfile(this.matrixProfile);
            if (this.newPoints > 1) {
                var firstQuery = newQuery(newBuffer, 0);
                var firstDistanceProfile = new MASS2<BaseWindowStatistic>().apply(
                    new DistanceProfileQuery<>(newStats, firstQuery, 0, winSize, fft, false, false)
                );
                firstProduct = Arrays.stream(firstDistanceProfile.product())
                    .skip(winSize - 1)
                    .limit(newMatrixProfile.profile().length)
                    .toArray();
            }
            var dropValue = 0.0d;
            for (int i = 0; i < this.newPoints; i++) {
                var startIdx = qIndex + i;
                var query = newQuery(newBuffer, startIdx);
                if (i == 0) {
                    var distanceProfile = new MASS2<BaseWindowStatistic>().apply(
                        new DistanceProfileQuery<>(newStats, query, 0, winSize, fft, false, false)
                    );
                    distProfile = distanceProfile.profile();
                    lastProduct = Arrays.stream(distanceProfile.product())
                        .skip(winSize - 1)
                        .limit(distProfile.length)
                        .toArray();

                } else {
                    var cache = lastProduct[0];
                    for (int j = 1; j < newBuffer.length - winSize; j++) {
                        var temp = lastProduct[j];
                        lastProduct[j] = cache - newBuffer[j - 1].x() *
                            dropValue + newBuffer[j + winSize - 1].x() *
                            query.getStatsBuffer().get(winSize - 1).x();
                        distProfile[j] = computeDistance(j, winSize, lastProduct[j], newStats, query);
                        cache = temp;
                    }
                    lastProduct[0] = firstProduct[startIdx];
                    distProfile[0] = computeDistance(0, winSize, lastProduct[0], newStats, query);
                }
                dropValue = query.x(0);
                var excSt = Math.max(0, (qIndex + i) - exclusionZone);
                var min = Double.POSITIVE_INFINITY;
                var minIdx = -1;
                for (int j = 0; j < distProfile.length; j++) {
                    if (distProfile[j] < 0) {
                        distProfile[j] = 0;
                    }
                    if (j >= excSt) {
                        distProfile[j] = Double.POSITIVE_INFINITY;
                    } else {
                        distProfile[j] = Math.sqrt(distProfile[j]);
                    }
                    // update matrix profile
                    if (distProfile[j] < newMatrixProfile.profile()[j]) {
                        newMatrixProfile.indexes()[j] = startIdx;
                        newMatrixProfile.profile()[j] = distProfile[j];
                    }
                    // find min distance and its index
                    if (distProfile[j] < min) {
                        min = distProfile[j];
                        minIdx = j;
                    }
                    if (j >= startIdx) {
                        // update left profile
                        if (distProfile[j] < newMatrixProfile.leftProfile()[j]) {
                            newMatrixProfile.leftProfile()[j] = distProfile[j];
                            newMatrixProfile.leftIndexes()[j] = startIdx;
                        }
                    } else {
                        // update right profile
                        if (distProfile[j] < newMatrixProfile.rightProfile()[j]) {
                            newMatrixProfile.rightProfile()[j] = distProfile[j];
                            newMatrixProfile.rightIndexes()[j] = startIdx;
                        }
                    }
                }
                newMatrixProfile.indexes()[startIdx] = minIdx;
                newMatrixProfile.profile()[startIdx] = min;
                newMatrixProfile.leftIndexes()[startIdx] = minIdx;
                newMatrixProfile.leftProfile()[startIdx] = min;
            }
            if (this.historySize > 0) {
                int offset = 0;
                while(history.size() + this.rollingStatistics().getStatsBuffer().size() > this.historySize) {
                    history.remove(0);
                    offset++;
                }
                if (offset > 0) {
                    newMatrixProfile = newMatrixProfileOfRange(newMatrixProfile, offset);
                }
            }
            this.matrixProfile = newMatrixProfile;
        }
        this.newPoints = 0;
        return this.matrixProfile;
    }

    private double computeDistance(int i, int winSize, double lastProduct,
        RollingWindowStatistics<BaseWindowStatistic> newStats,
        RollingWindowStatistics<BaseWindowStatistic> query) {
        var q = query.getStatsBuffer().getR(0);
        var dist = (lastProduct - winSize * newStats.mean(i) * q.mean()) / (newStats.stdDev(i) * q.stdDev());
        return 2 * (winSize - dist);
    }

    private RollingWindowStatistics<BaseWindowStatistic> newQuery(BaseWindowStatistic[] stats, int from) {
        return new BaseRollingWindowStatistics<>(
            rollingStatistics().windowSize(),
            Stream.of(stats)
                .skip(from)
                .limit(rollingStatistics().windowSize())
                .toArray(BaseWindowStatistic[]::new),
            true
        );
    }

    private MatrixProfile newMatrixProfile(MatrixProfile profile) {
        var newSize = profile.profile().length + newPoints;
        double[] newProfile = Arrays.copyOf(profile.profile(), newSize);
        int[] newIndexes = Arrays.copyOf(profile.indexes(), newSize);
        double[] newRightProfile = Arrays.copyOf(profile.rightProfile(), newSize);
        double[] newLeftProfile = Arrays.copyOf(profile.leftProfile(), newSize);
        int[] newRightIndexes = Arrays.copyOf(profile.rightIndexes(), newSize);
        int[] newLeftIndexes = Arrays.copyOf(profile.leftIndexes(),newSize);
        for (int i = profile.profile().length; i < newProfile.length; i++) {
            newProfile[i] = Double.POSITIVE_INFINITY;
            newIndexes[i] = -1;
            newRightProfile[i] = Double.POSITIVE_INFINITY;
            newLeftProfile[i] = Double.POSITIVE_INFINITY;
            newRightIndexes[i] = -1;
            newLeftIndexes[i] = -1;
        }
        return new BaseMatrixProfile(
            newProfile,
            newIndexes,
            newRightProfile,
            newLeftProfile,
            newRightIndexes,
            newLeftIndexes
        );
    }

    private MatrixProfile newMatrixProfileOfRange(MatrixProfile profile, int offset) {
        var size = profile.profile().length;
        double[] newProfile = Arrays.copyOfRange(profile.profile(), offset, size);
        int[] newIndexes = Arrays.copyOfRange(profile.indexes(), offset, size);
        double[] newRightProfile = Arrays.copyOfRange(profile.rightProfile(), offset, size);
        double[] newLeftProfile = Arrays.copyOfRange(profile.leftProfile(), offset, size);
        int[] newRightIndexes = Arrays.copyOfRange(profile.rightIndexes(), offset, size);
        int[] newLeftIndexes = Arrays.copyOfRange(profile.leftIndexes(), offset, size);
        for (int i = 0; i < newIndexes.length; i++) {
            newIndexes[i] -= offset;
            newRightIndexes[i] -= offset;
            newLeftIndexes[i] -= offset;
        }
        return new BaseMatrixProfile(
            newProfile,
            newIndexes,
            newRightProfile,
            newLeftProfile,
            newRightIndexes,
            newLeftIndexes
        );
    }
}
