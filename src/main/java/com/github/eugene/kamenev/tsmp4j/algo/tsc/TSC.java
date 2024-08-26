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

package com.github.eugene.kamenev.tsmp4j.algo.tsc;

import com.github.eugene.kamenev.tsmp4j.algo.mp.stomp.RangeIndexMatrixProfile;
import com.github.eugene.kamenev.tsmp4j.algo.mp.stomp.STOMP;
import com.github.eugene.kamenev.tsmp4j.stats.BaseWindowStatistic;
import com.github.eugene.kamenev.tsmp4j.stats.RollingWindowStatistics;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.function.Supplier;
import java.util.stream.IntStream;

/**
 * Robust Time Series Chain Discovery with Incremental Nearest Neighbors
 * <a href="https://arxiv.org/abs/2211.02146">Paper</a>
 * <a href="https://sites.google.com/view/robust-time-series-chain-22">Code</a>
 * <p>
 * Time series motif discovery has been a fundamental task to identify meaningful repeated patterns
 * in time series. Recently, time series chains were introduced as an expansion of time series
 * motifs to identify the continuous evolving patterns in time series data. Informally, a time
 * series chain (TSC) is a temporally ordered set of time series subsequences, in which every
 * subsequence is similar to the one that precedes it, but the last and the first can be arbitrarily
 * dissimilar. TSCs are shown to be able to reveal latent continuous evolving trends in the time
 * series, and identify precursors of unusual events in complex systems. Despite its promising
 * interpretability, unfortunately, we have observed that existing TSC definitions lack the ability
 * to accurately cover the evolving part of a time series: the discovered chains can be easily cut
 * by noise and can include non-evolving patterns, making them impractical in real-world
 * applications. Inspired by a recent work that tracks how the nearest neighbor of a time series
 * subsequence changes over time, we introduce a new TSC definition which is much more robust to
 * noise in the data, in the sense that they can better locate the evolving patterns while excluding
 * the non-evolving ones.
 * </p>
 */
public class TSC implements Supplier<int[][]> {

    private final STOMP stomp;

    public TSC(STOMP stomp) {
        this.stomp = stomp;
    }

    public TSC(RollingWindowStatistics<BaseWindowStatistic> stats) {
        this(new STOMP(stats, true));
    }

    public TSC(RollingWindowStatistics<BaseWindowStatistic> stats, double exclusionZone) {
        this(new STOMP(stats, exclusionZone, true));
    }

    public TSC(int winSize, int bufferSize) {
        this(new STOMP(winSize, bufferSize, true));
    }

    public TSC(int winSize, int bufferSize, double exclusionZone) {
        this(new STOMP(winSize, bufferSize, exclusionZone, true));
    }

    public void update(double value) {
        stomp.update(value);
    }

    @Override
    public int[][] get() {
        var mp = stomp.get();
        if (mp == null) {
            return null;
        }

        int winSize = mp.windowSize();
        int len = stomp.rollingStatistics().getStatsBuffer().getLength();
        int[][] chains = new int[len - winSize + 1][];

        for (int i = 0; i < chains.length; i++) {
            chains[i] = computeChain(i, winSize, mp);
            Arrays.sort(chains[i]);
        }

        return chains;
    }

    public static BestScore bestScore(double[] ts, int[][] chains, int winSize) {
        int maxScore = Integer.MIN_VALUE;
        double maxCorx = Double.NEGATIVE_INFINITY;
        int loc = -1;

        for (int i = 0; i < chains.length; i++) {
            var chain = chains[i];
            if (chain == null || chain.length == 0) continue;

            double maxDist = Double.NEGATIVE_INFINITY;
            double sumCorx = 0;

            var firstZ = zscore(Arrays.copyOfRange(ts, chain[0], chain[0] + winSize));
            var lastZ = zscore(Arrays.copyOfRange(ts, chain[chain.length - 1], chain[chain.length - 1] + winSize));
            double dist = calculateDistance(firstZ, lastZ);

            for (int j = 0; j < chain.length - 1; j++) {
                var curZ = zscore(Arrays.copyOfRange(ts, chain[j], chain[j] + winSize));
                double[] nextZ = zscore(Arrays.copyOfRange(ts, chain[j + 1], chain[j + 1] + winSize));

                // Update the max distance
                double currentDist = calculateDistance(curZ, nextZ);
                maxDist = Math.max(maxDist, currentDist);

                // Accumulate dot product sum for correlation
                sumCorx += dotProduct(curZ, nextZ);
            }

            int score = (int) Math.round(dist / maxDist);

            // Check for the maximum score and update loc
            if (score > maxScore || (score == maxScore && sumCorx > maxCorx)) {
                maxScore = score;
                maxCorx = sumCorx;
                loc = i;
            }
        }

        return new BestScore(loc, maxScore, chains[loc]);
    }

    private int[] computeChain(int idx, int winSize, RangeIndexMatrixProfile mp) {
        var chain = new ArrayList<Integer>();

        while (mp.leftIndexes()[idx] > -1) {
            chain.add(idx);
            idx = mp.leftIndexes()[idx];
        }

        if (chain.size() <= 1) {
            return new int[0];
        }

        var acceptChain = new ArrayList<Integer>();
        acceptChain.add(chain.get(0));
        var crit = new ArrayList<Integer>();

        for (int i = 0; i < chain.size() - 1; i++) {
            int anchor = chain.get(i + 1);
            int node = chain.get(i);
            if (isBSF(mp, anchor, node, winSize).isBSF()) {
                acceptChain.add(anchor);
                crit.add(node);
            } else if (!crit.isEmpty() && isBSF(mp, anchor, crit.get(crit.size() - 1), winSize).isBSF()) {
                acceptChain.add(anchor);
            } else {
                break;
            }
        }

        return acceptChain.stream().mapToInt(Integer::intValue).toArray();
    }

    private static double[] zscore(double[] data) {
        double mean = Arrays.stream(data).average().orElse(0);
        double std = Math.sqrt(Arrays.stream(data).map(x -> Math.pow(x - mean, 2)).average().orElse(0));
        return Arrays.stream(data).map(x -> (x - mean) / std).toArray();
    }

    private static double calculateDistance(double[] a, double[] b) {
        return Math.sqrt(IntStream.range(0, a.length).mapToDouble(i -> Math.pow(a[i] - b[i], 2)).sum());
    }

    private static double dotProduct(double[] a, double[] b) {
        return IntStream.range(0, a.length).mapToDouble(i -> a[i] * b[i]).sum();
    }

    private IsBSF isBSF(RangeIndexMatrixProfile mp, int x, int y, int winLen) {
        var len = winLen * 0.5d;
        var rIdx = mp.rightToLeft()[x];
        int minDiff = Math.abs(y - x);
        for (int i = rIdx.indexes().length - 1; i >= 0; i--) {
            var idx = rIdx.indexes()[i];
            var diff = Math.abs(y - idx);
            if (diff < minDiff) {
                minDiff = diff;
            }
        }
        return new IsBSF(minDiff < len, minDiff);
    }

    public record BestScore(int index, double score, int[] predIdxOurs) {

    }

    private record IsBSF(boolean isBSF, int idx) {

    }
}
