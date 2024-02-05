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

package com.github.eugene.kamenev.tsmp4j.algo.extras.tguw;

import static com.github.eugene.kamenev.tsmp4j.algo.extras.tguw.Utils.balanceNp;
import static com.github.eugene.kamenev.tsmp4j.algo.extras.tguw.Utils.balanceP;
import static com.github.eugene.kamenev.tsmp4j.algo.extras.tguw.Utils.orthmatrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;
import java.util.function.Function;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * Performs the bottom-up unbalanced wavelet decomposition.
 * Details of the TGUW transformation can be found in H. Maeng and P. Fryzlewicz (2023),
 * Detecting linear trend changes in data sequences.
 */
public class TGUW {

    private final MergeHistory mergeHistory;

    private final List<int[]> edges;

    private final double[] weightsConst;

    private double[] details;

    private final double[] weightsLin;

    private final double[] tsCoeffs;

    private final int n;

    private int numberOfEdges;

    private int stepsLeft;

    private int currentStep;

    private final TreeSet<Integer> idx;

    private final TreeSet<Integer> paired;

    /**
     * A vector indicating locations of the detail coefficients returned by Type 3 merges
     * (merging two sets of paired smooth coefficients)
     */
    private final List<Integer> twoTogether;

    private final Map<Integer, Integer> removableNodes = new HashMap<>();

    private final TGUWDecomposition decomposition;

    private final double p;

    private final double[] temp = new double[6];

    private final List<double[][]> tempFilter = new ArrayList<>(2);

    /**
     * @param x An input vector to be decomposed
     * @param p Proportion of all possible remaining merges which specifies the number of
     *          merges allowed in a single pass over the data. The default is 0.04.
     */
    public TGUW(double[] x, double p) {
        this.n = x.length;
        this.p = p;
        this.numberOfEdges = n - 2;
        this.weightsConst = new double[n];
        Arrays.fill(weightsConst, 1);
        this.tsCoeffs = Arrays.copyOf(x, x.length);
        this.weightsLin = IntStream.rangeClosed(1, n).asDoubleStream().toArray();
        this.idx = new TreeSet<>(IntStream.rangeClosed(1, n).boxed().toList());
        this.paired = new TreeSet<>();
        this.edges = new ArrayList<>(numberOfEdges);
        for (int i = 0; i < numberOfEdges; i++) {
            edges.add(new int[] {i + 1, i + 2, i + 3, 0});
        }
        this.mergeHistory = new MergeHistory(numberOfEdges);
        this.stepsLeft = n - 1;
        this.currentStep = 0;
        this.twoTogether = new ArrayList<>(n);
        this.decomposition = new TGUWDecomposition(mergeHistory, tsCoeffs, x, twoTogether);
    }

    /**
     * Call the forward method and then get the decomposition to operate on it
     * using its helper methods
     * @return decomposition result
     */
    public TGUWDecomposition getDecomposition() {
        return decomposition;
    }

    /**
     * This method performs the forward transformation of the Tail Greedy Unbalanced Wavelet (TGUW).
     * It iteratively processes the edges until there are none left.
     */
    public void forward() {
        // Loop until all edges are processed
        var nonZero = new ArrayList<Integer>();
        while (!edges.isEmpty()) {
            // Calculate the maximum number of steps for the current iteration
            int maxCurrentSteps = (int) Math.ceil(p * stepsLeft);

            // Compute the details for the non-zero edges
            this.computeDetails(nonZero);
            nonZero.clear();

            // Sort the details in ascending order of their absolute values
            var ordDet = IntStream.range(0, details.length)
                .boxed()
                .sorted(Comparator.comparingDouble(i -> Math.abs(details[i])))
                .mapToInt(ele -> ele)
                .toArray();

            // Remove edges based on the maximum number of steps and the sorted details
            var detailsIndexes = removeEdges(maxCurrentSteps, numberOfEdges, edges, ordDet);

            // Get the number of current steps
            int noOfCurrentSteps = detailsIndexes.size();
            // Initialize a list to hold the subset of edges
            var edgesSubset = new ArrayList<int[]>(noOfCurrentSteps);
            // Merge the history of the edges
            this.mergeHistory(noOfCurrentSteps, detailsIndexes, edgesSubset);
            // Update the paired set with the first two elements of each edge in the subset
            // and remove the third element
            for (var edge : edgesSubset) {
                paired.add(edge[0]);
                paired.add(edge[1]);
                paired.remove(edge[2]);
            }
            // Clear the edges list for the next iteration
            edges.clear();
            int[] tempEdge = null;
            // Populate the edges list with the new edges
            int m = 1;
            // Initialize an iterator for the index set
            var iter = idx.iterator();
            for (int i = 0; i < idx.size() - 2; i++) {
                var e = new int[4];
                if (tempEdge == null) {
                    e[0] = iter.next();
                    e[1] = iter.next();
                    e[2] = iter.next();
                } else {
                    e[0] = tempEdge[1];
                    e[1] = tempEdge[2];
                    e[2] = iter.next();
                }
                var match = new BitSet(3);
                for (int j = 0; j < 3; j++) {
                    if (paired.contains(e[j])) {
                        match.set(j);
                    }
                }
                if (match.cardinality() == 3) {
                    // add to top of the list but after previous inserted
                    edges.add(m-1, e);
                    e[3] = (int) Math.ceil((double) m / 2);
                    nonZero.add(m);
                    m++;
                } else if (!(match.cardinality() == 1 || (match.get(0) && !match.get(1) && match.get(2)))) {
                    edges.add(e);
                }
                tempEdge = e;
            }

            // Update the numbers for the next iteration
            numberOfEdges = edges.size();
            currentStep += noOfCurrentSteps;
            stepsLeft -= noOfCurrentSteps;

            if (numberOfEdges == 1) {
                edges.get(0)[3] = 0;
            }
        }
    }

    /**
     * This method computes the details of the edges based on the provided indexes.
     * It separates the indexes into even and odd groups, and then applies a filter operation on the edges.
     * The filter operation calculates the weights and coefficients for each edge.
     * The results are stored in the provided arrays.
     * If an onCompute function is provided, it is applied on each edge after the filter operation.
     *
     * @param indexes The list of indexes to be processed.
     */
    private void computeDetails(List<Integer> indexes) {
        if (!indexes.isEmpty()) {
            // Separate the indexes into even and odd groups
            var even = IntStream.range(0, indexes.size())
                .filter(i -> i % 2 == 0)
                .map(indexes::get).toArray();
            var odd = IntStream.range(0, indexes.size())
                .filter(i -> i % 2 != 0)
                .map(indexes::get).toArray();
            // Initialize arrays to hold the weights and coefficients for each edge
            var evenWc = new double[even.length][3];
            var evenWl = new double[even.length][3];
            var evenTc = new double[even.length][3];
            var evenCoeffs = new double[even.length][];
            var wc = new double[even.length][2];
            var wl = new double[even.length][2];
            var tc = new double[even.length][2];
            var evenDetails = new double[even.length];
            var oddDetails = new double[odd.length];
            var pDetail = new double[even.length + odd.length];
            // Apply a filter operation on the edges
            filter(edges, evenDetails, even, evenWc, evenWl, evenTc, evenCoeffs, (i) -> {
                // Compute the orthogonal matrix for the coefficients
                var _m0 = orthmatrix(evenCoeffs[i]);
                // Update the weights and coefficients for each edge
                for (int k = 1; k < evenWc[i].length; k++) {
                    wc[i][k - 1] += evenWc[i][0] * _m0[k][0] + evenWc[i][1] * _m0[k][1]
                        + evenWc[i][2] * _m0[k][2];
                    wl[i][k - 1] += evenWl[i][0] * _m0[k][0] + evenWl[i][1] * _m0[k][1]
                        + evenWl[i][2] * _m0[k][2];
                    tc[i][k - 1] += evenTc[i][0] * _m0[k][0] + evenTc[i][1] * _m0[k][1]
                        + evenTc[i][2] * _m0[k][2];
                }
                // Get the odd edge and its index
                var oddEdgeIndex = odd[i] - 1;
                var oddEdge = edges.get(oddEdgeIndex);
                var _2 = oddEdge[2] - 1;

                // Update the temporary array with the weights and coefficients
                temp[0] = wc[i][0];
                temp[1] = wc[i][1];
                temp[2] = weightsConst[_2];
                temp[3] = wl[i][0];
                temp[4] = wl[i][1];
                temp[5] = weightsLin[_2];
                // Apply a filter operation on the temporary array
                double[] _h = Utils.filter(temp, true);
                // Update the details for the odd edges
                oddDetails[i] += _h[0] * tc[i][0] + _h[1] * tc[i][1] + _h[2] * tsCoeffs[_2];
                int pIdx = i == 0 ? 0 : i * 2;
                // Update the pDetail array with the maximum absolute value of the even and odd details
                pDetail[pIdx] = Math.max(Math.abs(evenDetails[i]), Math.abs(oddDetails[i]));
                pDetail[pIdx + 1] = pDetail[pIdx];
                return null;
            });
            // If the size of the indexes list is not equal to the size of the edges list, apply a filter operation on the remaining edges
            if (indexes.size() != edges.size()) {
                var edgeRow1 = IntStream.rangeClosed(1, edges.size())
                    .filter(i -> !indexes.contains(i)).toArray();
                var _details = new double[edgeRow1.length];
                filter(edges, _details, edgeRow1, null, null, null,  null, null);
                // combine the details
                details = Stream.of(Arrays.stream(pDetail), Arrays.stream(_details))
                    .flatMapToDouble(x1 -> x1)
                    .toArray();
            } else {
                details = pDetail;
            }
        } else {
            // If the indexes list is empty, apply a filter operation on all edges
            var details = new double[edges.size()];
            filter(edges, details,null, null, null, null,  null, null);
            this.details = details;
        }
    }

    /**
     * This method is responsible for merging edges. It first checks if there are any zero edges. If
     * there are, it adds them to the 'twoTogether' list. Then it updates the weights and
     * coefficients of the edges subset using the 'updating' method. Finally, it updates the merge
     * history using the 'mergeHistory.update' method.
     *
     * @param zeroEdges   A list of zero edges. If this is not null, the edges are added to the
     *                    'twoTogether' list.
     * @param from        The starting index for the update in the merge history.
     * @param to          The ending index for the update in the merge history.
     * @param edgesSubset The subset of edges to be updated.
     * @param idx0        A collection of indices representing the current index set.
     * @param balanced    A 2D array representing the balanced edges.
     * @param step        The current step in the merging process.
     */
    private void mergeEdges(List<Integer> zeroEdges, int from, int to, List<int[]> edgesSubset, Collection<Integer> idx0, double[][] balanced, int step) {
        if (zeroEdges != null) twoTogether.addAll(zeroEdges);
        this.updating(edgesSubset, weightsConst, weightsLin, tsCoeffs, idx0);
        this.mergeHistory.update(from, to, edgesSubset, this.tempFilter.get(0),
            this.tempFilter.get(1), balanced, step);
    }

    /**
     * This method is responsible for merging the history of the edges.
     * It first separates the edges into two categories: non-zero edges and zero edges.
     * If all edges are zero edges, it balances them and merges them.
     * If there are non-zero edges, it creates two sets of edges (ee_p1 and ee_p2) and balances and merges them.
     * If the size of the edges subset is not equal to the total, it balances and merges the remaining edges.
     *
     * @param noOfCurrentSteps The number of current steps.
     * @param detailsMinInd The indices of the details.
     * @param edgesSubset The subset of edges to be merged.
     */
    private void mergeHistory(int noOfCurrentSteps, List<Integer> detailsMinInd, List<int[]> edgesSubset) {
        // Create a copy of the current index set
        var idx0 = new ArrayList<>(idx);
        // Initialize lists to hold non-zero and zero edges
        List<Integer> nonZeroEdges = new ArrayList<>();
        List<Integer> zeroEdges = new ArrayList<>();
        // Separate the edges into non-zero and zero edges
        for (Integer index : detailsMinInd) {
            var e = edges.get(index);
            edgesSubset.add(e);
            if (e[3] != 0) {
                nonZeroEdges.add(e[3]);
            } else {
                zeroEdges.add(e[3]);
            }
        }
        // If all edges are zero edges, balance and merge them
        if (zeroEdges.size() == edgesSubset.size()) {
            var balanced = balanceNp(paired, edgesSubset, idx0, currentStep + noOfCurrentSteps, n);
            mergeEdges(zeroEdges, currentStep,  currentStep + noOfCurrentSteps, edgesSubset, idx,
                balanced, 1);
        } else {
            // If there are non-zero edges, create two sets of edges (ee_p1 and ee_p2)
            var pr = new int[2][nonZeroEdges.size() / 2];
            var ee_p1 = new ArrayList<int[]>(nonZeroEdges.size() / 2);
            var ee_p2 = new ArrayList<int[]>(nonZeroEdges.size() / 2);
            int total = 0, j = 0;
            for (int i = 0; i < edgesSubset.size(); i++) {
                if (edgesSubset.get(i)[3] != 0) {
                    int _0 = j % 2, _1 = j / 2;
                    pr[_0][_1] = i;
                    if (_0 != 0) {
                        edgesSubset.get(i)[0] = edgesSubset.get(i - 1)[0];
                        edgesSubset.get(i)[1] = edgesSubset.get(i - 1)[1];
                    } else {
                        total += 2;
                    }
                    (_0 != 0 ? ee_p2 : ee_p1).add(edgesSubset.get(i));
                    j++;
                }
            }
            // Balance and merge the two sets of edges
            var balanced = balanceP(pr, ee_p1, idx0, ee_p2, n);
            this.mergeEdges(nonZeroEdges, currentStep, currentStep + total, ee_p1, idx, balanced, 2);
            this.mergeEdges(null, currentStep + 1, currentStep + total + 1, ee_p2, idx, balanced,
                2);
            // If the size of the edges subset is not equal to the total, balance and merge the remaining edges
            if (edgesSubset.size() != total) {
                var indicies = edgesSubset.stream().filter(arr -> arr[3] == 0).mapToInt(arr -> arr[3]).boxed().toList();
                var edgesSubsetNp = edgesSubset.stream().filter(arr -> arr[3] == 0).toList();
                balanced = balanceNp(paired, edgesSubsetNp, idx0, noOfCurrentSteps - total, n);
                this.mergeEdges(indicies, currentStep + total, currentStep + noOfCurrentSteps, edgesSubsetNp, idx, balanced, 1);
            }
        }
    }

    /**
     * This method applies a filter operation on the given edges. It calculates the weights and coefficients
     * for each edge and applies a filter function on them. The results are stored in the provided arrays.
     * If an onCompute function is provided, it is applied on each edge after the filter operation.
     *
     * @param edges The list of edges to be filtered.
     * @param details The array to store the details of each edge after the filter operation. If null, a new array is created.
     * @param row The indices of the edges to be filtered. If null, all edges are filtered.
     * @param wc The array to store the weights of each edge after the filter operation. If null, a new array is created.
     * @param wl The array to store the linear weights of each edge after the filter operation. If null, a new array is created.
     * @param tc The array to store the coefficients of each edge after the filter operation. If null, a new array is created.
     * @param h The array to store the filter results of each edge. If null, a new array is created.
     * @param onCompute A function to be applied on each edge after the filter operation. If null, no function is applied.
     */
    private void filter(List<int[]> edges, double[] details, int[] row,
        double[][] wc, double[][] wl, double[][] tc, double[][] h, Function<Integer, Void> onCompute) {

        int len = row == null ? edges.size() : row.length;
        if (wc == null) {
            wc = new double[len][3];
        }
        if (wl == null) {
            wl = new double[len][3];
        }
        if (tc == null) {
            tc = new double[len][3];
        }
        if (h == null) {
            h = new double[len][];
        }
        for (int i = 0, j = 0; i < len; i++) {
            int edgeIndex = row == null ? j++ : row[i] - 1;
            int[] e = edges.get(edgeIndex);
            int _0 = e[0] - 1,
                _1 = e[1] - 1,
                _2 = e[2] - 1;
            tc[i][0] = tsCoeffs[_0];
            tc[i][1] = tsCoeffs[_1];
            tc[i][2] = tsCoeffs[_2];
            temp[0] = wc[i][0] = weightsConst[_0];
            temp[1] = wc[i][1] = weightsConst[_1];
            temp[2] = wc[i][2] = weightsConst[_2];
            temp[3] = wl[i][0] = weightsLin[_0];
            temp[4] = wl[i][1] = weightsLin[_1];
            temp[5] = wl[i][2] = weightsLin[_2];
            h[i] = Utils.filter(temp, true);
            details[i] += h[i][0] * tc[i][0] + h[i][1] * tc[i][1] + h[i][2] * tc[i][2];
            if (onCompute != null) {
                onCompute.apply(i);
            }
        }
    }

    /**
     * This method updates the weights and coefficients of the edges in the given list.
     * It first initializes arrays to hold the weights and coefficients of each edge.
     * Then it applies a filter operation on the edges, which calculates the weights and coefficients for each edge.
     * After the filter operation, it updates the weights and coefficients of the edges using the results of the filter operation.
     * Finally, it clears the tempFilter list and adds the new weights and coefficients to it.
     *
     * @param edges The list of edges to be updated.
     * @param weightsConst The array of constant weights for each edge.
     * @param weightsLin The array of linear weights for each edge.
     * @param tsCoeffs The array of coefficients for each edge.
     * @param idx A collection of indices representing the current index set.
     */
    public void updating(List<int[]> edges, double[] weightsConst, double[] weightsLin,
        double[] tsCoeffs, Collection<Integer> idx) {
        int len = edges.size();
        var wc = new double[len][3];
        var wl = new double[len][3];
        var tc = new double[len][3];
        var h = new double[len][];
        var details = new double[len];
        filter(edges, details, null, wc, wl, tc, h, (i) -> {
            var _m0 = orthmatrix(h[i]);
            var e = edges.get(i);
            int _0 = e[0] - 1, _1 = e[1] - 1, _2 = e[2] - 1;
            double _wc0_0 = wc[i][0], _wc0_1 = wc[i][1], _wc0_2 = wc[i][2],
                _wl0_0 = wl[i][0], _wl0_1 = wl[i][1], _wl0_2 = wl[i][2],
                _tc0_0 = tc[i][0], _tc0_1 = tc[i][1], _tc0_2 = tc[i][2];

            for (int k = 0; k < 3; k++) {
                var wcValue = _wc0_0 * _m0[k][0] + _wc0_1 * _m0[k][1] + _wc0_2 * _m0[k][2];
                var wlValue = _wl0_0 * _m0[k][0] + _wl0_1 * _m0[k][1] + _wl0_2 * _m0[k][2];
                var tcValue = _tc0_0 * _m0[k][0] + _tc0_1 * _m0[k][1] + _tc0_2 * _m0[k][2];
                if (k == 0) {
                    weightsConst[_2] = wc[i][k] = wcValue;
                    weightsLin[_2] = wl[i][k] = wlValue;
                    tsCoeffs[_2] = tc[i][k] = tcValue;
                } else if (k == 1) {
                    weightsConst[_0] = wc[i][k] = wcValue;
                    weightsLin[_0] = wl[i][k] = wlValue;
                    tsCoeffs[_0] = tc[i][k] = tcValue;
                } else {
                    weightsConst[_1] = wc[i][k] = wcValue;
                    weightsLin[_1] = wl[i][k] = wlValue;
                    tsCoeffs[_1] = tc[i][k] = tcValue;
                }
            }

            idx.remove(e[2]);
            return null;
        });
        this.tempFilter.clear();
        this.tempFilter.add(h);
        this.tempFilter.add(tc);
    }

    /**
     * This method is responsible to provide removable edge indexes.
     * It iterates over the edges and checks if the sum of the removable nodes for each edge is equal to 3 or 6.
     * If the condition is met, the edge is removed and its index is added to the list of removed indexes.
     * The method continues until the size of the removed indexes list is less than the maximum current steps and
     * the edge index is less than the number of edges.
     *
     * @param maxCurrentSteps The maximum number of steps that can be taken.
     * @param noe The number of edges.
     * @param edges The list of edges.
     * @param ordDet The ordered details.
     * @return A list of integers representing the indexes of the removed edges.
     */
    private List<Integer> removeEdges(int maxCurrentSteps, int noe, List<int[]> edges,
        int[] ordDet) {
        removableNodes.clear();
        List<Integer> removedIndexes = new ArrayList<>(maxCurrentSteps);
        int tei = 0;
        while (removedIndexes.size() < maxCurrentSteps && tei < noe) {
            var index = ordDet[tei];
            int[] edge = edges.get(index);
            boolean isEdgeRemoved = edge[3] > 0;
            int sum = sumNodes(edge);
            var next = tei + 1;
            var nextEdge = next < edges.size() ? edges.get(next) : null;
            if (isEdgeRemoved && nextEdge != null) {
                sum += sumNodes(nextEdge);
            }

            if (sum == (isEdgeRemoved ? 6 : 3)) {
                removedIndexes.add(index);
                removeNodes(edge);
                if (isEdgeRemoved && nextEdge != null) {
                    index = ordDet[next];
                    removeNodes(nextEdge);
                    removedIndexes.add(index);
                }
            }
            tei += isEdgeRemoved ? 2 : 1;
        }
        return removedIndexes;
    }

    /**
     * Helper function to sum the nodes of an edge.
     * @param e
     * @return
     */
    private int sumNodes(int[] e) {
        return removableNodes.getOrDefault(e[0] - 1, 1) +
            removableNodes.getOrDefault(e[1] - 1, 1) +
            removableNodes.getOrDefault(e[2] - 1, 1);
    }

    /**
     * Helper function to remove the nodes of an edge.
     * @param e
     */
    private void removeNodes(int[] e) {
        removableNodes.put(e[0] - 1, 0);
        removableNodes.put(e[1] - 1, 0);
        removableNodes.put(e[2] - 1, 0);
    }
}
