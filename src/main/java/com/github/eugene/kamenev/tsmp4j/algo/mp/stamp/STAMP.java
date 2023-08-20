package com.github.eugene.kamenev.tsmp4j.algo.mp.stamp;

import com.github.eugene.kamenev.tsmp4j.algo.mp.BaseMatrixProfile;
import com.github.eugene.kamenev.tsmp4j.algo.mp.BaseMatrixProfileAlgorithm;
import com.github.eugene.kamenev.tsmp4j.algo.mp.DistanceProfileFunction;
import com.github.eugene.kamenev.tsmp4j.algo.mp.DistanceProfileFunction.DistanceProfileQuery;
import com.github.eugene.kamenev.tsmp4j.algo.mp.MatrixProfile;
import com.github.eugene.kamenev.tsmp4j.algo.mp.mass.MASS2;
import com.github.eugene.kamenev.tsmp4j.stats.BaseRollingWindowStatistics;
import com.github.eugene.kamenev.tsmp4j.stats.BaseWindowStatistic;
import com.github.eugene.kamenev.tsmp4j.stats.RollingWindowStatistics;
import com.github.eugene.kamenev.tsmp4j.stats.WindowStatistic;
import com.github.eugene.kamenev.tsmp4j.utils.Util;
import java.util.Arrays;
import org.apache.commons.math3.complex.Complex;

/**
 * STAMP: Scalable Time Series Anytime Matrix Profile Reference: Yeh CCM, Zhu Y, Ulanova L, Begum N,
 * Ding Y, Dau HA, et al. Matrix profile I: All pairs similarity joins for time series: A unifying
 * view that includes motifs, discords and shapelets. Proc - IEEE Int Conf Data Mining, ICDM.
 * 2017;1317-22. Reference: Zhu Y, Imamura M, Nikovski D, Keogh E. Matrix Profile VII: Time Series
 * Chains: A New Primitive for Time Series Data Mining. Knowl Inf Syst. 2018 Jun 2;1-27. Reference
 * Website: <a href="http://www.cs.ucr.edu/~eamonn/MatrixProfile.html">MatrixProfile.html</a>
 */
public class STAMP extends BaseMatrixProfileAlgorithm<BaseWindowStatistic, MatrixProfile> {

    public STAMP(RollingWindowStatistics<BaseWindowStatistic> rollingWindowStatistics) {
        super(rollingWindowStatistics);
    }

    public STAMP(int windowSize, int bufferSize) {
        this(new BaseRollingWindowStatistics<>(windowSize, bufferSize));
    }

    @Override
    public MatrixProfile get(RollingWindowStatistics<BaseWindowStatistic> query) {
        if (this.isReady()) {
            return stamp(this.rollingStatistics(), query);
        }
        return null;
    }

    @Override
    public MatrixProfile get() {
        if (this.isReady()) {
            return stamp(this.rollingStatistics(), null);
        }
        return null;
    }

    public static <S extends WindowStatistic> MatrixProfile stamp(RollingWindowStatistics<S> ts,
        RollingWindowStatistics<S> query, double exclusionZone, int sSize,
        DistanceProfileFunction<S> distFunc) {
        int windowSize = ts.windowSize();
        boolean isJoin = query != null;
        if (!isJoin) {
            query = ts;
        } else {
            exclusionZone = 0;
        }
        int exZone = (int) Math.round(windowSize * exclusionZone + Util.EPS);
        int dataSize = ts.dataSize();
        int querySize = query.dataSize();
        int mpSize = dataSize - windowSize + 1;
        int numQueries = querySize - windowSize + 1;
        if (querySize > dataSize) {
            throw new IllegalArgumentException(
                "Query must be smaller or the same size as reference data.");
        }
        if (windowSize < 4) {
            throw new IllegalArgumentException("Window size must be at least 4.");
        }

        double[] leftMatrixProfile, rightMatrixProfile;
        int[] leftProfileIndex, rightProfileIndex;

        var matrixProfile = new double[mpSize];
        var profileIndex = new int[mpSize];

        if (isJoin) {
            leftMatrixProfile = rightMatrixProfile = null;
            leftProfileIndex = rightProfileIndex = null;
        } else {
            leftMatrixProfile = new double[matrixProfile.length];
            rightMatrixProfile = new double[matrixProfile.length];
            leftProfileIndex = new int[matrixProfile.length];
            rightProfileIndex = new int[matrixProfile.length];
        }

        for (int i = 0; i < mpSize; i++) {
            matrixProfile[i] = Double.POSITIVE_INFINITY;
            profileIndex[i] = -1;
            if (!isJoin) {
                rightMatrixProfile[i] = Double.POSITIVE_INFINITY;
                leftMatrixProfile[i] = Double.POSITIVE_INFINITY;
                rightProfileIndex[i] = -1;
                leftProfileIndex[i] = -1;
            }
        }

        sSize = Math.min(sSize, numQueries);
        int[] order = new int[numQueries];
        for (int i = 0; i < numQueries; i++) {
            order[i] = i;
        }
        Util.shuffleArray(order);
        int[] sampledOrder = new int[sSize];
        System.arraycopy(order, 0, sampledOrder, 0, sSize);
        var fft = Util.forwardFft(ts, false, 0, Util.padSize(dataSize));
        var mp = new BaseMatrixProfile(matrixProfile, profileIndex, rightMatrixProfile,
            leftMatrixProfile, rightProfileIndex, leftProfileIndex);
        for (var i : sampledOrder) {
            computeAnytime(i, isJoin, windowSize, exZone, ts, query, mp, distFunc, fft);
        }
        return mp;
    }

    public static <S extends WindowStatistic> void computeAnytime(int index, boolean isJoin,
        int windowSize, int exZone, RollingWindowStatistics<S> ts, RollingWindowStatistics<S> query,
        BaseMatrixProfile matrixProfile, DistanceProfileFunction<S> distFunc, Complex[] fft) {

        var dist = distFunc.apply(new DistanceProfileQuery<>(ts, query, index, windowSize, fft))
            .profile();
        var mpSize = matrixProfile.profile().length;
        var profile = matrixProfile.profile();
        var profileIndex = matrixProfile.indexes();
        var leftMatrixProfile = matrixProfile.leftProfile();
        var leftProfileIndex = matrixProfile.leftIndexes();
        var rightMatrixProfile = matrixProfile.rightProfile();
        var rightProfileIndex = matrixProfile.rightIndexes();

        if (exZone > 0) {
            int excSt = Math.max(0, index - exZone);
            int excEd = Math.min(mpSize - 1, index + exZone);
            for (int k = excSt; k <= excEd; k++) {
                dist[k] = Double.POSITIVE_INFINITY;
            }
        }

        for (var k = 0; k < mpSize; k++) {
            if (ts.stdDev(k) < Util.EPS || ts.skip(k) || ts.skip(index)) {
                dist[k] = Double.POSITIVE_INFINITY;
            }
        }

        if (query.stdDev(index) < Util.EPS) {
            dist[index] = Double.POSITIVE_INFINITY;
        }

        if (!isJoin) {
            // left matrixProfile
            for (int k = index; k < mpSize; k++) {
                if (dist[k] < leftMatrixProfile[k]) {
                    leftMatrixProfile[k] = dist[k];
                    leftProfileIndex[k] = index;
                }
            }

            // right matrixProfile
            for (int k = 0; k <= index; k++) {
                if (dist[k] < rightMatrixProfile[k]) {
                    rightMatrixProfile[k] = dist[k];
                    rightProfileIndex[k] = index;
                }
            }
        }

        // normal matrixProfile
        for (int k = 0; k < mpSize; k++) {
            if (dist[k] < profile[k]) {
                profile[k] = dist[k];
                profileIndex[k] = index;
            }
        }
    }

    public static <S extends WindowStatistic> MatrixProfile stamp(RollingWindowStatistics<S> ts,
        RollingWindowStatistics<S> query) {
        return stamp(ts, query, 0.5d, Integer.MAX_VALUE, new MASS2<>());
    }

    public static MatrixProfile of(double[] ts, double[] query, int windowSize) {
        var dataS = new BaseRollingWindowStatistics<>(windowSize, ts.length);
        var queryS = new BaseRollingWindowStatistics<>(windowSize, query.length);
        Arrays.stream(ts).forEach(dataS::apply);
        Arrays.stream(query).forEach(queryS::apply);
        return stamp(dataS, queryS);
    }

    public static MatrixProfile of(double[] ts, int windowSize) {
        var dataS = new BaseRollingWindowStatistics<>(windowSize, ts.length);
        Arrays.stream(ts).forEach(dataS::apply);
        return stamp(dataS, null);
    }
}
