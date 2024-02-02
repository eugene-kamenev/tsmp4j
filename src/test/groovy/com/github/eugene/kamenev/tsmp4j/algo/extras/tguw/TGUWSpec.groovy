package com.github.eugene.kamenev.tsmp4j.algo.extras.tguw

import com.github.eugene.kamenev.tsmp4j.BaseSpec

class TGUWSpec extends BaseSpec {

    public static double ERROR = Math.pow(10, -9);

    private static double[] data = [
            10000.0, 10340.12, 10351.42, 10391.9, 10375.58, 10391.63, 10391.63, 10391.63,
            10296.04, 10297.89, 10316.3, 10304.32, 10256.09, 10149.47, 10147.84, 10248.89,
            10253.25, 10414.6, 10356.75, 10322.37, 10337.55, 10314.26, 10255.44, 10263.33,
            10255.58, 10308.87, 10233.85, 10289.48, 10358.69, 10307.0, 10338.61, 10328.53,
            10281.43, 10307.25, 10305.21, 10262.85, 10265.96, 10272.46, 10256.08, 10248.67,
            10240.32, 10254.56, 10254.63, 10241.79, 10208.34, 10179.1, 10225.4, 10213.69,
            10198.27, 10026.56, 10062.53, 10108.06, 10089.7, 10102.02, 10142.52, 10069.75,
            10073.2, 10093.88, 10056.42, 10028.63, 9995.45, 10067.46, 10005.18, 10016.05,
            10095.66, 10070.53, 10092.05, 10036.64, 10030.18, 9991.84, 10044.87, 10054.43,
            10068.2, 10077.6, 10102.87, 10161.51, 10146.21, 10159.55, 10155.22, 10124.55,
            10111.76, 10113.43, 10078.37, 10079.86, 10066.37, 10078.96, 10147.32, 10121.65,
            10153.96, 10170.93, 10239.21, 10343.84, 10305.07, 10326.58, 10269.94, 10295.8,
            10296.24, 10308.8, 10320.73, 10388.87, 10385.45, 10415.13, 10387.34, 10335.68,
            10365.43, 10364.64, 10336.83, 10291.89, 10295.54, 10284.03, 10278.48, 10306.43,
            10272.91, 10306.52, 10323.98, 10307.01, 10286.2, 10296.57, 10296.3, 10226.61,
            10220.93, 10210.25, 10253.79, 10287.11, 10277.15, 10341.34, 10309.03, 10298.45
    ]

    public static final double[][] merging_hist_1 = load2D('merging_hist_1.csv', TGUWSpec)

    public static final double[][] merging_hist_2 = load2D('merging_hist_2.csv', TGUWSpec)

    public static final double[][] merging_hist_3 = load2D('merging_hist_3.csv', TGUWSpec)

    public static final double[][] merging_hist_4 = load2D('merging_hist_4.csv', TGUWSpec)

    public static final double[] tsCoeffs = load1D('ts_coeffs.csv', TGUWSpec)

    public static final double[] twotogether = load1D('twotogether.csv', TGUWSpec)

    public static final double[] est1 = load1D('est1.csv', TGUWSpec)
    public static final double[] est2 = load1D('est2.csv', TGUWSpec)
    public static final double[] est3 = load1D('est3.csv', TGUWSpec)
    public static final double[] est4 = load1D('est4.csv', TGUWSpec)
    public static final double[] est5 = load1D('est5.csv', TGUWSpec)
    public static final double[] est6 = load1D('est6.csv', TGUWSpec)



    def 'test tail greedy unbalance wavelet decomposition'() {
        given:
        var mergingHist = new double[][][]{merging_hist_1, merging_hist_2, merging_hist_3, merging_hist_4};
        var tguw = new TGUW(data, 0.04);

        when:
        tguw.forward()
        var decomp2 = tguw.getDecomposition();

        then:
        equals(decomp2.twotogether().stream().mapToDouble(e -> e.doubleValue()).toArray(), twotogether)
        equals(decomp2.transformed(), tsCoeffs)

        for (int i = 0; i < mergingHist.length; i++) {
            var actual = mergingHist[i];
            double[][] toCheck = switch (i) {
                case 0 -> Arrays.stream(decomp2.mergeHistory().indexes()).map(x -> Arrays.stream(x).mapToDouble(y -> (double) y).toArray()).toArray(double[][]::new);
                case 1 -> decomp2.mergeHistory().details();
                case 2 -> decomp2.mergeHistory().smooths();
                case 3 -> decomp2.mergeHistory().balanced();
                default -> throw new IllegalStateException("Unexpected value: " + i);
            };
            for (int j = 0; j < actual.length; j++) {
                equals(actual[j], toCheck[j])
                println "Merging hist for ${i} and ${j} is OK"
            }
        }
    }

    def 'test trend segmentr estimates and change points'() {
        given:
        var changePoints1 = [1, 17, 49, 75, 86, 101]
        var changePoints2 = [49, 86, 101]

        when:
        var r1 = Trend.segment(data);
        var r2 = Trend.segment(data, true, true, true, 0.04, -1, 0, -1, false);
        var r3 = Trend.segment(data, false, false, true, 0.04, -1, 0, -1, false);
        var r4 = Trend.segment(data, false, true, true, 0.04, -1, 0, -1, false);
        var r5 = Trend.segment(data, true, false, false, 0.04, -1, 0, -1, false);
        var r6 = Trend.segment(data, true, false, true, 0.04, -1, 0, -1, false);

        then:
        r1.changePoints() == changePoints2
        r2.changePoints() ==  changePoints1
        r3.changePoints() == changePoints2
        r4.changePoints() == changePoints2
        r5.changePoints() == changePoints1
        r6.changePoints() == changePoints1
        equals(r1.estimates(), est1)
        equals(r2.estimates(), est2)
        equals(r3.estimates(), est3)
        equals(r4.estimates(), est4)
        equals(r5.estimates(), est5)
        equals(r6.estimates(), est6)
    }

    @Override
    boolean equals(double[] a, double[] b) {
        return super.equals(a, b, ERROR)
    }

}
