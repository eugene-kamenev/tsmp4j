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

    private static double[] data2 = [
            36.995930447650764, 35.06311360448808, 51.92107995846314, 43.78283712784589,
            53.879310344827594, 69.4815606627472, 48.85197850512946, 29.717682020802386,
            58.247286205983585, 78.03790412486065, 56.73758865248227, 61.13870844478411,
            42.643923240938165, 44.84304932735426, 34.22313483915127, 50.98572399728077,
            62.47559547051933, 66.97578567748583, 83.33333333333333, 100.0, 100.0, 100.0,
            25.000000000000014, 28.57142857142857, 69.1358024691358, 91.34948096885813,
            96.15384615384616, 97.54562617998741, 100.0, 62.50000000000001, 66.66666666666666,
            60.0, 47.05882352941176, 54.229934924078094, 54.229934924078094, 45.45454545454546,
            37.73584905660377, 59.25925925925926, 64.51612903225806, 85.33333333333333,
            45.80852038479157, 75.25083612040135, 100.0, 22.22222222222223,
            20.0, 27.270248159258244, 51.20702267739576, 23.679848448969935, 42.678047479327816,
            52.770448548812666, 100.0, 80.0, 37.495313085864275, 65.73541495480691, 81.19079837618403,
            82.78145695364239, 67.93478260869566, 85.57457212713936, 100.0, 25.000000000000014,
            48.27586206896552, 37.495313085864275, 70.57163020465774, 92.83135636926251,
            28.01512816921138, 29.45508100147275, 23.150827642088203, 14.328700386874914, 18.234865061998548,
            30.845157310302284, 18.10938065918147, 45.33091568449683, 57.012542759407076, 47.01457451810061,
            69.7612020391736, 78.003120124805, 100.0, 57.142857142857146, 77.61194029850746, 59.25925925925926,
            44.44444444444444, 54.55537370430987, 29.629629629629633, 50.26704366949419, 67.58832565284177,
            65.8544616397761, 80.94186902133922, 89.46726311508743, 100.0, 100.0, 100.0, 100.0, 100.0, 40.0,
            100.0, 100.0, 100.0, 100.0, 75.0, 44.44444444444444, 33.33333333333333, 33.33333333333333, 40.0,
            68.9655172413793, 47.755491881566385, 48.309178743961354, 43.440486533449175, 64.82982171799028, 76.64380798709158,
            72.15007215007215, 83.82229673093042, 100.0, 50.0, 50.0, 70.58823529411765, 86.48648648648648, 88.88888888888889,
            94.16195856873823, 100.0, 100.0, 100.0, 66.66666666666666, 100.0, 100.0, 100.0, 100.0, 100.0,
            8.695652173913047
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

    def 'test trend segmentr no changepoints result'() {
        when:
        var r = Trend.segment(data2);

        then:
        r.changePoints().isEmpty()
    }

    @Override
    boolean equals(double[] a, double[] b) {
        return super.equals(a, b, ERROR)
    }

}
