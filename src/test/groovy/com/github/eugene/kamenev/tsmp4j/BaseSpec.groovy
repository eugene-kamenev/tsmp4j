package com.github.eugene.kamenev.tsmp4j


import spock.lang.Specification

import java.util.function.Function

import static org.hamcrest.MatcherAssert.assertThat
import static org.hamcrest.Matchers.closeTo
import static org.hamcrest.Matchers.equalTo

class BaseSpec extends Specification {

    /**
     * Note that due to using double type for calculations,
     * we experience some error in calculations
     */
    public static double ERROR = Math.pow(10, -11);

    public static double[][] load2D(String file, Class clazz, int skip = 1) {
        clazz.getResource(file).readLines().stream()
            .skip(skip)
            .filter(s -> !s.trim().isEmpty())
            .map(it -> Arrays.stream(it.split(','))
            .mapToDouble(s -> Double.parseDouble(s)).toArray())
            .toArray(double[][]::new);
    }

    public static double[] load1D(String file, Class clazz, int skip = 1) {
        clazz.getResource(file).readLines().stream()
            .skip(skip)
            .filter(s -> !s.trim().isEmpty())
            .mapToDouble(it -> Double.parseDouble(it))
            .toArray();
    }

    public static final List<ToyData> data = loadData('mp_toy_data.csv', (rows) -> {
        return Arrays.stream(rows).map(s -> {
            return new ToyData(s[0] as double, s[1] as double, s[2] as double)
        }).toList()
    }, BaseSpec)

    boolean equals(double[] a, double[] b, double th = ERROR) {
        if (a.length != b.length) {
            throw new IllegalStateException("Arrays have different length");
        }
        for (int i = 0; i < a.length; i++) {
            if (!a[i].equals(b[i])) {
                assertThat("On index ${i}", a[i], closeTo(b[i], th))
            }
        }
        return true
    }

    boolean equals(int[] a, int[] b) {
        if (a.length != b.length) {
            throw new IllegalStateException("Arrays have different length");
        }
        for (int i = 0; i < a.length; i++) {
            assertThat("On index ${i}", a[i], equalTo(b[i]))
        }
        return true
    }

    static <T> T loadData(String file, Function<String[][], T> transform, Class clazz = this.getClass()) {
        String[][] data = clazz.getResourceAsStream(file)
                .readLines()
                .stream()
                .skip(1)
                .map(s -> s.split(','))
                .toArray(String[][]::new)
        return transform.apply(data)
    }

    static MP loadMP(String file, Class clazz = this.getClass()) {
        loadData(file, (rows) -> {
            double[] mp = new double[rows.length]
            int[] pi = new int[rows.length]
            int counter = 0
            Arrays.stream(rows).forEach(row -> {
                mp[counter] = row[0] as double
                pi[counter++] = (row[1] as int) - 1 // in R and Matlab indexing starts from 1
            })
            return new MP(mp, pi)
        }, clazz)
    }

    static record ToyData(double x, double y, double z) {}

    static record MP(double[] mp, int[] pi) {}

    static record MPDist(double[] x) {
        static MPDist load(String path, Class clazz) {
            return loadData(path, (rows) -> {
                double[] x = new double[rows.length]
                int counter = 0
                Arrays.stream(rows).forEach(row -> {
                    x[counter++] = row[0] as double
                })
                return new MPDist(x)
            }, clazz)
        }
    }

    static double parseDouble(String value) {
        if (value == 'Inf') {
            return Double.POSITIVE_INFINITY
        } else if (value == '-Inf') {
            return Double.NEGATIVE_INFINITY
        } else {
            return value as double
        }
    }

    static int parseInt(String value) {
        if (value == '-Inf') {
            return -1
        }
        return value as int
    }
}
