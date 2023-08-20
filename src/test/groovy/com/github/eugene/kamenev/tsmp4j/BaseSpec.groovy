package com.github.eugene.kamenev.tsmp4j

import spock.lang.Shared
import spock.lang.Specification

import java.util.function.Function

import static org.hamcrest.MatcherAssert.assertThat
import static org.hamcrest.Matchers.closeTo
import static org.hamcrest.Matchers.equalTo

class BaseSpec extends Specification {

    public static double ERROR = 0.0000001d

    @Shared
    protected List<ToyData> data

    def setupSpec() {
        this.data = loadData('mp_toy_data.csv', (rows) -> {
            return Arrays.stream(rows)
                    .map(s -> {
                        return new ToyData(s[0] as double, s[1] as double, s[2] as double)
                    }).toList()
        }, BaseSpec)
    }


    double[] getTimeSeries() {
        return new double[]{
                36.3, 55.57, 56.91, 57.48, 61.13, 74.48, 77.75, 63.1, 67.98, 69.13, 49.98, 47.72, 35.07, 41.95,
                60.76, 55.02, 42.88, 46.02, 42.8, 38.26, 37.86, 23.48, 23.45, 37.65, 27.14, 31.55, 48.78, 42.91,
                38.04, 45.37, 54.68, 53.94, 44.3, 48.53, 63.52, 48.78, 35.56, 41.58, 27.52, 18.24, 4.16, 4.09,
                23.2, 22.75, 23.21, 24.34, 25.39, 21.66, 20.67, 33.89, 33.16, 29.5, 38.93, 75.18, 79.44, 74.13,
                78.61, 88.71, 89.83, 79.87, 81.25, 81.61, 68.68, 60.17, 65.88, 63.51, 52.7, 50.86, 48.77, 43.64,
                38.32, 49.64, 51.77, 68.38, 58.36, 51.4, 48.75, 33.83, 30.24, 47.32, 54.83, 47.32, 35.34, 46.15,
                54.79, 52.17, 51.32, 41.17, 37.89, 28.14, 27.46, 15.66, 36.84, 26.96, 23.88, 21.77, 37.32, 25.93,
                31.87, 34.94, 49.85, 43.42, 33.83, 37.45, 46.01, 40.1, 34.77, 38.01, 45.47, 35.45, 22.09, 23.45,
                19.28, 17.96, 45.1, 47.15, 68.33, 74.77, 62.24, 61.74, 61.17, 67.89, 50.45, 44.3, 46.25, 56.41,
                55.84, 42.2
        }
    }

    double[] getQuery() {
        return new double[]{38.04, 45.37, 54.68, 53.94, 44.3, 48.53}
    }

    boolean equals(double[] a, double[] b, double th = ERROR) {
        if (a.length != b.length) {
            return false
        }
        for (int i = 0; i < a.length; i++) {
            if (a[i] != b[i]) {
                assertThat("On index ${i}", a[i], closeTo(b[i], th))
            }
        }
        return true
    }

    boolean equals(int[] a, int[] b) {
        if (a.length != b.length) {
            return false
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
            double[] mp = new double[rows.size()]
            int[] pi = new int[rows.size()]
            int counter = 0
            Arrays.stream(rows).forEach(row -> {
                mp[counter] = row[0] as double
                pi[counter++] = (row[1] as int) - 1 // in R indexing starts from 1
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
