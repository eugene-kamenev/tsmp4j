package com.github.eugene.kamenev.tsmp4j;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;
import java.util.List;
import java.util.function.Function;

import static org.junit.jupiter.api.Assertions.assertEquals;

public abstract class BaseTest {

    /**
     * Note that due to using double type for calculations,
     * we experience some error in calculations
     */
    public static final double ERROR = Math.pow(10, -11);

    public static final List<ToyData> data = loadData("mp_toy_data.csv", rows -> Arrays.stream(rows)
        .map(s -> new ToyData(Double.parseDouble(s[0]), Double.parseDouble(s[1]), Double.parseDouble(s[2])))
        .toList(), BaseTest.class);

    public static double[][] load2D(String file, Class<?> clazz) {
        return load2D(file, clazz, 1);
    }

    public static double[][] load2D(String file, Class<?> clazz, int skip) {
        return readLines(file, clazz).stream()
            .skip(skip)
            .filter(s -> !s.trim().isEmpty())
            .map(it -> Arrays.stream(it.split(","))
                .mapToDouble(Double::parseDouble)
                .toArray())
            .toArray(double[][]::new);
    }

    public static double[] load1D(String file, Class<?> clazz) {
        return load1D(file, clazz, 1);
    }

    public static double[] load1D(String file, Class<?> clazz, int skip) {
        return readLines(file, clazz).stream()
            .skip(skip)
            .filter(s -> !s.trim().isEmpty())
            .mapToDouble(Double::parseDouble)
            .toArray();
    }

    public static <T> T loadData(String file, Function<String[][], T> transform, Class<?> clazz) {
        String[][] rows = readLines(file, clazz).stream()
            .skip(1)
            .map(s -> s.split(","))
            .toArray(String[][]::new);
        return transform.apply(rows);
    }

    public static MP loadMP(String file, Class<?> clazz) {
        return loadData(file, rows -> {
            double[] mp = new double[rows.length];
            int[] pi = new int[rows.length];
            for (int i = 0; i < rows.length; i++) {
                mp[i] = Double.parseDouble(rows[i][0]);
                pi[i] = Integer.parseInt(rows[i][1]) - 1; // in R and Matlab indexing starts from 1
            }
            return new MP(mp, pi);
        }, clazz);
    }

    public static double parseDouble(String value) {
        if ("Inf".equals(value)) {
            return Double.POSITIVE_INFINITY;
        } else if ("-Inf".equals(value)) {
            return Double.NEGATIVE_INFINITY;
        }
        return Double.parseDouble(value);
    }

    public static int parseInt(String value) {
        if ("-Inf".equals(value)) {
            return -1;
        }
        return Integer.parseInt(value);
    }

    /**
     * Tolerance used by {@link #equals(double[], double[])}. Subclasses may narrow or widen it.
     *
     * @return allowed deviation between expected and actual values
     */
    protected double error() {
        return ERROR;
    }

    protected void equals(double[] a, double[] b) {
        equals(a, b, error());
    }

    protected void equals(double[] a, double[] b, double th) {
        if (a.length != b.length) {
            throw new IllegalStateException("Arrays have different length");
        }
        for (int i = 0; i < a.length; i++) {
            if (!Double.valueOf(a[i]).equals(Double.valueOf(b[i]))) {
                assertEquals(b[i], a[i], th, "On index " + i);
            }
        }
    }

    protected void equals(int[] a, int[] b) {
        if (a.length != b.length) {
            throw new IllegalStateException("Arrays have different length");
        }
        for (int i = 0; i < a.length; i++) {
            assertEquals(b[i], a[i], "On index " + i);
        }
    }

    private static List<String> readLines(String file, Class<?> clazz) {
        var resource = clazz.getResource(file);
        if (resource == null) {
            throw new IllegalArgumentException("Test resource not found: " + file + " for " + clazz.getName());
        }
        try (var reader = new BufferedReader(new InputStreamReader(resource.openStream(), StandardCharsets.UTF_8))) {
            return reader.lines().toList();
        } catch (IOException e) {
            throw new UncheckedIOException(e);
        }
    }

    public record ToyData(double x, double y, double z) { }

    public record MP(double[] mp, int[] pi) { }

    public record MPDist(double[] x) {

        public static MPDist load(String path, Class<?> clazz) {
            return loadData(path, rows -> {
                double[] x = new double[rows.length];
                for (int i = 0; i < rows.length; i++) {
                    x[i] = Double.parseDouble(rows[i][0]);
                }
                return new MPDist(x);
            }, clazz);
        }
    }
}
