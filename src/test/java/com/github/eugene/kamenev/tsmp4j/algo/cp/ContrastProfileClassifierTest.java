package com.github.eugene.kamenev.tsmp4j.algo.cp;

import java.util.Random;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

class ContrastProfileClassifierTest extends com.github.eugene.kamenev.tsmp4j.BaseTest {

    private static final int WINDOW = 20;

    private static final double AMP = 5.0d;

    private static final double NOISE = 0.1d;

    // Each class is defined by a distinct waveform: a dip, a peak and a rectangular plateau.
    private static final double[][] SHAPES = {
        shape("dip"), shape("peak"), shape("plateau")
    };

    @Test
    void contrastProfileClassifiesEachWaveform() {
        var training = trainSeries();
        var classifier = new ContrastProfileClassifier(WINDOW, training);

        assertEquals(3, classifier.numClasses());
        assertEquals(WINDOW, classifier.windowSize());

        int correct = 0;
        int total = 0;
        for (int cls = 0; cls < SHAPES.length; cls++) {
            for (double[] test : testSeries(cls)) {
                total++;
                var distances = classifier.distances(test);
                int predicted = classifier.classify(test);
                correct += predicted == cls ? 1 : 0;

                assertEquals(cls, predicted,
                    "Waveform " + cls + " misclassified, distances=" + describe(distances));
                assertTrue(distances[cls] < 3.0d,
                    "Match to own class should be tight, got " + distances[cls]);
                for (int other = 0; other < distances.length; other++) {
                    if (other != cls) {
                        assertTrue(distances[cls] < distances[other],
                            "Own class distance must beat class " + other);
                    }
                }
            }
        }
        assertEquals(total, correct);
    }

    @Test
    void relativeFrequencyContrastProfileClassifiesEachWaveform() {
        var training = trainSeries();
        var classifier =
            new RelativeFrequencyContrastProfileClassifier(WINDOW, 3, training);

        assertEquals(3, classifier.numClasses());
        assertEquals(3, classifier.maxFreq());

        int correct = 0;
        int total = 0;
        for (int cls = 0; cls < SHAPES.length; cls++) {
            for (double[] test : testSeries(cls)) {
                total++;
                var distances = classifier.distances(test);
                int predicted = classifier.classify(test);
                correct += predicted == cls ? 1 : 0;

                assertEquals(cls, predicted,
                    "Waveform " + cls + " misclassified, distances=" + describe(distances));
                assertTrue(distances[cls] < 3.0d,
                    "Match to own class should be tight, got " + distances[cls]);
                for (int other = 0; other < distances.length; other++) {
                    if (other != cls) {
                        assertTrue(distances[cls] < distances[other],
                            "Own class distance must beat class " + other);
                    }
                }
            }
        }
        assertEquals(total, correct);
    }

    @Test
    void platoRepresentsAndDiscriminatesItsClass() {
        var training = trainSeries();
        var classifier = new ContrastProfileClassifier(WINDOW, training);

        for (int c = 0; c < SHAPES.length; c++) {
            var plato = classifier.plato(c);
            assertEquals(WINDOW, plato.length);

            // The plato extracted for class c must be reproducible inside class c's own data:
            // class c's training series should be closest to plato c.
            var distancesToOwnData = classifier.distances(training[c]);
            assertEquals(c, argMin(distancesToOwnData),
                "Plato " + c + " is not the best match for its own class: "
                    + describe(distancesToOwnData));
            assertTrue(distancesToOwnData[c] < 3.0d,
                "Plato " + c + " does not match its own class closely");

            // ...and must clearly beat the platos of every other class on that same data.
            for (int other = 0; other < distancesToOwnData.length; other++) {
                if (other != c) {
                    assertTrue(distancesToOwnData[c] < distancesToOwnData[other],
                        "Plato " + c + " is not discriminative against class " + other);
                }
            }
        }
    }

    private static int argMin(double[] values) {
        int best = 0;
        for (int i = 1; i < values.length; i++) {
            if (values[i] < values[best]) {
                best = i;
            }
        }
        return best;
    }

    private static double[][] trainSeries() {
        var series = new double[SHAPES.length][];
        for (int c = 0; c < SHAPES.length; c++) {
            series[c] = build(200, 100 + c, new int[]{40, 120}, SHAPES[c]);
        }
        return series;
    }

    private static double[][] testSeries(int cls) {
        var seeds = new int[]{1001, 1002, 1003};
        var positions = new int[]{30, 60, 90};
        var series = new double[seeds.length][];
        for (int i = 0; i < seeds.length; i++) {
            series[i] = build(140, seeds[i] + cls, new int[]{positions[i]}, SHAPES[cls]);
        }
        return series;
    }

    private static double[] build(int length, int seed, int[] embedPositions, double[] shape) {
        var random = new Random(seed);
        var ts = new double[length];
        for (int i = 0; i < length; i++) {
            ts[i] = (random.nextDouble() * 2 - 1) * NOISE;
        }
        for (int pos : embedPositions) {
            for (int k = 0; k < WINDOW; k++) {
                ts[pos + k] += AMP * shape[k];
            }
        }
        return ts;
    }

    private static double[] shape(String kind) {
        var s = new double[WINDOW];
        double center = (WINDOW - 1) / 2.0d;
        for (int i = 0; i < WINDOW; i++) {
            s[i] = switch (kind) {
                case "dip" -> -gaussian(i, center, 3.0d);
                case "peak" -> gaussian(i, center, 3.0d);
                default -> (i >= 5 && i <= 14) ? 1.0d : 0.0d;
            };
        }
        return s;
    }

    private static double gaussian(double x, double mu, double sigma) {
        return Math.exp(-((x - mu) * (x - mu)) / (2 * sigma * sigma));
    }

    private static String describe(double[] distances) {
        var sb = new StringBuilder("[");
        for (int i = 0; i < distances.length; i++) {
            sb.append(i == 0 ? "" : ", ").append(String.format("%.3f", distances[i]));
        }
        return sb.append(']').toString();
    }
}
