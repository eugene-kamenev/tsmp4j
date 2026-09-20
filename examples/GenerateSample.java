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

import com.github.eugene.kamenev.tsmp4j.algo.mp.MatrixProfile;
import com.github.eugene.kamenev.tsmp4j.algo.mp.stomp.STOMP;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Random;

/**
 * Generates a reproducible sample time series and its matrix profile with tsmp4j, then writes both as CSV
 * so {@code plot.py} can render the README figure.
 *
 * <p>The series is built from two motifs that repeat several times — a sine wave (motif A) and a sawtooth
 * (motif B) — with a single, unique "discord" waveform spliced into one repetition. Because motif instances
 * recur, their matrix-profile distances are near zero, while the discord, which matches nothing, shows up as
 * the profile's maximum. A fixed seed keeps the output byte-for-byte reproducible.
 *
 * <p>Run (see the "Generating the sample figure" section of the README):
 * <pre>{@code
 * ./gradlew jar
 * CMATH=$(find "$HOME/.gradle" -name commons-math3-3.6.1.jar | head -1)
 * CP="build/libs/tsmp4j-1.2.jar:$CMATH"
 * javac -cp "$CP" -d /tmp/examples-classes examples/GenerateSample.java
 * java -cp "$CP:/tmp/examples-classes" GenerateSample examples/data
 * python3 examples/plot.py
 * }</pre>
 */
public final class GenerateSample {

    private static final int WINDOW_SIZE = 30;
    private static final int MOTIF_LENGTH = 60;
    private static final int REPEATS = 8;
    private static final int DISCORD_REPETITION = 4;
    private static final double NOISE_STD = 0.05;
    private static final long SEED = 42L;

    private GenerateSample() {
    }

    public static void main(String[] args) throws IOException {
        Path outDir = Path.of(args.length > 0 ? args[0] : "examples/data");
        Files.createDirectories(outDir);

        int period = 2 * MOTIF_LENGTH;
        int length = REPEATS * period;
        double[] series = new double[length];

        for (int r = 0; r < REPEATS; r++) {
            int startA = r * period;
            fillSine(series, startA, MOTIF_LENGTH, 1.0);
            fillSawtooth(series, startA + MOTIF_LENGTH, MOTIF_LENGTH, 1.5);
        }
        fillDiscord(series, DISCORD_REPETITION * period, MOTIF_LENGTH);

        Random random = new Random(SEED);
        for (int i = 0; i < length; i++) {
            series[i] += random.nextGaussian() * NOISE_STD;
        }

        MatrixProfile mp = STOMP.of(series, WINDOW_SIZE);
        double[] profile = mp.profile();
        int[] indexes = mp.indexes();

        writeSeries(outDir.resolve("timeseries.csv"), series);
        writeProfile(outDir.resolve("matrix_profile.csv"), profile, indexes);

        int discord = argMax(profile);
        System.out.printf(
            "series length=%d, window=%d, profile length=%d%n",
            length, mp.windowSize(), profile.length);
        System.out.printf(
            "discord at subsequence %d (distance %.3f, matched to %d)%n",
            discord, profile[discord], indexes[discord]);
    }

    private static void fillSine(double[] series, int start, int len, double amp) {
        for (int k = 0; k < len; k++) {
            series[start + k] += amp * Math.sin(2.0 * Math.PI * k / len);
        }
    }

    private static void fillSawtooth(double[] series, int start, int len, double amp) {
        for (int k = 0; k < len; k++) {
            series[start + k] += amp * (2.0 * k / len - 1.0);
        }
    }

    private static void fillDiscord(double[] series, int start, int len) {
        for (int k = 0; k < len; k++) {
            series[start + k] += Math.exp(-k / (len / 3.0)) * Math.sin(2.0 * Math.PI * k / 8.0);
        }
    }

    private static int argMax(double[] values) {
        int best = 0;
        for (int i = 1; i < values.length; i++) {
            if (values[i] > values[best]) {
                best = i;
            }
        }
        return best;
    }

    private static void writeSeries(Path file, double[] series) throws IOException {
        StringBuilder sb = new StringBuilder("index,value\n");
        for (int i = 0; i < series.length; i++) {
            sb.append(i).append(',').append(series[i]).append('\n');
        }
        Files.writeString(file, sb, StandardCharsets.UTF_8);
    }

    private static void writeProfile(Path file, double[] profile, int[] indexes) throws IOException {
        StringBuilder sb = new StringBuilder("index,profile,match_index\n");
        for (int i = 0; i < profile.length; i++) {
            sb.append(i).append(',').append(profile[i]).append(',').append(indexes[i]).append('\n');
        }
        Files.writeString(file, sb, StandardCharsets.UTF_8);
    }
}
