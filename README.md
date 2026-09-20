# TSMP4J

**Matrix Profile algorithms for time series mining, in pure Java.**

[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)
![Java](https://img.shields.io/badge/Java-21-orange.svg)
![Version](https://img.shields.io/badge/version-1.1-informational.svg)

TSMP4J is a Java implementation of the [Matrix Profile](https://www.cs.ucr.edu/~eamonn/MatrixProfile.html)
family of algorithms and its downstream primitives. It is an attempt to port the algorithms collected
in [tsmp](https://github.com/matrix-profile-foundation/tsmp) (the reference R implementation of the
Matrix Profile Foundation) with a single design goal: every algorithm should consume the series as a
**stream**, not as an array you must materialise first.

The matrix profile — for every subsequence of a time series, the distance to its nearest non-overlapping
match plus that match's location — is a unifying primitive. Once you have it, motifs, discords, all-pairs
similarity joins, shapelets, semantic segmentation, anomaly scores and density estimates all become cheap
lookups on the same array.

## Contents

- [Highlights](#highlights)
- [Requirements](#requirements)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Reading the result](#reading-the-result)
- [Algorithms](#algorithms)
- [Design notes](#design-notes)
- [Building and testing](#building-and-testing)
- [Project layout](#project-layout)
- [References](#references)
- [License](#license)

## Highlights

- **Streaming first.** All matrix profile algorithms sit on top of a rolling-window statistics buffer, so
  data can be scored point by point as it arrives — no need to hold the whole series in memory.
- **Anytime by construction.** STAMP scores a random sample of queries, so a partial profile is available
  as soon as the first distance profiles are computed, and it sharpens as more queries are scored.
- **Incremental profiles.** STOMPI extends an existing profile as new points arrive instead of recomputing it.
- **Beyond the profile.** Contrast profiles, relative-frequency profiles, FLUSS segmentation, robust time
  series chains, contrast-profile classifiers, window-size selection and trend changepoint detection.
- **No native or GPU dependencies.** Pure Java 21 plus [Apache Commons Math](https://commons.apache.org/math/)
  for FFT and distributions.

## Requirements

| Requirement | Detail |
|---|---|
| Language level | Java 21 (a Gradle [toolchain](https://docs.gradle.org/current/userguide/toolchains.html) is configured in `build.gradle`) |
| Build tool | Gradle — the wrapper (`./gradlew`) is committed, pinned to Gradle 8.10 |
| Runtime dependencies | `org.apache.commons:commons-math3` (declared as `api`, transitively resolved) |

## Installation

### From your build tool

`build.gradle` publishes the library as `com.github.eugene-kamenev:tsmp4j:1.1` to Maven Local and to GitHub
Packages (`maven.pkg.github.com/eugene-kamenev/tsmp4j`). Resolve it from whichever of those you use.

Gradle:

```gradle
repositories {
    mavenCentral()
    maven {
        url = "https://maven.pkg.github.com/eugene-kamenev/tsmp4j"
        credentials {
            username = System.getenv("USERNAME")
            password = System.getenv("TOKEN")
        }
    }
}

dependencies {
    implementation "com.github.eugene-kamenev:tsmp4j:1.1"
}
```

Maven:

```xml
<repositories>
    <repository>
        <id>github</id>
        <url>https://maven.pkg.github.com/eugene-kamenev/tsmp4j</url>
    </repository>
</repositories>

<dependency>
    <groupId>com.github.eugene-kamenev</groupId>
    <artifactId>tsmp4j</artifactId>
    <version>1.1</version>
</dependency>
```

### From source

```bash
git clone https://github.com/eugene-kamenev/tsmp4j.git
cd tsmp4j
./gradlew publishToMavenLocal
```

This installs `com.github.eugene-kamenev:tsmp4j:1.1` (jar, sources jar and POM) into `~/.m2`, which you can
then resolve with `mavenLocal()`.

## Quick start

### Single batch

When the whole series is available up front, the static factories do everything:

```java
double[] data = ...;                 // your series
int windowSize = 10;

MatrixProfile mp = STAMP.of(data, windowSize);

System.out.println(mp.profile().length);   // data.length - windowSize + 1
```

### Streaming

Feed the algorithm point by point; the profile is computed from whatever has been buffered when you ask
for it:

```java
var windowSize = 10;
var bufferSize = 1024;                    // how much history the profile is computed over
var stamp = new STAMP(windowSize, bufferSize);

DoubleStream stream = ...;                // your data stream
stream.forEach(stamp::update);            // keeps rolling statistics current, does not compute MP yet

MatrixProfile mp = stamp.get();           // runs the MP algorithm over the buffered statistics
```

`update(...)` is available on every matrix profile algorithm; `stamp.get()` returns `null` until the
statistics buffer is full — see [Design notes](#design-notes).

### Similarity join (AB-join)

Pass a query series to find, for every query subsequence, its nearest neighbour in the reference series:

```java
MatrixProfile join = STAMP.of(reference, query, windowSize);
```

For AAMP the query is another statistics buffer:

```java
var aamp = new AAMP(new RollingWindowWithoutStatistics(windowSize, reference.length), 2.0);
Arrays.stream(reference).forEach(aamp::update);

var queryStats = new RollingWindowWithoutStatistics(windowSize, query.length);
Arrays.stream(query).forEach(queryStats::apply);

MatrixProfile join = aamp.get(queryStats);
```

### Incremental profile

```java
var initialStats = new BaseRollingWindowStatistics<BaseWindowStatistic>(windowSize, 200);
for (double v : firstBatch) {
    initialStats.apply(v);
}

var stompi = new STOMPI(initialStats, /* historySize */ 0);
stompi.update(nextPoint);                 // extend the profile instead of recomputing it
OnlineMatrixProfile mp = stompi.get();
```

## Reading the result

`MatrixProfile` exposes the profile itself plus the left/right variants used by chains, arcs and
segmentation:

| Accessor | Meaning |
|---|---|
| `profile()` / `indexes()` | matrix profile distance and the index of the matching subsequence |
| `leftProfile()` / `leftIndexes()` | nearest match strictly to the left |
| `rightProfile()` / `rightIndexes()` | nearest match strictly to the right |
| `windowSize()`, `exclusionZone()` | parameters the profile was computed with |

`OnlineMatrixProfile` additionally exposes `offset()`, and can be grown or trimmed with
`OnlineMatrixProfile.extend(...)` and `OnlineMatrixProfile.offset(...)`.

## Algorithms

### Matrix profile — z-normalised Euclidean distance

| Algorithm | Description |
|---|---|
| [STAMP](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/mp/stamp/STAMP.java) | Anytime matrix profile: random query sampling, usable profile after the first distance profile |
| [STOMP](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/mp/stomp/STOMP.java) | Ordered, FFT-based diagonal scan; optionally returns range indexes |
| [STOMPI](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/mp/stompi/STOMPI.java) | Real-time STOMP: extends an existing profile as new points arrive |
| [SKIMP](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/pmp/SKIMP.java) | Pan matrix profile — one profile per window size in a single pass |
| [MPX](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/mp/mpx/MPX.java) | Direct matrix profile without FFT |
| [MASS2](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/mp/mass/MASS2.java) | MP-DIST / MASS_V2: fast distance profile for a single query subsequence |

### Matrix profile — pure (non-normalised) Euclidean distance

| Algorithm | Description |
|---|---|
| [AAMP](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/mp/aamp/AAMP.java) | Approximate AAMP; the `p` parameter selects the distance (2.0 Euclidean, 1.0 Manhattan) |

### Contrast profiles and classification

| Algorithm | Description |
|---|---|
| [ContrastProfile](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/cp/ContrastProfileAlgorithm.java) | Contrast profile: closeness to repeated behaviour in a positive series versus distance from a negative one |
| [PanContrastProfile](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/cp/PanContrastProfileAlgorithm.java) | Contrast profile across several window sizes |
| [RelativeFrequencyMatrixProfile](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/cp/RelativeFrequencyMatrixProfileAlgorithm.java) | Relative-frequency matrix profile |
| [RelativeFrequencyContrastProfile](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/cp/RelativeFrequencyContrastProfileAlgorithm.java) | Relative-frequency contrast profile |
| [ContrastProfileClassifier](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/cp/ContrastProfileClassifier.java) | Nearest-plato classifier over each class' discriminative motif |
| [RelativeFrequencyContrastProfileClassifier](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/cp/RelativeFrequencyContrastProfileClassifier.java) | Same scheme over the relative-frequency contrast profile |

```java
int windowSize = 20;
double[][] trainingSeries = { class0Series, class1Series, class2Series };

var classifier = new ContrastProfileClassifier(windowSize, trainingSeries);

int predicted = classifier.classify(querySubsequence);
double[] distances = classifier.distances(querySubsequence);   // distance to every class' plato
double[] plato = classifier.plato(predicted);                  // the discriminative motif of that class
```

Each class' *plato* is extracted with its own series as the positive series and all other classes
concatenated as the negative series, so it is the subsequence that repeats inside its class yet stays far
from every other class. With three waveform classes (dip, peak, plateau, amplitude 5 on noise 0.1), a query
drawn from class 1 reports `[4.400, 0.166, 2.283]` and is predicted as class 1.

### Built on top of the profile

| Algorithm | Description |
|---|---|
| [FLUSS](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/fluss/FLUSS.java) | Domain-agnostic semantic segmentation: corrected arc counts and changepoints from a matrix profile |
| [TSC](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/tsc/TSC.java) | Robust time series chains — temporally ordered, continuously evolving motifs |

```java
MatrixProfile mp = STOMP.of(data, windowSize);
FLUSSCP segmentation = new FLUSS(windowSize, /* numSegments */ 3).apply(mp);

int[] changePoints = segmentation.changePoints();   // selected boundaries, strongest first
double[] arcCounts = segmentation.cac();            // corrected arc counts behind them
```

```java
var tsc = new TSC(windowSize, data.length, 0.25d);
Arrays.stream(data).forEach(tsc::update);
int[][] chains = tsc.get();
TSC.BestScore best = TSC.bestScore(data, chains, windowSize);
```

### Extras (not matrix-profile algorithms)

| Algorithm | Description |
|---|---|
| [MWF](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/extras/windowfinder/MWF.java) | Multi-Window-Finder: domain-agnostic subsequence length selection |
| [trendSegmentR](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/extras/tguw/Trend.java) | Detection of linear trend changes in a univariate series |
| [TGUW](src/main/java/com/github/eugene/kamenev/tsmp4j/algo/extras/tguw/TGUW.java) | Tail-greedy unbalanced Haar wavelet decomposition |

```java
int windowSize = MWF.mwf(series, /* lbound */ 10, /* ubound */ series.length / 4);

Trend trend = Trend.segment(series);          // trend.changePoints(), trend.estimates()
```

`ubound` must stay below the length of the series. For `Trend`, note the Javadoc convention: a reported
changepoint `p` corresponds to position `p - 1` in the input.

## Design notes

**Rolling statistics.** Every algorithm is built around
[`RollingWindowStatistics`](src/main/java/com/github/eugene/kamenev/tsmp4j/stats/RollingWindowStatistics.java),
which maintains the mean, standard deviation and value of the current window in a circular buffer using the
shifted-data variance algorithm. Algorithms are generic over the statistic type, so specialised statistics
(`MPXStatistic`, `NoStatistic` for non-normalised distances) plug into the same machinery.

**Buffer size bounds the profile.** `new STAMP(windowSize, bufferSize)` keeps at most `bufferSize`
subsequences; the returned profile has length `min(bufferSize, pointsFed) - windowSize + 1`, i.e. the profile
describes the most recent `bufferSize` points, not the whole stream. Choose `bufferSize` deliberately: it is
the memory/time budget of the computation.

**`get()` returns `null` until the buffer is full.** Check `isReady()` before consuming a profile, or feed at
least `bufferSize` points.

**Exclusion zone.** The exclusion zone is a fraction of the window size and defaults to `0.5`
(`floor(windowSize * 0.5)` subsequences on each side). Pass an explicit value to the constructor when a
different neighbourhood is required. In an AB-join the exclusion zone is disabled — reference and query are
independent series.

**Anytime control.** `STAMP.of(...)` and `get()` score every subsequence as a query, giving the exact
profile. To trade accuracy for latency, call the static overload
`STAMP.stamp(data, query, exclusionZone, exclusionZoneSize, sSize, distanceFunction)` with a bounded sample
size `sSize`: the profile then reflects the first `sSize` distance profiles and can be refined by scoring
more queries. Queries are visited in random order, so when several matches are equally close the reported
index can differ between runs.

**Degenerate windows.** Flat (zero-variance) windows and `NaN`/infinite input are marked as skipped and their
distance is reported as `+Infinity`, so they never become a nearest neighbour. Window sizes below 4 are
rejected.

## Building and testing

```bash
./gradlew build      # compile, run tests, assemble the jar
./gradlew test       # JUnit 5 test suite only
./gradlew publishToMavenLocal
```

Tests compare the computed profiles against committed reference profiles under `src/test/resources` (CSV).

## Project layout

```text
src/main/java/com/github/eugene/kamenev/tsmp4j/
├── algo/
│   ├── cp/        contrast profiles, relative-frequency profiles, classifiers
│   ├── extras/    MWF window finder, TGUW and trend segmentation
│   ├── fluss/     semantic segmentation
│   ├── mp/        matrix profile: stamp, stomp, stompi, mpx, mass, aamp
│   ├── pmp/       pan matrix profile (SKIMP)
│   └── tsc/       robust time series chains
├── stats/         rolling window statistics and window statistic types
└── utils/         FFT, buffers, shared helpers
```

## References

The implementations follow the published Matrix Profile papers and reference code; each class documents its
own citation.

- C.-C. M. Yeh *et al.* — *Matrix Profile I: All Pairs Similarity Joins for Time Series*, ICDM 2017 (STAMP).
- Y. Zhu *et al.* — *Matrix Profile II: Exploiting a Novel Algorithm and GPUs to Break the One Hundred Million
  Barrier for Time Series Motifs and Joins*, ICDM 2016 (STOMP).
- S. Gharghabi *et al.* — *Matrix Profile VIII: Domain Agnostic Online Semantic Segmentation at Superhuman
  Performance Levels*, ICDM 2017 (FLUSS).
- Y. Zhu, M. Imamura, D. Nikovski, E. Keogh — *Matrix Profile VII: Time Series Chains*, Knowledge and
  Information Systems 2018.
- *Matrix Profile XXIII: Contrast Profile: A Novel Time Series Primitive that Allows Real World
  Classification*, ICDM 2021, DOI [10.1109/ICDM51629.2021.00151](https://doi.org/10.1109/ICDM51629.2021.00151).
- A. Mueen — *MASS_V2: Mueen's Algorithm for Similarity Search*
  ([reference page](https://www.cs.unm.edu/~mueen/FastestSimilaritySearch.html)).
- *Pan Matrix Profile / SKIMP* ([project site](https://sites.google.com/view/pan-matrix-profile/home)).
- T. Mondal, R. Akbarinia, F. Masseglia — *Efficient Algorithms for Knowledge Discovery from Time Series*
  (AAMP), [journal version](https://mondal-tanmoy.github.io/files/pdf/journal/AAMP_Journal.pdf).
- S. Imani *et al.* — *Multi-Window-Finder: Domain Agnostic Window Size for Time Series Data*, MiLeTS 2021.
- H. Maeng, P. Fryzlewicz — *Detecting Linear Trend Changes in Data Sequences* (TGUW),
  [arXiv:1906.01939](https://arxiv.org/abs/1906.01939); port of [`trendsegmentR`](https://cran.r-project.org/web/packages/trendsegmentR/).
- *Robust Time Series Chain Discovery with Incremental Nearest Neighbors*,
  [arXiv:2211.02146](https://arxiv.org/abs/2211.02146).
- Reference R implementations: [tsmp](https://github.com/matrix-profile-foundation/tsmp), notably
  [`mpx.R`](https://github.com/matrix-profile-foundation/tsmp/blob/master/R/mpx.R).

## License

Apache License 2.0 — see [LICENSE](LICENSE).
