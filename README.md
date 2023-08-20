## Time Series with Matrix Profile in Java

This repository contains Matrix Profile algorithms implemented in Java.
It is an attempt to port algorithms presented
in [tsmp](https://github.com/matrix-profile-foundation/tsmp).

[The Matrix Profile](https://www.cs.ucr.edu/~eamonn/MatrixProfile.html), has the potential to
revolutionize time series data mining because of its generality,
versatility, simplicity and scalability.
In particular, it has implications for time series motif discovery, time series joins, shapelet
discovery (classification), density estimation, semantic segmentation, visualization, rule
discovery, clustering etc.

This library includes the following algorithms:

1. [STAMP](/src/main/java/com/github/eugene/kamenev/tsmp4j/algo/mp/stamp/STAMP.java)
2. [STOMP](/src/main/java/com/github/eugene/kamenev/tsmp4j/algo/mp/stomp/STOMP.java)
3. [SKIMP](/src/main/java/com/github/eugene/kamenev/tsmp4j/algo/pmp/SKIMP.java)
4. [MPX](/src/main/java/com/github/eugene/kamenev/tsmp4j/algo/mp/mpx/MPX.java)
5. [MP-DIST (MASS2)](/src/main/java/com/github/eugene/kamenev/tsmp4j/algo/mp/mass/MASS2.java)

More algorithms will be added in the future.

## Usage
All algorithms here are built around [RollingWindowStatistics](/src/main/java/com/github/eugene/kamenev/tsmp4j/stats/RollingWindowStatistics.java) object. 
It simply computes statistics required to run Matrix Profile algorithms on the fly, in a circular buffer manner, hence allows streaming data processing out of the box.

### Streaming case
```java
DoubleStream stream = ... // your data stream
var bs = 1024; // statistics buffer size
var w = 10; // window for MP algorithm
var stamp = new STAMP(w, bs);

stream.forEach(stamp::update); // it keeps statistics updated, not the matrix profile

var matrixProfile = stamp.get(); // execute MP algorithm for statistics collected
```

### Single batch case
```java
double[] data = ... // your data
var w = 10; // window for MP algorithm    
var matrixProfile = STAMP.of(data, w);

```
