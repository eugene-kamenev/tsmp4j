## Time Series with Matrix Profile in Java

This repository contains Matrix Profile algorithms implemented in Java.
It is an attempt to port algorithms presented in [tsmp](https://github.com/matrix-profile-foundation/tsmp)
and [stumpy](https://github.com/matrix-profile-foundation/tsmp).

[The Matrix Profile](https://www.cs.ucr.edu/~eamonn/MatrixProfile.html), has the potential to revolutionize time series data mining because of its generality, 
versatility, simplicity and scalability. 
In particular, it has implications for time series motif discovery, time series joins, shapelet discovery (classification), density estimation, semantic segmentation, visualization, rule discovery, clustering etc.

This library includes the following algorithms:
1. [STAMP](/src/main/java/com/github/eugene/kamenev/tsmp4j/STAMP.java)
2. [SKIMP](/src/main/java/com/github/eugene/kamenev/tsmp4j/SKIMP.java)
3. [MPX](/src/main/java/com/github/eugene/kamenev/tsmp4j/MPX.java)
4. [MP-DIST (MASS2)](/src/main/java/com/github/eugene/kamenev/tsmp4j/MASS2.java)

More algorithms will be added in the future.
