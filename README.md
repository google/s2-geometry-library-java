# S2 Geometry Library

## Overview

This is a package for manipulating geometric shapes. Unlike many geometry
libraries, S2 is primarily designed to work with _spherical geometry_, i.e.,
shapes drawn on a sphere rather than on a planar 2D map. This makes it
especially suitable for working with geographic data.

If you want to learn more about the library, start by reading the
[overview](http://s2geometry.io/about/overview) and [quick start
document](http://s2geometry.io/devguide/cpp/quickstart), then read the
introduction to the [basic types](http://s2geometry.io/devguide/basic_types).

S2 documentation can be found on [s2geometry.io](http://s2geometry.io).

## Build and Install

You may either download the source as a ZIP archive, or [clone the git
repository](https://help.github.com/articles/cloning-a-repository/).
The Java packages are built and tested using [Maven](https://maven.apache.org/).

In the directory containing the ```pom.xml``` file, use Maven to
compile the package, run tests, and install the package. For example:

```
mvn clean
mvn compile
mvn test
mvn package
mvn install
```

### Benchmarks

After building packages, a "benchmarks.jar" file may be found in
```benchmarks/target/```. Benchmarks can then be run with the command
```java -jar benchmarks/target/benchmarks.jar```
Parameters for the benchmarks can be passed on the command line. For example, to
run just the S2Loop benchmarks,
```java -jar benchmarks/target/benchmarks.jar S2Loop```

## S2 implementations

The S2 library has implementations in several languages. In addition to this
Java implementation, Google provides:

* [C++](https://github.com/google/s2geometry) (The reference implementation
  and the most full featured).
* [Go](https://github.com/golang/geo) (Approximately 40% complete.)
* [Python](https://github.com/google/s2geometry/tree/master/src/python)

We (the S2 developers) aim to provide similar classes and APIs across all
implementations, while adapting to language idioms where that makes sense, and
changing APIs where required for good performance. The implementations have
varying degrees of completeness and maturity. This Java implementation is
heavily used within Google and is generally mature, aside from the newest
features, but is not as complete as C++.

## 2022 Q4 Release Highlights

Many improvements have been made to the Java implementation of S2 since the last
release in September 2021. Some highlights:

*   Implementations of S2ClosestEdgeQuery and S2FurthestEdgeQuery, as well as
    S2ChainInterpolationQuery and S2HausdorffDistanceQuery.

*   New tools for clustering and sharding data: S2RegionSharder and
    S2DensityTree.

*   Support for map projections in MercatorProjection, PlateCarreeProjection,
    and S2EdgeTesselator.

*   S2Coder implementations are now complete and should be interoperable
    between Go, C++ and Java.

*   New Java implementation benchmarks based on JMH, the [Java Microbenchmark
    Harness](https://github.com/openjdk/jmh).

*   Addition of S2IndexingHelper, which supports spatial indexing of regions
    in documents and finding document regions that intersect a query region.

*   Many new methods on existing classes adding features which had formerly been
    missing in Java compared to the C++ implementation.

*   Many Javadoc and other comment improvements.

*   Many implementation cleanups. Most should not cause visible changes, but
    there were some bug fixes for improved accuracy.

*   @CheckReturnValue is a package-wide default now, set in the newly added
    "package-info.java".

*   The "testdatagenerator" subpackage supports generation of test data for both
    benchmarks and unit tests. Note that the APIs in this package are not
    part of the external S2 API and are subject to change without notice.

### Breaking API changes:

*   We are now using Java 11, updating from Java 8, and beginning to use Java
    11 language features.

*   Class Matrix3x3 was renamed Matrix.

*   S2.simpleCrossing(a, b, c, d) has been removed as redundant; S2EdgeUtil
    provides the same API along with robustCrossing() and the EdgeCrosser class.

*   A bug fix for EdgeCrosser.robustCrossing() may break clients that were
    relying on robustCrossing() to detect both repeated vertices in polylines
    and self-intersecting polylines. Specifically, for degenerate edges that
    don't have shared endpoints, robustCrossing now returns -1 to indicate "No
    Crossing" while previously it would return 0. Now, a return value of 0 is
    used only to indicate that two points from different edges are the same,
    i.e. edges touch. This now behaves exactly like the equivalent
    CrossingSign() in C++ for degenerate edges.

*   Removed S2CellId.toTokenOld().

*   The nested class RangeIterator was moved from S2ShapeIndex to S2ShapeUtil.

*   The nested classes AreaMeasure, CentroidMeasure, and AreaCentroidMeasure
    were removed as they were redundant. Instead, use the methods in
    S2ShapeMeasures to compute centroids and areas for S2Shapes.

*   Support for GWT has been removed. Support for J2CL is in progress but
    incomplete. The JS API is not final.

## Disclaimer

This is not an official Google product.
