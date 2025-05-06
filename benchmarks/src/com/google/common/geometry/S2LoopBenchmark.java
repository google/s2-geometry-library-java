/*
 * Copyright 2021 Google Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package com.google.common.geometry.benchmarks;

import static com.google.common.geometry.TestDataGenerator.DEFAULT_LOOP_RADIUS;
import static com.google.common.geometry.TestDataGenerator.DEFAULT_NESTED_LOOP_GAP;
import static com.google.common.geometry.TestDataGenerator.indexLoop;
import static com.google.common.geometry.TestDataGenerator.makeDisjointLoopPair;
import static com.google.common.geometry.TestDataGenerator.makeNestedLoopPair;
import static com.google.common.geometry.TestDataGenerator.metersToAngle;
import static java.lang.Math.scalb;
import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2LatLngRect;
import com.google.common.geometry.S2Loop;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2ShapeIndex;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.List;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

/** Benchmarks for {@link S2Loop}. */
@CheckReturnValue
public class S2LoopBenchmark {
  private S2LoopBenchmark() {}

  /** Set the given array to the indexed loops in the given list. */
  private static void setIndexedLoops(S2Loop[] results, List<S2Loop> loops) {
    results[0] = indexLoop(loops.get(0));
    results[1] = indexLoop(loops.get(1));
  }

  /**
   * Benchmark state consisting of a circular loop of DEFAULT_RADIUS at a "random" center point,
   * which is always the same as the random number generator is reset in setup().
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class LoopState extends S2BenchmarkBaseState {
    @Param({"4", "8", "64", "512", "4096", "32768", "262144"})
    int numVertices;

    private S2Loop loop;

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      loop =
          indexLoop(
              S2Loop.makeRegularLoop(data.getRandomPoint(), DEFAULT_LOOP_RADIUS, numVertices));
    }

    /** Measures the average time to construct an S2Loop from a list of vertices. */
    @Benchmark
    public S2Loop constructor() {
      return new S2Loop(loop.vertices());
    }

    /** Measures the average time to check validity of a loop. */
    @Benchmark
    public boolean isValid() {
      return loop.isValid();
    }

    /**
     * Measures the average time to create an S2ShapeIndex, add a single loop to it, and update the
     * index.
     */
    @Benchmark
    public boolean index() {
      S2ShapeIndex index = new S2ShapeIndex();
      index.add(loop);
      index.applyUpdates();
      return index.isFresh();
    }
  }

  /**
   * Benchmark state consisting of a circular loop of DEFAULT_LOOP_RADIUS and @Param numVertices
   * points, as well as an array of points sampled from the rectangular bound of the loop. The
   * "random" center point and "randomly" sampled points are actually the same each time.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class LoopWithSampledPointsState extends S2BenchmarkBaseState {
    @Param({"4", "8", "64", "512", "4096", "32768", "262144"})
    int numVertices;

    /** Number of sample points to use. */
    private static final int NUM_QUERIES_PER_LOOP = 64;

    private S2Loop loop;
    private S2Point[] samplePoints = new S2Point[NUM_QUERIES_PER_LOOP];
    private int index;

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      loop = indexLoop(
          S2Loop.makeRegularLoop(data.getRandomPoint(), DEFAULT_LOOP_RADIUS, numVertices));

      S2LatLngRect bound = loop.getRectBound();
      samplePoints = new S2Point[NUM_QUERIES_PER_LOOP];
      for (int i = 0; i < NUM_QUERIES_PER_LOOP; i++) {
        samplePoints[i] = data.samplePoint(bound);
      }
      index = 0;
    }

    /** Looping iteration through the sample points. */
    public S2Point nextSamplePoint() {
      index = (index + 1) % NUM_QUERIES_PER_LOOP;
      return samplePoints[index];
    }

    /**
     * Measures the average time to compute contains(S2Point), rotating continuously through the
     * samplePoints for a mix of in/out results.
     */
    @Benchmark
    public boolean containsPoint() {
      return loop.contains(nextSamplePoint());
    }
  }

  /**
   * Benchmark state providing a pair of nested loops. Both are centered at the same "random" point,
   * with @Param 'numVertices' vertices. The outer loop has DEFAULT_RADIUS radius, and the inner
   * loop is inset by DEFAULT_NESTED_LOOP_GAP times the edge length of the outer loop. If the edge
   * length is not small enough, this would cause the inner radius to be negative, so it is at
   * clamped to at least 1% of the outer radius. The random number generator is reset in setup(), so
   * the center point is always the same.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class NestedLoopPairState extends S2BenchmarkBaseState {
    @Param({"4", "8", "64", "512", "4096", "32768", "262144"})
    int numVertices;

    private final S2Loop[] nested = new S2Loop[2];

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      setIndexedLoops(
          nested,
          makeNestedLoopPair(
              DEFAULT_LOOP_RADIUS, DEFAULT_NESTED_LOOP_GAP, numVertices, data.getRandomPoint()));
    }

    /** Measures the average time to compute compareBoundary() for nested loops. */
    @Benchmark
    public int compareBoundaryNested() {
      return nested[0].compareBoundary(nested[1]);
    }

    /** Measures the average time to compute contains() for nested loops. */
    @Benchmark
    public boolean containsNested() {
      return nested[0].contains(nested[1]);
    }

    /** Measures the average time to compute intersects() for nested loops. */
    @Benchmark
    public boolean intersectsNested() {
      return nested[0].intersects(nested[1]);
    }
  }

  /**
   * Benchmark state providing a pair of disjoint loops. The outer one is shaped somewhat like the
   * outline of the letter "C", while the inner one is nested inside the C. The center point is
   * "random" but always the same, as the random number generator is reset in setup().
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class DisjointLoopPairState extends S2BenchmarkBaseState {
    @Param({"4", "8", "64", "512", "4096", "32768", "262144"})
    int numVertices;

    private final S2Loop[] disjoint = new S2Loop[2];

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      setIndexedLoops(
          disjoint,
          makeDisjointLoopPair(
              DEFAULT_LOOP_RADIUS, DEFAULT_NESTED_LOOP_GAP, numVertices, data.getRandomPoint()));
    }

    /** Measures the average time to compute compareBoundary() for disjoint loops. */
    @Benchmark
    public int compareBoundaryDisjoint() {
      return disjoint[0].compareBoundary(disjoint[1]);
    }

    /** Measures the average time to compute contains() for disjoint loops. */
    @Benchmark
    public boolean containsDisjoint() {
      return disjoint[0].contains(disjoint[1]);
    }

    /** Measures the average time to compute intersects() for disjoint loops. */
    @Benchmark
    public boolean intersectsDisjoint() {
      return disjoint[0].intersects(disjoint[1]);
    }
  }

  /**
   * Benchmark state providing a pair of crossing loops. The first loop will have its center at a
   * "random" point 'a' with DEFAULT_RADIUS, while the second will have radius 10 times larger, and
   * a center along the arc between 'a' and another random point 'b'. Both loops have @Param
   * numVertices vertices. The random number generator is reset in setup() so the loop locations are
   * always the same.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class CrossingLoopPairState extends S2BenchmarkBaseState {
    @Param({"4096"})
    int numVertices;

    // A pair of loops.
    private final S2Loop[] crosses = new S2Loop[2];

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      setIndexedLoops(
          crosses,
          data.makeCrossingLoopPairDefault(
              numVertices, data.getRandomPoint(), data.getRandomPoint()));
    }

    /** Measures the average time to compute compareBoundary() for crossing loops. */
    @Benchmark
    public int compareBoundaryCrosses() {
      return crosses[0].compareBoundary(crosses[1]);
    }

    /** Measures the average time to compute contains() for crossing loops. */
    @Benchmark
    public boolean containsCrosses() {
      return crosses[0].contains(crosses[1]);
    }

    /** Measures the average time to compute intersects() for crossing loops. */
    @Benchmark
    public boolean intersectsCrosses() {
      return crosses[0].intersects(crosses[1]);
    }
  }

  /**
   * Benchmark state providing a pair of nested loops. Both are centered at the same "random" point,
   * with @Param numVertices vertices. The outer loop has DEFAULT_RADIUS radius, and the inner loop
   * is inset by @Param gapRatio times the edge length of the outer loop. If the edge length is not
   * small enough, this could cause the inner radius to be negative, so it is at clamped to at least
   * 1% of the outer radius. The random number generator is reset in setup(), so the center point is
   * always the same.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class NestedVariableGapLoopState extends S2BenchmarkBaseState {
    @Param({"4", "8", "64", "512", "4096", "32768", "262144"})
    int numVertices;

    @Param({"0", "1", "8", "64"})
    int gapRatio;

    // A pair of loops.
    private final S2Loop[] nestedGap = new S2Loop[2];

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      setIndexedLoops(
          nestedGap,
          makeNestedLoopPair(DEFAULT_LOOP_RADIUS, gapRatio, numVertices, data.getRandomPoint()));
    }

    /** Measures the average time to compute contains() for nested loops with various gaps. */
    @Benchmark
    public boolean containsNestedGap() {
      return nestedGap[0].contains(nestedGap[1]);
    }
  }

  /**
   * Benchmark state providing a pair of nested loops with variable radii. Both are centered at the
   * same "random" point, with @Param numVertices vertices. The outer loop has @Param meters radius,
   * and the inner loop is inset by DEFAULT_NESTED_LOOP_GAP times the edge length of the outer loop.
   * If the edge length is not small enough, this would cause the inner radius to be negative, so it
   * is at least 1% of the outer radius. The random number generator is reset in setup(), so the
   * center point is always the same.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class NestedVariableRadiiLoopState extends S2BenchmarkBaseState {
    @Param({"4", "8", "64", "512", "4096", "32768", "262144"})
    int numVertices;

    @Param({"1", "8", "64", "1024", "65535"})
    int meters; // Loop radii

    // A pair of loops.
    private final S2Loop[] loops = new S2Loop[2];

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      setIndexedLoops(
          loops,
          makeNestedLoopPair(
              metersToAngle(meters), DEFAULT_NESTED_LOOP_GAP, numVertices, data.getRandomPoint()));
    }

    /** Measures the average time to compute contains() for nested loops with various radii. */
    @Benchmark
    public boolean containsVsRadius() {
      return loops[0].contains(loops[1]);
    }
  }

  /**
   * Benchmark state providing two crossing loops, with variable log ratios in their radii. The
   * first loop has radius DEFAULT_RADIUS and a "random" center at 'a'. The second has its center
   * along the arc from 'a' to a different random point, and will have radius DEFAULT_RADIUS scaled
   * by 2 ^ @Param radiusScaleFactor. Both have @Param numVertices vertices. The random number
   * generator is reset in setup(), so the loops are always at the same locations.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class CrossingVariableLogRatiosLoopState extends S2BenchmarkBaseState {
    @Param({"4096"})
    int numVertices;

    /** The second loop radius is scaled by 2 to the following scale factors: */
    @Param({"-10", "-5", "0", "5", "10"})
    int radiusScaleFactor;

    // A pair of loops.
    private final S2Loop[] loops = new S2Loop[2];

    /** Returns the given log scale of {@link DEFAULT_RADIUS}. */
    private static S1Angle scaledDefaultRadius(int log) {
      return S1Angle.radians(scalb(DEFAULT_LOOP_RADIUS.radians(), log));
    }

    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();
      setIndexedLoops(
          loops,
          data.makeCrossingLoopPair(
              DEFAULT_LOOP_RADIUS,
              scaledDefaultRadius(radiusScaleFactor),
              numVertices,
              data.getRandomPoint(),
              data.getRandomPoint()));
    }

    /**
     * Measures the average time to compute intersects() on crossing loops as a function of the
     * relative sizes of the two loops (e.g., one loop radius much larger than the other).
     * Performance of spatial indexing can degrade when one loop is much larger than the other.
     */
    @Benchmark
    public boolean intersectsCrossesLogNeg10() {
      return loops[0].intersects(loops[1]);
    }
  }
}
