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

import static com.google.common.geometry.TestDataGenerator.generateEvenlySpacedIds;
import static java.util.concurrent.TimeUnit.NANOSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.R2Vector;
import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2LatLng;
import com.google.common.geometry.S2Point;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;
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

/** Benchmarks for {@link S2CellId}. */
public class S2CellIdBenchmark {
  private S2CellIdBenchmark() {}

  /** Benchmark state which supports repeated iteration through evenly spaced leaf cell ids. */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class EvenlySpacedLeafCellIdsState {
    private static final int NUM_IDS = 1 << 8;
    private static final int MASK_IDS = NUM_IDS - 1;

    protected S2CellId[] ids = generateEvenlySpacedIds(NUM_IDS);
    protected int index;

    @Setup(Level.Trial)
    public void setup() {
       index = 0;
    }

    /**
     * Measures the amount of time it takes to convert from leaf cell ids to non-normalized points.
     * Should be comparable to BM_ToPointRaw in the C++ implementation.
     */
    @Benchmark
    public S2Point toPointRaw() {
      return ids[index++ & MASK_IDS].toPointRaw();
    }

    /**
     * Measures the amount of time it takes to convert from leaf cell ids to points.
     * <p>Should be comparable to BM_ToPoint in the C++ implementation.
     */
    @Benchmark
    public S2Point toPoint() {
      return ids[index++ & MASK_IDS].toPoint();
    }

    /**
     * Measures the amount of time it takes to obtain the (i,j) coordinates and orientation of the
     * i,j axes for leaf cells.
     */
    @Benchmark
    public long toIJOrientation() {
      // toIJOrientation() is not public, so this benchmark uses getOrientation(), which is a single
      // bit mask operation on the toIJOrientation() result.
      return ids[index++ & MASK_IDS].getOrientation();
    }

    /** Measures the amount of time it takes to compute edge neighbors of leaf cells. */
    @Benchmark
    public S2CellId[] getEdgeNeighbors() {
      S2CellId[] neighbors = new S2CellId[4];
      ids[index++ & MASK_IDS].getEdgeNeighbors(neighbors);
      return neighbors;
    }

    /** Measures the amount of time it takes to compute vertex neighbors of leaf cells. */
    @Benchmark
    public List<S2CellId> getVertexNeighbors() {
      List<S2CellId> neighbors = new ArrayList<>();
      ids[index++ & MASK_IDS].getVertexNeighbors(S2CellId.MAX_LEVEL, neighbors);
      return neighbors;
    }

    /** Measures the amount of time it takes to compute all neighbors of leaf cells. */
    @Benchmark
    public List<S2CellId> getAllNeighbors() {
      List<S2CellId> neighbors = new ArrayList<>();
      ids[index++ & MASK_IDS].getAllNeighbors(S2CellId.MAX_LEVEL, neighbors);
      return neighbors;
    }
  }

  /**
   * State for benchmarks which iterate repeatedly through an array of 128 random cell ids at
   * a single level, for several parameterized levels.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class RandomCellIdsAtLevelState extends S2BenchmarkBaseState {
    private static final int NUM_IDS = 1 << 7;
    private static final int MASK_IDS = NUM_IDS - 1;

    @Param({"10", "20", "30"})
    int level;

    protected S2CellId[] ids;
    protected String[] tokens;
    protected int index;
    protected int anotherIndex;

    @Setup(Level.Trial)
    @Override
    public void setup() throws IOException {
      super.setup();
      ids = new S2CellId[NUM_IDS];
      tokens = new String[NUM_IDS];
      for (int i = 0; i < NUM_IDS; i++) {
        ids[i] = data.getRandomCellId(level);
        tokens[i] = ids[i].toToken();
      }

      index = 0;
      anotherIndex = 3;
    }

    /**
     * Measures the amount of time it takes to simply return the cell id.  This is essentially a
     * measure of the overhead for iterating through the S2CellIds of this benchmark state. Subtract
     * this score from other benchmark scores to get the actual time spent in the benchmarked
     * methods.
     */
    @Benchmark
    public long benchmarkStateOverhead() {
      return ids[index++ & MASK_IDS].id();
    }

    /**
     * Measures the amount of time it takes to compute the child position of a single cell id.
     */
    @Benchmark
    public int childPosition() {
      return ids[index++ & MASK_IDS].childPosition();
    }

    /**
     * Measures the performance of 1000 childPosition() calls on different ids in a tight loop.
     * Should be comparable to BM_1000_child_position in the C++ implementation.
     */
    @Benchmark
    public int childPosition1000() {
      int sum = 0;
      for (int i = 0; i < 1000; i++) {
        sum += ids[index++ & MASK_IDS].childPosition();
      }
      return sum;
    }

    /** Measures the amount of time it takes to compare one S2CellId to another. */
    @Benchmark
    public int compareTo() {
      return ids[index++ & MASK_IDS].compareTo(ids[(anotherIndex += 3) & MASK_IDS]);
    }

    /** Measures the amount of time it takes to check containment of one S2CellId in another. */
    @Benchmark
    public boolean contains() {
      return ids[index++ & MASK_IDS].contains(ids[(anotherIndex += 3) & MASK_IDS]);
    }

    /** Measures the amount of time it takes to check intersection of one S2CellId with another. */
    @Benchmark
    public boolean intersects() {
      return ids[index++ & MASK_IDS].intersects(ids[(anotherIndex += 3) & MASK_IDS]);
    }

    /** Measures the amount of time it takes to get the parent of an S2CellId. */
    @Benchmark
    public long parent() {
      return ids[index++ & MASK_IDS].parent().id();
    }

    /**
     * Measures the amount of time it takes to get the lowestOnBit for an S2CellId's id. This is a
     * heavily-used low-level method.
     */
    @Benchmark
    public long lowestOnBit() {
      return ids[index++ & MASK_IDS].lowestOnBit();
    }

    /**
     * Measures the amount of time it takes to get the (s, t) coordinates of the center of a cell.
     */
    @Benchmark
    public R2Vector getCenterST() {
      return ids[index++ & MASK_IDS].getCenterST();
    }

    /**
     * Measures the amount of time it takes to convert an S2CellId to an S2LatLng. This is a
     * convenience method currently implemented as passing toPointRaw() to the S2LatLng(S2Point)
     * constructor, so perhaps doesn't need its own benchmark, but it is heavily used.
     */
    @Benchmark
    public S2LatLng toLatLng() {
      return ids[index++ & MASK_IDS].toLatLng();
    }

    /**
     * Measures the amount of time it takes to run toToken() on random cell ids. Should be
     * comparable to BM_ToToken in the C++ implementation.
     */
    @Benchmark
    public String toToken() {
      return ids[index++ & MASK_IDS].toToken();
    }

    /**
     * Measures the amount of time it takes to run fromToken(). Should be comparable to BM_FromToken
     * in the C++ implementation.
     */
    @Benchmark
    public S2CellId fromToken() {
      return S2CellId.fromToken(tokens[index++ & MASK_IDS]);
    }

    /**
     * Measures the performance of level() in a tight loop. Should be comparable to BM_1000_level in
     * the C++ implementation.
     */
    @Benchmark
    public int level1000() {
      int sum = 0;
      for (int i = 0; i < 1000; i++) {
        // level() only depends on how many zero bits are at the end of the cell id, so we do not
        // actually need to test it on different cells at the same level, but this does anyway to
        // ensure Java can't over-optimize.
        sum += ids[index++ & MASK_IDS].level();
      }
      return sum;
    }

    /** Measures the amount of time it takes to compute the level of a single cell id. */
    @Benchmark
    public int level() {
      // level() only depends on how many zero bits are at the end of the cell id, so we do not need
      // to test it on different cells at the same level.
      return ids[1].level();
    }
  }

  /**
   * Benchmark state that supports repeatedly iterating through 20 points in a helical spiral around
   * the z-axis.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class HelicalSpiralState {
    // Number of points in the spiral.
    private static final int NUM_POINTS = 20;

    // Sample points used for benchmarking fromPoint().
    protected ArrayList<S2Point> points;
    protected int index;

    @Setup(Level.Trial)
    public void setup() throws IOException {
      points = new ArrayList<>(NUM_POINTS);
      // The sample points follow a spiral curve that completes one revolution around the z-axis
      // every 1/dt samples.  The z-coordinate increases from -4 to +4 over 'count' samples.
      S2Point p = new S2Point(1, 0, -4);
      double dz = (-2 * p.getZ()) / NUM_POINTS;
      double dt = 0.05;
      for (int r = NUM_POINTS; r > 0; --r) {
        points.add(p);
        // Cheap revolution around the z-axis
        p = p.add(new S2Point(-dt * p.getY(), dt * p.getX(), dz));
      }
      index = 0;
    }

  /**
   * Measures the cost of calling {@link S2CellId#fromPoint(S2Point)} on 20 points in a helical
   * spiral around the z-axis. Should be comparable to BM_FromPoint in the C++ implementation.
   */
    @Benchmark
    public S2CellId fromPoint() {
      S2CellId c = S2CellId.fromPoint(points.get(index));
      index = (index + 1) % NUM_POINTS;
      return c;
    }
  }

  /**
   * Benchmark state that supplies a cube face (range 0..5) and i- and j-coordinates in the range
   * [0, 2^30 - 1].
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class FaceIJState {
    protected int face;
    protected int i;
    protected int j;

    @Setup(Level.Trial)
    public void setup() {
      // These do not need to be varied, because fromFaceIJ's performance is not input dependant.
      Random random = new Random(123455);
      face = random.nextInt(6);
      // i, j are "random" integers in [0, 2^30 - 1].
      i = random.nextInt(1 << 30);
      j = random.nextInt(1 << 30);
    }

    /** Measures the amount of time it takes to return a leaf cell id with fromFaceIJ. */
    @Benchmark
    public S2CellId fromFaceIJ() {
      return S2CellId.fromFaceIJ(face, i, j);
    }
  }
}
