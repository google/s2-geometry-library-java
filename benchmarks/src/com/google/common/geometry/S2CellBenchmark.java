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

import static com.google.common.geometry.S2Projections.PROJ;
import static com.google.common.geometry.TestDataGenerator.generateEvenlySpacedIds;
import static java.util.concurrent.TimeUnit.MILLISECONDS;
import static java.util.concurrent.TimeUnit.NANOSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.R2Rect;
import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S1ChordAngle;
import com.google.common.geometry.S2Cap;
import com.google.common.geometry.S2Cell;
import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2LatLngRect;
import com.google.common.geometry.S2Point;
import java.io.IOException;
import java.util.ArrayList;
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
import org.openjdk.jmh.infra.Blackhole;

/** Benchmarks for {@link S2Cell}. */
public class S2CellBenchmark {
  private S2CellBenchmark() {}

  /**
   * Benchmark state that supports iterating through fixed arrays of S2Cells at leaf level
   * or level 18, S2CellIds, and nearby points.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class LeafCellIdsState extends S2BenchmarkBaseState {
    private static final int NUM_IDS = 1 << 8;
    private static final int MASK_IDS = NUM_IDS - 1;

    // Leaf cell ids, evenly spaced in the space of leaf cell ids
    protected S2CellId[] leafCellIds;
    // Leaf cells
    protected S2Cell[] leafCells;

    // Non leaf cell ids, randomly selected at nonLeafLevel.
    protected S2Cell[] nonLeafCells;
    private static final int NON_LEAF_LEVEL = 18; // average edge length 35m

    // Random points within about 35 meters of the centers of the non-leaf cell with the same index.
    protected S2Point[] pointA;
    protected S2Point[] pointB;

    // Temporary, allocated once and modified repeatedly.
    protected final S2Cell[] children = {
      S2Cell.fromFace(0), S2Cell.fromFace(0), S2Cell.fromFace(0), S2Cell.fromFace(0)
    };

    protected int index;

    @Setup(Level.Trial)
    @Override
    public void setup() throws IOException {
      super.setup();
      leafCellIds = generateEvenlySpacedIds(NUM_IDS);
      leafCells = new S2Cell[NUM_IDS];

      nonLeafCells = new S2Cell[NUM_IDS];
      pointA = new S2Point[NUM_IDS];
      pointB = new S2Point[NUM_IDS];

      for (int i = 0; i < NUM_IDS; i++) {
        leafCells[i] = new S2Cell(leafCellIds[i]);

        nonLeafCells[i] = new S2Cell(data.getRandomCellId(NON_LEAF_LEVEL));
        S2Cap bound = nonLeafCells[i].getCapBound();
        S2Cap sampleCap = S2Cap.fromAxisChord(
            bound.axis(), S1ChordAngle.add(bound.radius(), bound.radius()));
        pointA[i] = data.samplePoint(sampleCap);
        pointB[i] = data.samplePoint(sampleCap);
      }

      index = 0;
    }

    /**
     * Measures the amount of time it takes to construct {@link S2Cell} instances from leaf {@link
     * S2CellId} instances.
     */
    @Benchmark
    public S2Cell construct() {
      return new S2Cell(leafCellIds[index++ & MASK_IDS]);
    }

    /**
     * Measures the amount of time it takes to get corner vertices of {@link S2Cell} instances.
     */
    @Benchmark
    public void getVerticesRaw(Blackhole bh) {
      S2Cell c = leafCells[index++ & MASK_IDS];
      bh.consume(c.getVertexRaw(0));
      bh.consume(c.getVertexRaw(1));
      bh.consume(c.getVertexRaw(2));
      bh.consume(c.getVertexRaw(3));
    }

    /**
     * Measures the amount of time it takes to get the inward-facing normals of the great circles
     * passing through the edge from vertex k to vertex k+1 (mod 4), as returned from getEdgeRaw().
     */
    @Benchmark
    public void getEdgesRaw(Blackhole bh) {
      S2Cell c = leafCells[index++ & MASK_IDS];
      bh.consume(c.getEdgeRaw(0));
      bh.consume(c.getEdgeRaw(1));
      bh.consume(c.getEdgeRaw(2));
      bh.consume(c.getEdgeRaw(3));
    }

    @Benchmark
    public boolean subdivideRandom(Blackhole bh) {
      S2Cell c = nonLeafCells[index++ & MASK_IDS];
      boolean nonLeaf = c.subdivide(children);
      bh.consume(children[0].id());
      bh.consume(children[1].id());
      bh.consume(children[2].id());
      bh.consume(children[3].id());
      return nonLeaf;
    }

    /**
     * Measures the amount of time it takes to get the UV Bounds of the cell, including allocating
     * and returning the R2Rect.
     */
    @Benchmark
    public R2Rect getBoundsUV() {
      return nonLeafCells[index++ & MASK_IDS].getBoundUV();
    }

    /**
     * Measures the amount of time it takes to get the approximate area of the cell.
     */
    @Benchmark
    public double approxArea() {
      return nonLeafCells[index++ & MASK_IDS].approxArea();
    }

    /**
     * Measures the amount of time it takes to get the exact area of a cell.
     */
    @Benchmark
    public double exactArea() {
      return nonLeafCells[index++ & MASK_IDS].exactArea();
    }

    /**
     * Measures the amount of time it takes to get the S2Cap bound of a cell.
     */
    @Benchmark
    public S2Cap getCapBound() {
      return nonLeafCells[index++ & MASK_IDS].getCapBound();
    }

    /**
     * Measures the amount of time it takes to get the rectangular lat,lng bound of a cell.
     */
    @Benchmark
    public S2LatLngRect getRectBound() {
      return nonLeafCells[index++ & MASK_IDS].getRectBound();
    }

    /** Measures the time it takes to determine if a nearby point is contained by a cell. */
    @Benchmark
    public boolean containsPoint() {
      boolean result = nonLeafCells[index].contains(pointA[index]);
      index = index++ & MASK_IDS;
      return result;
    }

    /** Measures the time it takes to compute the distance between an S2Cell and a nearby point. */
    @Benchmark
    public S1ChordAngle getDistanceToPoint() {
      S1ChordAngle result = nonLeafCells[index].getDistance(pointA[index]);
      index = index++ & MASK_IDS;
      return result;
    }

    /**
     * Measures the time it takes to compute the maximum distance between an S2Cell and a nearby
     * point.
     */
    @Benchmark
    public S1ChordAngle getMaxDistanceToPoint() {
      S1ChordAngle result = nonLeafCells[index].getMaxDistance(pointA[index]);
      index = index++ & MASK_IDS;
      return result;
    }

    /**
     * Measures the amount of time it takes to get the distance to an edge of an S2Cell from a
     * nearby edge.
     */
    @Benchmark
    public S1ChordAngle getDistanceToEdge() {
      S1ChordAngle result = nonLeafCells[index].getDistanceToEdge(pointA[index], pointB[index]);
      index = index++ & MASK_IDS;
      return result;
    }

    /**
     * Measures the amount of time it takes to get the maximum distance to an edge of an S2Cell from
     * a nearby edge.
     */
    @Benchmark
    public S1ChordAngle getMaxDistanceToEdge() {
      S1ChordAngle result = nonLeafCells[index].getMaxDistance(pointA[index], pointB[index]);
      index = index++ & MASK_IDS;
      return result;
    }

    /**
     * Measures the amount of time it takes to get the distance from one S2Cell to another, both
     * randomly selected and at the same level.
     */
    @Benchmark
    public S1ChordAngle getDistanceToCell() {
      S2Cell a = nonLeafCells[index];
      index = index++ & MASK_IDS;
      return a.getDistance(nonLeafCells[index]);
    }

    /**
     * Measures the amount of time it takes to get the maximum distance from one S2Cell to another,
     * both randomly selected and at the same level.
     */
    @Benchmark
    public S1ChordAngle getMaxDistanceToCell() {
      S2Cell a = nonLeafCells[index];
      index = index++ & MASK_IDS;
      return a.getMaxDistance(nonLeafCells[index]);
    }
  }

  /** Benchmark state with storage for four child S2Cells at each of multiple levels. */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MILLISECONDS)
  @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
  public static class CellChildrenState {
    // Every level deeper takes about four times longer, as expected.
    @Param({"11"})
    int maxSubdivisionLevel;

    // Pre-allocated space for the subdivided children at each level.
    private S2Cell[][] children;

    @Setup(Level.Trial)
    public void setup() {
      children = new S2Cell[maxSubdivisionLevel][4];
      S2CellId id = new S2CellId(0);
      for (int level = 0; level < maxSubdivisionLevel; ++level) {
        children[level][0] = new S2Cell(id);
        children[level][1] = new S2Cell(id);
        children[level][2] = new S2Cell(id);
        children[level][3] = new S2Cell(id);
      }
    }

    private void expandChildren1(S2Cell cell, int level, Blackhole bh) {
      bh.consume(cell.subdivide(children[level]));
      if (children[level][0].level() < maxSubdivisionLevel) {
        for (int pos = 0; pos < 4; ++pos) {
          expandChildren1(children[level][pos], level + 1, bh);
        }
      }
    }

    private S2CellId expandChildren2(S2Cell cell, int level, Blackhole bh) {
      S2CellId id = cell.id().childBegin();
      for (int pos = 0; pos < 4; ++pos, id = id.next()) {
        S2Cell child = new S2Cell(id);
        if (child.level() < maxSubdivisionLevel) {
          bh.consume(expandChildren2(child, level + 1, bh));
        }
      }
      return id;
    }

    /**
     * Measures the amount of time it takes to recursively subdivide a cell, starting from a face
     * and descending to maxSubdivisionLevel. Equivalent to BM_Subdivide in the C++ implementation.
     */
    @Benchmark
    public void subdivide(Blackhole bh) {
      S2Cell cell = S2Cell.fromFace(0);
      expandChildren1(cell, 0, bh);
    }

    /**
     * Measures the amount of time it takes to recursively construct S2Cells, starting from a face
     * and descending to maxSubdivisionLevel. Equivalent to BM_Constructor in the C++
     * implementation.
     */
    @Benchmark
    public void constructor(Blackhole bh) {
      S2Cell cell = S2Cell.fromFace(0);
      bh.consume(expandChildren2(cell, 0, bh));
    }
  }

  /**
   * Benchmark state with preconfigured lists of S2Cells and S2Points designed for benchmarking the
   * distance methods on S2Cell.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  public static class GetDistanceState extends S2BenchmarkBaseState {
    private static final int NUM_CASES = 5000;

    private final ArrayList<S2Cell> cells1 = new ArrayList<>();
    private final ArrayList<S2Cell> cells2 = new ArrayList<>();
    private final ArrayList<S2Point> points1 = new ArrayList<>();
    private final ArrayList<S2Point> points2 = new ArrayList<>();

    private int index;

    /**
     * Each test case consists of two random S2CellIDs and two S2Points somewhere in the vicinity of
     * the first S2CellId.
     */
    @Setup(Level.Trial)
    @Override
    public void setup() throws IOException {
      super.setup();
      for (int i = 0; i < NUM_CASES; ++i) {
        S2CellId cellId = data.getRandomCellId();
        cells1.add(new S2Cell(cellId));
        cells2.add(new S2Cell(data.getRandomCellId()));
        S1Angle radius = S1Angle.radians(2 * PROJ.maxDiag.getValue(cellId.level()));
        S2Cap cap = S2Cap.fromAxisAngle(cellId.toPoint(), radius);
        points1.add(data.samplePoint(cap));
        points2.add(data.samplePoint(cap));
      }
      index = 0;
    }

    /**
     * Measures the amount of time it takes to compute the distance from a cell to a nearby edge.
     */
    @Benchmark
    public void getDistanceToEdge(Blackhole bh) {
      bh.consume(cells1.get(index).getDistanceToEdge(points1.get(index), points2.get(index)));
      index = (index + 1) % NUM_CASES;
    }

    /** Measures the amount of time it takes to compute the distance from one cell to another. */
    @Benchmark
    public void getDistanceToCell(Blackhole bh) {
      bh.consume(cells1.get(index).getDistance(cells2.get(index)));
      index = (index + 1) % NUM_CASES;
    }
  }
}
