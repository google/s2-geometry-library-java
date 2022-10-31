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

import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.MILLISECONDS;
import static java.util.concurrent.TimeUnit.NANOSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.base.Preconditions;
import com.google.common.geometry.LittleEndianInput;
import com.google.common.geometry.LittleEndianOutput;
import com.google.common.geometry.S1Angle;
import com.google.common.geometry.S2CellId;
import com.google.common.geometry.S2Loop;
import com.google.common.geometry.S2Point;
import com.google.common.geometry.S2Polygon;
import com.google.common.geometry.S2RegionCoverer;
import com.google.common.geometry.S2Shape;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.TestDataGenerator.PolygonFactory;
import com.google.common.geometry.benchmarks.S2BenchmarkBaseState.ContainsPointBaseState;
import com.google.common.geometry.benchmarks.S2BenchmarkBaseState.PolygonAndComplementListState;
import com.google.common.geometry.benchmarks.S2BenchmarkBaseState.PolygonAndCoveringListState;
import com.google.common.geometry.benchmarks.S2BenchmarkBaseState.PolygonListState;
import com.google.errorprone.annotations.CheckReturnValue;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
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
import org.openjdk.jmh.infra.BenchmarkParams;
import org.openjdk.jmh.infra.Blackhole;
import org.openjdk.jmh.infra.Control;

/** Benchmarks for {@link S2Polygon}. */
@CheckReturnValue
public class S2PolygonBenchmark {
  private S2PolygonBenchmark() {}

  /** Constructs an S2Polygon from a list of lists of vertices, via constructing S2Loops. */
  protected static S2Polygon constructFromVertices(List<List<S2Point>> loopVertices) {
    int numLoops = loopVertices.size();

    List<S2Loop> loops = new ArrayList<>(numLoops);
    for (int i = 0; i < numLoops; ++i) {
      loops.add(new S2Loop(loopVertices.get(i)));
    }
    return new S2Polygon(loops);
  }

  /**
   * Benchmark state which provides a single S2Loop, as well as a List with one List of S2Points
   * containing the same points as the S2Loop.
   */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class ConstructorSingleLoopState extends S2BenchmarkBaseState {
    @Param({"8", "1024", "1048576"})
    int totalNumVertices;

    // The loop which will be repeatedly copied and used to construct a polygon.
    private S2Loop testLoop;

    // A list of one list of points which will be repeatedly used to construct a polygon.
    private List<List<S2Point>> testVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        totalNumVertices = 8;
      }
      super.setup();

      S1Angle radius = S1Angle.degrees(1);
      testLoop = S2Loop.makeRegularLoop(data.getRandomPoint(), radius, totalNumVertices);

      // Prepare the input for the "full constructor" benchmark.
      testVertices = new ArrayList<>();
      List<S2Point> pointList = new ArrayList<>();
      for (int j = 0; j < testLoop.numVertices(); ++j) {
        pointList.add(testLoop.vertex(j));
      }
      testVertices.add(pointList);
    }

    /**
     * Measures the time to construct an S2Polygon of a single loop. Includes the time needed to
     * make a copy of the test loop, which is necessary because the S2Polygon constructor clears its
     * input list.
     */
    @Benchmark
    public S2Polygon constructorSingleLoop(Control control) {
      List<S2Loop> loops = new ArrayList<>();
      loops.add(new S2Loop(testLoop));
      return new S2Polygon(loops);
    }

    /**
     * Measures the time to construct a single-loop polygon from a list of list of vertices,
     * including the time required to construct the S2Loop.
     */
    @Benchmark
    public S2Polygon fullConstructorSingleLoop() {
      return constructFromVertices(testVertices);
    }
  }

  /**
   * Benchmark state which provides an S2Polygon from the LoopGridPolygonFactory, as well as a list
   * of lists of points that form the same loops as the polygon.
   */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class ConstructorLoopGridState extends S2BenchmarkBaseState {
    // Use enough vertices that more than one S2ShapeIndex cell is required.
    static final int VERTICES_PER_LOOP = 20;

    @Param({"2", "64", "4096"})
    int numLoops;

    // The single test polygon that will be repeatedly reconstructed.
    private S2Polygon testPolygon;
    private List<List<S2Point>> testVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 2;
      }
      super.setup();

      testPolygon = PolygonFactory.LOOP_GRID.newPolygon(data, numLoops, VERTICES_PER_LOOP);

      // Get a copy of the test polygon loops.
      ArrayList<S2Loop> inputLoops = new ArrayList<>();
      S2Polygon copy = new S2Polygon(testPolygon);
      copy.release(inputLoops);

      testVertices = new ArrayList<>();
      for (int i = 0; i < numLoops; ++i) {
        List<S2Point> loop = new ArrayList<>();
        S2Loop inputLoop = inputLoops.get(i);
        for (int j = 0; j < inputLoop.numVertices(); ++j) {
          loop.add(inputLoop.vertex(j));
        }
        testVertices.add(loop);
      }
    }

    /**
     * Measures S2Polygon constructor performance for various sizes of grids of small loops.
     * Includes the time needed to make a copy of the input polygon and get its loops, which is
     * necessary because the S2Polygon constructor clears its input list.
     */
    @Benchmark
    public S2Polygon constructorLoopGrid() {
      S2Polygon copy = new S2Polygon(testPolygon);
      ArrayList<S2Loop> loops = new ArrayList<>();
      copy.release(loops);
      return new S2Polygon(loops);
    }

    /**
     * Measures the time to construct a polygon from a list of a list of vertices, where the polygon
     * is a grid of small loops, including the time required to construct the S2Loops.
     */
    @Benchmark
    public S2Polygon fullConstructorLoopGrid() {
      return constructFromVertices(testVertices);
    }
  }

  /**
   * Benchmark state extending PolygonListState, and providing a list of random preindexed polygons
   * from the ConcentricLoopsPolygonFactgory. Supports iterating through this list of polygons.
   */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class IsValidConcentricLoopsState extends PolygonListState {
    // TODO(eengle): Re-enable when isValid is faster.
    // @Param({"1", "2", "256", "4096"})
    @Param({"1", "2", "10"})
    int numLoops;

    @Param({"4", "16", "256", "4096", "65536"})
    int totalNumVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 1;
        totalNumVertices = 4;
      }
      setupPolygons(data, PolygonFactory.CONCENTRIC_LOOPS, numLoops, totalNumVertices, true);
    }

    /** Measures the time required to check if concentric loops polygons are valid. */
    @Benchmark
    public void isValidConcentricLoops() {
      Preconditions.checkState(currentPolygon.isValid());
      advancePolygon();
    }
  }

  /**
   * Benchmark state which supports iterating through a list of random polygons from the
   * LoopGridPolygonFactory, with a list of query points for each, sampled from the polygon bound.
   */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  public static class ContainsPointLoopGridState extends ContainsPointBaseState {
    // TODO(eengle): Re-enable when initialization is faster.
    // @Param({"4", "64", "4096"})
    @Param({"1", "2", "10"})
    int numLoops;

    @Param({"8", "1024", "262144"})
    int totalNumVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 1;
        totalNumVertices = 8;
      }
      setupPolygonsAndQueries(
          data, PolygonFactory.LOOP_GRID, numLoops, totalNumVertices, false);
    }

    /**
     * Measures the time required to run S2Polygons.contains(S2Point) on polygons which are grids of
     * regular loops, and points sampled from the polygon bounds.
     */
    @Benchmark
    public boolean containsPoint() {
      advanceQuery();
      return currentPolygon.contains(currentQuery());
    }
  }

  /**
   * Benchmark state which supports iterating through a list of nested fractal polygons, with a list
   * of query points for each, sampled from the polygon bound.
   */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class ContainsPointNestedFractalsState extends ContainsPointBaseState {
    // TODO(eengle): Re-enable when initialization is faster.
    // @Param({"2", "4", "8", "16"})
    @Param({"1", "2", "10"})
    int numLoops;

    @Param({"64", "4096", "262144"})
    int totalNumVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 1;
        totalNumVertices = 8;
      }
      setupPolygonsAndQueries(
          data, PolygonFactory.NESTED_FRACTALS, numLoops, totalNumVertices, false);
    }

    /**
     * Measures the time required to run S2Polygons.contains(S2Point) on polygons which are grids of
     * regular loops, and points sampled from the polygon bounds.
     */
    @Benchmark
    public boolean containsPoint() {
      advanceQuery();
      return currentPolygon.contains(currentQuery());
    }
  }

  /**
   * Benchmark state supporting iteration through polygons with either two nested loops or one
   * fractal loop, and an S2RegionCoverer.
   */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class CoveringState extends PolygonListState {
    @Param({"8", "64", "512"})
    int maxCells;

    @Param({"8", "128", "2048", "32768", "524288"})
    int totalNumVertices;

    @Param({"NESTED_LOOPS", "NESTED_FRACTALS"})
    PolygonFactory polygonFactory;

    private S2RegionCoverer coverer;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        maxCells = 8;
        totalNumVertices = 8;
      }
      setupPolygons(
          data,
          polygonFactory,
          (polygonFactory == PolygonFactory.NESTED_FRACTALS) ? 1 : 2,
          totalNumVertices,
          true);
      coverer = S2RegionCoverer.builder().setMaxCells(maxCells).build();
    }

    /** Measures the time to get the covering of a polygon. */
    @Benchmark
    public ArrayList<S2CellId> getCovering() {
      ArrayList<S2CellId> covering = new ArrayList<>();
      coverer.getCovering(currentPolygon, covering);
      advancePolygon();
      return covering;
    }
  }

  /** Benchmark state supporting decoding polygons. */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class DecodePolygonState extends PolygonListState {
    @Param({"1", "2"})
    int numLoops;

    @Param({"16", "64", "256", "1024", "4096"})
    int totalNumVertices;

    protected ArrayList<byte[]> uncompressedPolygonData = new ArrayList<>(NUM_POLYGONS);
    protected ArrayList<byte[]> encodedPolygonData = new ArrayList<>(NUM_POLYGONS);

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 1;
        totalNumVertices = 16;
      }
      setupPolygons(data, PolygonFactory.NESTED_LOOPS, numLoops, totalNumVertices, true);

      ByteArrayOutputStream uncompressedOut = new ByteArrayOutputStream();
      ByteArrayOutputStream out = new ByteArrayOutputStream();
      LittleEndianOutput uncompressedEncoder = new LittleEndianOutput(uncompressedOut);

      for (int i = 0; i < NUM_POLYGONS; i++) {
        // Encode each factory polygon using the internal API for the uncompressed encoding.
        polygons.get(i).encodeUncompressed(uncompressedEncoder);
        uncompressedPolygonData.add(uncompressedOut.toByteArray());
        uncompressedOut.reset();

        // Snap the factory polygon to leaf cells and encode with the public API. This will use a
        // compressed encoding.
        S2Polygon snapped = new S2Polygon();
        snapped.initToSnapped(polygons.get(i), S2CellId.MAX_LEVEL);
        snapped.encode(out);
        encodedPolygonData.add(out.toByteArray());
        out.reset();
      }
    }

    /**
     * Measures the time to decode a polygon from bytes using the uncompressed encoding.
     */
    @Benchmark
    public S2Polygon decodeUncompressed() throws IOException {
      LittleEndianInput decoder =
          new LittleEndianInput(
              new ByteArrayInputStream(uncompressedPolygonData.get(polygonIndex)));

      byte unusedVersion = decoder.readByte();
      S2Polygon decodedPolygon = S2Polygon.decodeUncompressed(decoder);
      advancePolygon();
      return decodedPolygon;
    }

    /**
     * Measures the time to decode a snapped polygon from bytes using the compressed encoding.
     */
    @Benchmark
    public S2Polygon decodeCompressed() throws IOException {
      ByteArrayInputStream inputStream =
          new ByteArrayInputStream(encodedPolygonData.get(polygonIndex));

      S2Polygon decodedPolygon = S2Polygon.decode(inputStream);
      advancePolygon();
      return decodedPolygon;
    }
  }

  /** Benchmark state supporting iteration through preindexed loop grid polygons. */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MILLISECONDS)
  public static class ContainsSelfLoopGridState extends PolygonListState {
    // TODO(eengle): Re-enable when initialization is faster.
    // @Param({"4", "128", "4096"})
    @Param({"1", "2", "10"})
    int numLoops;

    @Param({"64", "512", "4096", "32768"})
    int totalNumVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 1;
        totalNumVertices = 64;
      }
      setupPolygons(data, PolygonFactory.LOOP_GRID, numLoops, totalNumVertices, true);
    }

    /** Measures the time to check containment of a loop grid polygon with a copy of itself. */
    @Benchmark
    public void containsSelfLoopGrid() {
      Preconditions.checkState(currentPolygon.contains(currentPolygonCopy));
      advancePolygon();
    }
  }

  /** Benchmark state supporting iteration through preindexed nested fractal polygons. */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class ContainsSelfNestedFractalsState extends PolygonListState {
    // If the parameter space for this benchmark was the same as ContainsSelfLoopGridState, the
    // factory could be parameterized to produce a single ContainsSelf benchmark state.
    @Param({"2"})
    int numLoops;

    @Param({"8", "64", "512", "4096", "32768", "262144"})
    int totalNumVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 2;
        totalNumVertices = 8;
      }
      setupPolygons(data, PolygonFactory.NESTED_FRACTALS, numLoops, totalNumVertices, true);
    }

    /** Measures the time to check containment of a nested fractal polygon with a copy of itself. */
    @Benchmark
    public void containsSelfNestedFractals() {
      Preconditions.checkState(currentPolygon.contains(currentPolygonCopy));
      advancePolygon();
    }
  }

  /**
   * Benchmark state which supports iteration through preindexed loop grid polygons with their
   * complements.
   */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class IntersectsCmplLoopGridState extends PolygonAndComplementListState {
    @Param({"4", "64", "512", "4096"})
    int numLoops;

    @Param({"64", "512", "4096", "32768"})
    int totalNumVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 4;
        totalNumVertices = 64;
      }
      setupPolygons(data, PolygonFactory.LOOP_GRID, numLoops, totalNumVertices, true);
    }

    /**
     * Measures the time to intersect a preindexed polygon of loops forming a grid, with its
     * complement.
     */
    @Benchmark
    public void intersectsCmplLoopGrid() {
      Preconditions.checkState(!currentPolygon.intersects(currentPolygonComplement));
      advancePolygon();
    }
  }

  /**
   * Benchmark state supporting iteration through preindexed, nested fractal polygons with their
   * complements.
   */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class IntersectsCmplNestedFractalsState extends PolygonAndComplementListState {
    // If the parameter space for this benchmark was the same as IntersectsCmplLoopGridState, the
    // factory could be parameterized to produce a single IntersectsCmpl benchmark state.
    @Param({"2"})
    int numLoops;

    @Param({"512", "262144"})
    int totalNumVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 2;
        totalNumVertices = 512;
      }
      setupPolygons(data, PolygonFactory.NESTED_FRACTALS, numLoops, totalNumVertices, true);
    }

    /** Measures the time to intersect a polygon with its complement. */
    @Benchmark
    public void intersectsCmplNestedFractals() {
      Preconditions.checkState(!currentPolygon.intersects(currentPolygonComplement));
      advancePolygon();
    }
  }

  /**
   * Benchmark state supporting iteration through preindexed, fractal polygons of a single loop,
   * with their respective cell coverings, and through the cells of the coverings and the polygons.
   * Intersection is the only interesting operation to benchmark, since for union and difference the
   * output polygon includes almost all of the input polygon.
   */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 15, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 15, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class IntersectFractalWithCoveringState extends PolygonAndCoveringListState {
    @Param({"1", "2", "32", "128", "4096"})
    int maxCells;

    @Param({"128", "32768"})
    int totalNumVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        maxCells = 1;
        totalNumVertices = 128;
      }
      setupPolygonsAndCoverings(
          data, PolygonFactory.NESTED_FRACTALS, 1, totalNumVertices, true, maxCells);
    }

    /**
     * Measures the time required to intersect a preindexed, fractal polygon with the individual
     * cells of its covering. For each invocation, intersects another cell with the current polygon,
     * and when all cells of the current polygon's covering are done, advances to the next polygon.
     */
    @Benchmark
    public S2Polygon intersectFractalWithCovering() {
      S2Polygon result = new S2Polygon();
      result.initToIntersection(currentPolygon, constructIndexedPolygonForCurrentCell());
      advanceCovering();
      return result;
    }
  }

  /**
   * Benchmark state providing a list of loop grid polygons with their respective cell coverings,
   * and supporting iteration through the cells of the coverings and the polygons. These can be very
   * slow operations, taking many seconds each when there are a large number of loops.
   */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 30, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 30, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MILLISECONDS)
  public static class IntersectLoopGridWithCovering extends PolygonAndCoveringListState {
    static final int NUM_VERTICES_PER_LOOP = 4;

    @Param({"1", "2", "32", "128", "4096"})
    int maxCells;

    @Param({"128", "32768"})
    int totalNumVertices;

    @Setup(Level.Iteration)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        maxCells = 1;
        totalNumVertices = 128;
      }
      int numLoops = totalNumVertices / NUM_VERTICES_PER_LOOP;
      setupPolygonsAndCoverings(
          data, PolygonFactory.LOOP_GRID, numLoops, totalNumVertices, true, maxCells);
    }

    /**
     * Measures the time to initialize a new S2Polygon as the intersection of a Loop Grid polygon
     * and one of the cells of its covering.
     */
    @Benchmark
    public S2Polygon intersectLoopGridWithCovering() {
      S2Polygon result = new S2Polygon();
      result.initToIntersection(currentPolygon, constructIndexedPolygonForCurrentCell());
      advanceCovering();
      return result;
    }
  }

  /** Benchmark state supporting iteration through preindexed nested fractal polygons. */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class UnionNestedFractalWithSelfState extends PolygonListState {
    @Param({"1", "16"})
    int numLoops;

    @Param({"32768"})
    int totalNumVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 1;
        totalNumVertices = 512;
      }
      setupPolygons(data, PolygonFactory.NESTED_FRACTALS, numLoops, totalNumVertices, true);
    }

    @Benchmark
    public S2Polygon unionNestedFractalWithSelf() {
      S2Polygon result = new S2Polygon();
      result.initToUnion(currentPolygon, currentPolygonCopy);
      advancePolygon();
      return result;
    }

    @Benchmark
    public S2Polygon subtractNestedFractalFromSelf() {
      S2Polygon result = new S2Polygon();
      result.initToDifference(currentPolygon, currentPolygonCopy);
      advancePolygon();
      return result;
    }
  }

  /** Benchmark state providing iteration through preindexed loop grid polygons. */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class UnionLoopGridWithBoundState extends PolygonListState {
    @Param({"1", "8", "64", "512", "4096"})
    int numLoops;

    @Param({"32768"})
    int totalNumVertices;

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 1;
        totalNumVertices = 512;
      }
      setupPolygons(data, PolygonFactory.LOOP_GRID, numLoops, totalNumVertices, true);
    }

    /**
     * Measures the time to initialize a new S2Polygon to the union of a Loop Grid polygon and its
     * bound.
     */
    @Benchmark
    public S2Polygon unionLoopGridWithBound() {
      S2Polygon result = new S2Polygon();
      result.initToUnion(currentPolygon, currentPolygonBound);
      advancePolygon();
      return result;
    }

    /**
     * Measures the time to initialize a new S2Polygon to the difference of a Loop Grid polygon and
     * its bound.
     */
    @Benchmark
    public S2Polygon subtractBoundFromLoopGrid() {
      S2Polygon result = new S2Polygon();
      result.initToDifference(currentPolygon, currentPolygonBound);
      advancePolygon();
      return result;
    }
  }

  /** Benchmark state supporting iteration through concentric loop polygons, and a MutableEdge. */
  @State(Scope.Thread)
  @Warmup(iterations = 3, time = 5, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 5, timeUnit = SECONDS)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  public static class PolygonShapeState extends PolygonListState {
    static final int LOOP_SIZE = 6;
    static final int GET_EDGE_REPS = 1000;

    @Param({"1", "3", "4", "5", "6", "15"})
    int numLoops;

    private final MutableEdge edge = new MutableEdge();

    @Setup(Level.Trial)
    public void setup(BenchmarkParams params) throws IOException {
      super.setup();
      setUnitTestModeFromParams(params);
      if (isUnitTest()) {
        numLoops = 1;
      }
      setupPolygons(data, PolygonFactory.CONCENTRIC_LOOPS, numLoops, LOOP_SIZE, false);
    }

    /**
     * Measures the time to get an S2Shape from the current polygon, and then do GET_EDGE_REPS
     * consecutive calls to getEdge() on it.
     */
    @Benchmark
    public void getManyShapeEdgesPerPolygon(Blackhole bh) {
      S2Shape shape = currentPolygon.shape();
      int reps = GET_EDGE_REPS;
      for (int e = 0; reps >= 0; --reps) {
        if (--e < 0) {
          e += LOOP_SIZE * numLoops;
        }
        shape.getEdge(e, edge);
        bh.consume(edge);
      }
      advancePolygon();
    }
  }
}
