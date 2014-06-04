package com.google.common.geometry;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.caliper.Benchmark;
import com.google.caliper.Param;
import com.google.common.collect.Lists;

import java.util.ArrayList;
import java.util.List;

/**
 * A benchmark for {@link S2Polygon}. Caliper ignores returned values, but these are computed and
 * returned so that corresponding for loops don't get optimized away.
 */
public class S2PolygonBenchmark {
  /** 
   * Number of random polygons or polygon pairs to use when averaging over various positions on the
   * sphere.
   */
  private static final int NUM_POLYGON_SAMPLES = 10;

  /** Default radius for loops. */
  private static final double DEFAULT_RADIUS_KM = 10.0;

  /** Default dimension for fractal loops. */
  private static final double LOOP_FRACTAL_DIMENSION = 1.2;

  /** Spacing between adjacent loops (e.g., in a grid), as a fraction of the loop radius. */
  private static final double LOOP_SEPARATION_FRACTION = 0.1;

  static class ConstructorSingleLoop {
    @Param({"8", "1024", "1048576"})
    int numVertices;

    GeometryTestCase testUtils;

    void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /** Measures constructor performance for single loops of different sizes. */
    @Benchmark void constructorSingleLoop(int reps) {
      setUp();
      S2Polygon input = constructorBenchmarkSingleLoop(numVertices, testUtils);
      benchmarkConstructor(reps, input);
    }

    /** Measures the total constructor time, including S2Loop construction. */
    @Benchmark void fullConstructorSingleLoop(int reps) {
      setUp();
      S2Polygon input = constructorBenchmarkSingleLoop(numVertices, testUtils);
      benchmarkFullConstructor(reps, input);
    }
  }

  static class ConstructorLoopGrid {
    @Param({"2", "64", "4096"})
    int numLoops;

    GeometryTestCase testUtils;

    void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /** Measures performance of constructing polygons from an NxN grid of small loops. */
    @Benchmark void constructorLoopGrid(int reps) {
      setUp();
      S2Polygon input = constructorBenchmarkLoopGrid(numLoops, testUtils);
      benchmarkConstructor(reps, input);
    }

    /** Measures the total constructor time, including S2Loop construction. */
    @Benchmark void fullConstructorLoopGrid(int reps) {
      setUp();
      S2Polygon input = constructorBenchmarkLoopGrid(numLoops, testUtils);
      benchmarkFullConstructor(reps, input);
    }
  }

  static class IsValidConcentricLoops {
    // TODO(eengle): Re-enable when isValid is faster.
    // @Param({"1", "2", "256", "4096"})
    @Param({"1", "2", "10"})
    int numLoops;

    @Param({"4", "16", "256", "4096", "65536"})
    int totalNumVertices;

    @Benchmark void isValidConcentricLoops(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      int numVertices = Math.max(totalNumVertices / numLoops, 3);
      S2Polygon polygon = testUtils.concentricLoopsPolygon(testUtils.randomPoint(), numLoops,
          numVertices);
      for (int r = reps; r > 0; --r) {
        assertTrue(polygon.isValid());
      }
    }
  }

  static class ContainsPointLoopGrid {
    // TODO(eengle): Re-enable when initialization is faster.
    // @Param({"4", "64", "4096"})
    @Param({"1", "2", "10"})
    int numLoops;

    @Param({"8", "1024", "262144"})
    int totalNumVertices;

    /** Benchmarks containsPoint() on an NxN grid of regular loops. */
    @Benchmark void containsPointLoopGrid(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      // Ensure that each loop is non-degenerate, i.e. has at least 3 vertices.
      totalNumVertices = Math.max(totalNumVertices, 3 * numLoops);
      LoopGridPolygonFactory factory = new LoopGridPolygonFactory(numLoops, totalNumVertices);
      benchmarkContainsPoint(reps, factory, testUtils);
    }
  }

  static class ContainsPointNestedFractals {
    // TODO(eengle): Re-enable when initialization is faster.
    // @Param({"2", "4", "8", "16"})
    @Param({"1", "2", "10"})
    int numLoops;

    @Param({"64", "4096", "262144"})
    int totalNumVertices;

    /** Benchmarks containsPoint() on nested fractal loops. */
    @Benchmark void containsPointNestedFractals(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      // Ensure that each loop is non-degenerate, i.e. has at least 3 vertices.
      NestedFractalsPolygonFactory factory =
          new NestedFractalsPolygonFactory(numLoops, totalNumVertices);
      benchmarkContainsPoint(reps, factory, testUtils);
    }
  }

  /** Mainly measures contains(S2Cell) and MayIntersect(S2Cell). */
  protected static void benchmarkCovering(int reps, int maxCells, PolygonFactory factory,
      GeometryTestCase testUtils) {
    S2RegionCoverer coverer = new S2RegionCoverer();
    coverer.setMaxCells(maxCells);
    S2Polygon polygon = new S2Polygon();
    // Bresenham-type algorithm for polygon sampling.
    int delta = 0;
    for (int r = reps; r > 0; --r) {
      delta -= NUM_POLYGON_SAMPLES;
      if (delta < 0) {
        delta += reps;
        polygon = new S2Polygon(factory.newPolygon(testUtils));
      }

      ArrayList<S2CellId> covering = new ArrayList<>();
      coverer.getCovering(polygon, covering);
    }
  }

  static class Covering {
    @Param({"8", "64", "512"})
    int maxCells;

    @Param({"8", "128", "2048", "32768", "524288"})
    int numVertices;

    GeometryTestCase testUtils;

    void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    @Benchmark void fractalCovering(int reps) {
      setUp();
      NestedFractalsPolygonFactory factory = new NestedFractalsPolygonFactory(1, numVertices);
      benchmarkCovering(reps, maxCells, factory, testUtils);
    }
    @Benchmark void annulusCovering(int reps) {
      setUp();
      NestedLoopsPolygonFactory factory = new NestedLoopsPolygonFactory(2, numVertices);
      benchmarkCovering(reps, maxCells, factory, testUtils);
    }
  }

  static class ContainsSelfLoopGrid {
    // TODO(eengle): Re-enable when initialization is faster.
    // @Param({"4", "128", "4096"})
    @Param({"1", "2", "10"})
    int numLoops;

    @Param({"64", "512", "4096", "32768"})
    int totalNumVertices;

    @Benchmark void containsSelfLoopGrid(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      LoopGridPolygonFactory factory = new LoopGridPolygonFactory(numLoops, totalNumVertices);
      benchmarkContainsSelf(reps, factory, testUtils);
    }
  }

  static class ContainsSelfNestedFractals {
    @Param({"2"})
    int numLoops;

    @Param({"8", "64", "512", "4096", "32768", "262144"})
    int totalNumVertices;

    @Benchmark void containsSelfNestedFractals(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      NestedFractalsPolygonFactory factory =
          new NestedFractalsPolygonFactory(numLoops, totalNumVertices);
      benchmarkContainsSelf(reps, factory, testUtils);
    }
  }

  static class IntersectsCmplLoopGrid {
    // TODO(eengle): Re-enable when initialization is faster.
    // @Param({"4", "64", "512", "4096"})
    @Param({"1", "2", "10"})
    int numLoops;

    @Param({"64", "512", "4096", "32768"})
    int totalNumVertices;

    /**
     * Benchmarks whether a loop grid polygon intersects() its complement (obtained by subtracting
     * the polygon from a simple bound).
     */
    @Benchmark void intersectsCmplLoopGrid(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      LoopGridPolygonFactory factory = new LoopGridPolygonFactory(numLoops, totalNumVertices);
      benchmarkIntersectsCmpl(reps, factory, testUtils);
    }
  }

  static class IntersectsCmplNestedFractals {
    @Param({"2"})
    int numLoops;

    @Param({"512", "262144"})
    int totalNumVertices;

    /** Similar to ContainsSelfNestedFractals, so we just do some spot checks. */
    @Benchmark void intersectsCmplNestedFractals(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      NestedFractalsPolygonFactory factory =
          new NestedFractalsPolygonFactory(numLoops, totalNumVertices);
      benchmarkIntersectsCmpl(reps, factory, testUtils);
    }
  }

  private interface PolygonOperation {
    S2Polygon operation(S2Polygon a, S2Polygon b);
  }

  private static class InitToIntersection implements PolygonOperation {
    @Override public S2Polygon operation(S2Polygon a, S2Polygon b) {
      S2Polygon c = new S2Polygon();
      c.initToIntersection(a,  b);
      return c;
    }
  }

  private static class InitToUnion implements PolygonOperation {
    @Override public S2Polygon operation (S2Polygon a, S2Polygon b) {
      S2Polygon c = new S2Polygon();
      c.initToUnion(a, b);
      return c;
    }
  }

  private static class InitToDifference implements PolygonOperation {
    @Override public S2Polygon operation (S2Polygon a, S2Polygon b) {
      S2Polygon c = new S2Polygon();
      c.initToDifference(a, b);
      return c;
    }
  }

  static class IntersectFractalWithCovering {
    @Param({"1", "2", "32", "128", "4096"})
    int maxCells;

    @Param({"128", "32768"})
    int totalNumVertices;

    GeometryTestCase testUtils;

    void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /** 
     * The only interesting case is intersection, since for union and difference the output polygon
     * includes almost all of the input polygon.
     */
    @Benchmark void intersetFractalWithCovering(int reps) {
      setUp();
      NestedFractalsPolygonFactory factory = new NestedFractalsPolygonFactory(1, totalNumVertices);
      benchmarkOpWithCovering(reps, maxCells, factory, new InitToIntersection(), testUtils);
    }

    @Benchmark void intersectLoopGridWithCovering(int reps) {
      setUp();
      int numVerticesPerLoop = 4;
      int numLoops = totalNumVertices / numVerticesPerLoop;
      LoopGridPolygonFactory factory = new LoopGridPolygonFactory(numLoops, totalNumVertices);
      benchmarkOpWithCovering(reps, maxCells, factory, new InitToIntersection(), testUtils);
    }
  }



  static class UnionNestedFractalWithSelf {
    @Param({"1", "16"})
    int numLoops;

    @Param({"32768"})
    int totalNumVertices;

    GeometryTestCase testUtils;

    void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    @Benchmark void unionNestedFractalWithSelf(int reps) {
      setUp();
      NestedFractalsPolygonFactory factory 
      = new NestedFractalsPolygonFactory(numLoops, totalNumVertices);
      benchmarkOpWithSelf(reps, factory, new InitToUnion(), testUtils);
    }

    @Benchmark void subtractNestedFractalFromSelf(int reps) {
      setUp();
      NestedFractalsPolygonFactory factory =
          new NestedFractalsPolygonFactory(numLoops, totalNumVertices);
      benchmarkOpWithSelf(reps, factory, new InitToDifference(), testUtils);
    }
  }

  static class UnionLoopGridWithBound {
    @Param({"1", "8", "64", "512", "4096"})
    int numLoops;

    @Param({"32768"})
    int totalNumVertices;

    GeometryTestCase testUtils;

    void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    @Benchmark void unionLoopGridWithBound(int reps) {
      setUp();
      LoopGridPolygonFactory factory = new LoopGridPolygonFactory(numLoops, totalNumVertices);
      benchmarkOpWithBound(reps, factory, new InitToUnion(), testUtils);
    }

    @Benchmark void subtractBoundFromLoopGrid(int reps) {
      setUp();
      LoopGridPolygonFactory factory = new LoopGridPolygonFactory(numLoops, totalNumVertices);
      benchmarkOpWithBound(reps, factory, new InitToDifference(), testUtils);
    }
  }

  /** A factory of a particular kind of S2Polygon to benchmark. */
  private interface PolygonFactory {
    S2Polygon newPolygon(GeometryTestCase testUtils);
  }

  /** 
   * Constructs a polygon consisting of nested regular loops.  Note that 'numLoops' should be
   * fairly small (<20), otherwise the edges can get so small that loops become invalid.  (This is
   * checked.)
   */
  private static class NestedLoopsPolygonFactory implements PolygonFactory {
    private int numLoops;
    private int totalNumVertices;

    NestedLoopsPolygonFactory(int numLoops, int totalNumVertices) {
      this.numLoops = numLoops;
      this.totalNumVertices = totalNumVertices;
    }

    /** Gets new S2Polygon consisting of nested regular loops. */
    @Override
    public S2Polygon newPolygon(GeometryTestCase testUtils) {
      List<S2Loop> loops = Lists.newArrayList();
      int numVerticesPerLoop = Math.max(3, totalNumVertices / numLoops);
      // Each loop is smaller than the previous one by a fixed ratio.  The first factor below
      // compensates for loop edges cutting across the interior of the circle, while the second
      // accounts for the separation requested by the user (which is measured relative to the radius
      // of the inner loop.
      double scale = Math.cos(Math.PI / numVerticesPerLoop) / (1 + LOOP_SEPARATION_FRACTION);
      S2Point center = testUtils.randomPoint();
      S1Angle radius = GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM);
      for (int i = 0; i < numLoops; ++i) {
        S2Loop loop = S2Loop.makeRegularLoop(center, radius, totalNumVertices);
        checkEdgeLength(loop);
        loops.add(loop);
        radius = S1Angle.radians(radius.radians() * scale);
      }

      return new S2Polygon(loops);
    }
  }

  /** 
   * Constructs a polygon consisting of nested fractal loops.  Note that 'numLoops' should be
   * fairly small (<20), otherwise the edges can get so small that loops become invalid.  (This is
   * checked.)
   */
  private static class NestedFractalsPolygonFactory implements PolygonFactory {
    private int numLoops;
    private int totalNumVertices;

    NestedFractalsPolygonFactory(int numLoops, int totalNumVertices) {
      this.numLoops = numLoops;
      this.totalNumVertices = totalNumVertices;
    }

    /** Gets new S2Polgyon consisting of nested fractals. */
    @Override
    public S2Polygon newPolygon(GeometryTestCase testUtils) {
      List<S2Loop> loops = Lists.newArrayList();
      int numVerticesPerLoop = Math.max(3, totalNumVertices / numLoops);
      S2FractalBuilder fractal = new S2FractalBuilder(testUtils.rand);
      fractal.setFractalDimension(LOOP_FRACTAL_DIMENSION);
      fractal.setLevelForApproxMaxEdges(numVerticesPerLoop);

      // Scale each loop so that it fits entirely within the previous loop, with a gap equal to some
      // fraction of the inner loop's radius.
      double scale = (fractal.minRadiusFactor() / fractal.maxRadiusFactor()
          / (1 + LOOP_SEPARATION_FRACTION));
      S2Point center = testUtils.randomPoint();
      S1Angle radius = GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM);
      for (int i = 0; i < numLoops; ++i) {
        S2Loop loop = fractal.makeLoop(Matrix3x3.fromCols(S2.getFrame(center)), radius);
        checkEdgeLength(loop);
        loops.add(loop);
        radius = S1Angle.radians(radius.radians() * scale);
      }
      return new S2Polygon(loops);
    }
  }

  /** Constructs a polygon consisting of a grid of regular loops. */
  private static class LoopGridPolygonFactory implements PolygonFactory {
    private int numLoops;
    private int totalNumVertices;

    LoopGridPolygonFactory(int numLoops, int totalNumVertices) {
      this.numLoops = numLoops;
      this.totalNumVertices = totalNumVertices;
    }

    /** Gets new S2Polygon consisting of a grid of regular loops. */
    @Override
    public S2Polygon newPolygon(GeometryTestCase testUtils) {
      // Reduce the loop radius if necessary to limit distortion caused by the spherical projection
      // (the total diameter of the grid is clamped to one-quarter of the sphere).  We then
      // reduce the loop radius by the maximum amount of distortion at the edges of the grid to
      // ensure that the loops do not intersect.
      int sqrtNumLoops = (int) Math.ceil(Math.sqrt(numLoops));
      double spacingMultiplier = 1 + LOOP_SEPARATION_FRACTION;
      double maxAngle = Math.min(Math.PI / 4, sqrtNumLoops * spacingMultiplier
          * GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM).radians());
      double spacing = 2 * maxAngle / sqrtNumLoops;
      double radius = 0.5 * spacing * Math.cos(maxAngle) / spacingMultiplier;

      // If 'numLoops' is not a perfect square, make the grid slightly larger and leave some
      // locations empty.
      int left = numLoops;
      int numVerticesPerLoop = Math.max(3, totalNumVertices / numLoops);
      List<S2Loop> loops = Lists.newArrayList();
      for (int i = 0; i < sqrtNumLoops && left > 0; ++i) {
        for (int j = 0; j < sqrtNumLoops && left > 0; ++j, --left) {
          S2Point center = new S2Point(Math.tan((i + 0.5) * spacing - maxAngle), 
              Math.tan((j + 0.5) * spacing - maxAngle), 1.0);
          S2Loop loop = S2Loop.makeRegularLoop(center, S1Angle.radians(radius),
              numVerticesPerLoop);
          checkEdgeLength(loop);
          loops.add(loop);
        }
      }

      return new S2Polygon(loops);
    }
  }

  /** 
   * Since our indexing is based on S2CellIds, it becomes inefficient when loop edges are much
   * smaller than the leaf cell size (which is about 1cm).  Since edges smaller than 1mm are not
   * typically needed for geographic data, we have a sanity check so that we don't accidentally
   * generate such data in the benchmarks.
   */
  private static void checkEdgeLength(S2Loop loop) {
    S1Angle minEdgeLength = GeometryTestCase.metersToAngle(0.001);
    for (int j = 0; j < loop.numVertices(); ++j) {
      assertTrue(new S1Angle(loop.vertex(j), loop.vertex(j + 1)).greaterThan(minEdgeLength));
    }
  }

  private static S2Polygon
      constructorBenchmarkSingleLoop(int numVertices, GeometryTestCase testUtils) {
    S1Angle radius = S1Angle.degrees(1);
    S2Loop loop = S2Loop.makeRegularLoop(testUtils.randomPoint(), radius, numVertices);
    return new S2Polygon(loop);
  }

  private static S2Polygon constructorBenchmarkLoopGrid(int numLoops, GeometryTestCase testUtils) {
    // More than one S2ShapeIndex cell.
    int verticesPerLoop = 20;
    LoopGridPolygonFactory factory =
        new LoopGridPolygonFactory(numLoops, verticesPerLoop * numLoops);
    return factory.newPolygon(testUtils);
  }

  private static int benchmarkConstructor(int reps, S2Polygon input) {
    List<S2Loop> loops = Lists.newArrayList();
    // Make a copy so that we can release the loops.
    S2Polygon copy = new S2Polygon();
    copy.copy(input);
    copy.release(loops);
    for (int r = reps; r > 0; --r) {
      S2Polygon p = new S2Polygon(loops);
      p.release(loops);
    }
    return loops.size();
  }

  /**
   * Benchmarks S2Polygon construction including its index and *also* including the time required to
   * construct the S2Loops and their indices.
   */
  private static int benchmarkFullConstructor(int reps, S2Polygon input) {
    int numLoops = input.numLoops();
    List<List<S2Point>> loopVertices = Lists.newArrayList();
    for (int i = 0; i < numLoops; ++i) {
      List<S2Point> loop = Lists.newArrayList();
      S2Loop inputLoop = input.loop(i);
      for (int j = 0; j < inputLoop.numVertices(); ++j) {
        loop.add(inputLoop.vertex(j));
      }
      loopVertices.add(loop);
    }

    int totalVertices = 0;
    for (int r = reps; r > 0; --r) {
      List<S2Loop> loops = Lists.newArrayList();
      for (int i = 0; i < numLoops; ++i) {
        loops.add(new S2Loop(loopVertices.get(i)));
      }
      S2Polygon p = new S2Polygon(loops);
      totalVertices += p.getNumVertices();
    }

    return totalVertices;
  }

  /** Benchmarks containsPoint() on polygons generated by the specific factory. */
  private static int benchmarkContainsPoint(int reps, PolygonFactory factory,
      GeometryTestCase testUtils) {
    int numQuerySamples = 100;
    List<S2Point> queries = new ArrayList<S2Point>();
    S2Polygon polygon = new S2Polygon();
    // Bresenham-type algorithm for polygon sampling.
    int delta = 0;
    // Prevent loop from being optimized away.
    int count = 0;
    for (int r = reps; r > 0; --r) {
      delta -= NUM_POLYGON_SAMPLES;
      if (delta < 0) {
        delta += reps;
        polygon = factory.newPolygon(testUtils);
        for (int i = 0; i < numQuerySamples; ++i) {
          queries.add(testUtils.samplePoint(polygon.getRectBound()));
        }
      }
      count += polygon.contains(queries.get(r % numQuerySamples)) ? 1 : 0;
    }

    return count;
  }

  /** Helps benchmark whether a polygon generated by 'factory' contains() itself. */
  private static void benchmarkContainsSelf(int reps, PolygonFactory factory,
      GeometryTestCase testUtils) {
    S2Polygon polygon = factory.newPolygon(testUtils);
    // Make a separate copy to avoid any special checks for A.contains(A).
    S2Polygon copy = new S2Polygon();
    copy.copy(polygon);
    for (int r = reps; r > 0; --r) {
      assertTrue(polygon.contains(copy));
    }
  }

  /** Returns a single loop polygon that contains the given polygon. */
  private static S2Polygon getBoundingPolygon(S2Polygon polygon) {
    S2Cap cap = polygon.getCapBound();
    List<S2Loop> loops = new ArrayList<>(1);
    int numCapEdges = 8;
    loops.add(S2Loop.makeRegularLoop(cap.axis(),
        S1Angle.radians(cap.angle().radians() / Math.cos(Math.PI / numCapEdges)), numCapEdges));
    return new S2Polygon(loops);
  }

  /** 
   * Benchmarks whether a polygon generated by 'factory' intersects() its 'complement' (obtained
   * by subtracting the polygon from a simple bound).
   */
  private static void benchmarkIntersectsCmpl(int reps, PolygonFactory factory,
      GeometryTestCase testUtils) {
    S2Polygon polygon = factory.newPolygon(testUtils);
    S2Polygon bound = getBoundingPolygon(polygon);
    S2Polygon complement = new S2Polygon();
    complement.initToDifference(bound,  polygon);
    for (int r = reps; r > 0; --r) {
      assertFalse(polygon.intersects(complement));
    }
  }

  /**
   * Benchmarks an operation (initToIntersection, etc.) between a polygon generated by 'factory'
   * and every cell in a covering of that polygon.
   */
  private static int benchmarkOpWithCovering(int reps, int maxCells, PolygonFactory factory,
      PolygonOperation operation, GeometryTestCase testUtils) {
    S2RegionCoverer coverer = new S2RegionCoverer();
    coverer.setMaxCells(maxCells);
    S2Polygon polygon = new S2Polygon();
    S2Polygon result = new S2Polygon();
    // Initialize the ArrayList covering to have some size greater than 0.
    ArrayList<S2CellId> covering = new ArrayList<>(16);
    // Bresenham-type algorithm for polygon sampling.
    int delta = 0;
    for (int r = reps, icover = 0; r > 0; --r, ++icover) {
      delta -= NUM_POLYGON_SAMPLES;
      if (delta < 0) {
        delta += reps;
        polygon = factory.newPolygon(testUtils);
        coverer.getCovering(polygon, covering);
      }
      if (icover >= covering.size()) {
        icover = 0;
      }
      S2Cell cell = new S2Cell(covering.get(icover));
      S2Polygon cellPoly = new S2Polygon(cell);
      // Suppress the warning here and elsewhere because we just want to time the operation.
      result = operation.operation(polygon, cellPoly);
    }
    
    return result.numLoops();
  }

  /** Benchmarks an operation between a polygon generated by 'factory' and itself. */
  private static int benchmarkOpWithSelf(int reps, PolygonFactory factory,
      PolygonOperation operation, GeometryTestCase testUtils) {
    S2Polygon polygon = new S2Polygon();
    S2Polygon copy = new S2Polygon();
    S2Polygon result = new S2Polygon();
    // Bresenham-type algorithm for polygon sampling.
    int delta = 0;
    for (int r = reps; r > 0; --r) {
      delta -= NUM_POLYGON_SAMPLES;
      if (delta < 0) {
        delta += reps;
        polygon = factory.newPolygon(testUtils);
        // Make a separate copy to avoid any special optimizations for operations between a polygon
        // and itself.
        copy.copy(polygon);
      }
      result = operation.operation(polygon, copy);
    }
    
    return result.numLoops();
  }

  /**
   * Benchmarks an operation between a polygon generated by 'factory' and a simple bounding
   * polygon.
   */
  private static int benchmarkOpWithBound(int reps, PolygonFactory factory,
      PolygonOperation operation, GeometryTestCase testUtils) {
    S2Polygon polygon = new S2Polygon();
    S2Polygon bound = new S2Polygon();
    S2Polygon result = new S2Polygon();
    // Bresenham-type algorithm for polygon-sampling.
    int delta = 0;
    for (int r = reps; r > 0; --r) {
      delta -= NUM_POLYGON_SAMPLES;
      if (delta < 0) {
        delta += reps;
        polygon = factory.newPolygon(testUtils);
        bound = getBoundingPolygon(polygon);
      }
      result = operation.operation(polygon, bound);
    }
    
    return result.numLoops();
  }
}