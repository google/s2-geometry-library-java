package com.google.common.geometry;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.caliper.BeforeExperiment;
import com.google.caliper.Benchmark;
import com.google.caliper.Param;
import com.google.caliper.api.BeforeRep;
import com.google.caliper.api.Macrobenchmark;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;

import java.util.ArrayList;
import java.util.List;

/**
 * A benchmark for {@link S2Polygon}. Caliper ignores returned values, but these are computed and
 * returned so that corresponding for loops don't get optimized away.
 */
public class S2PolygonBenchmark {
  /** Default radius for loops. */
  private static final double DEFAULT_RADIUS_KM = 10.0;

  /** Default dimension for fractal loops. */
  private static final double LOOP_FRACTAL_DIMENSION = 1.2;

  /** Spacing between adjacent loops (e.g., in a grid), as a fraction of the loop radius. */
  private static final double LOOP_SEPARATION_FRACTION = 0.1;

  static class ConstructorSingleLoop {
    @Param({"8", "1024", "1048576"})
    int numVertices;

    private GeometryTestCase testUtils;

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

    private GeometryTestCase testUtils;

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
      S2Polygon polygon = index(
          GeometryTestCase.concentricLoopsPolygon(testUtils.randomPoint(), numLoops, numVertices));
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
  
  static class Covering {
    @Param({"8", "64", "512"})
    int maxCells;

    @Param({"8", "128", "2048", "32768", "524288"})
    int numVertices;

    private GeometryTestCase testUtils; 
    private NestedLoopsPolygonFactory nestedLoopsPolygonFactory;
    private NestedFractalsPolygonFactory nestedFractalsPolygonFactory;
    private S2Polygon nestedLoopsPolygon;
    private S2Polygon nestedFractalsPolygon;
    private S2RegionCoverer coverer;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
      nestedLoopsPolygonFactory = new NestedLoopsPolygonFactory(2, numVertices);
      nestedFractalsPolygonFactory = new NestedFractalsPolygonFactory(1, numVertices);
      coverer = new S2RegionCoverer();
      coverer.setMaxCells(maxCells);
    }
    
    /** We create both polygons, although only one is needed for each macrobenchmark. */
    @BeforeRep void setUpRep() {
      nestedLoopsPolygon = index(nestedLoopsPolygonFactory.newPolygon(testUtils));
      nestedFractalsPolygon = index(nestedFractalsPolygonFactory.newPolygon(testUtils));
    }
    
    @Macrobenchmark void annulusCovering() {
      ArrayList<S2CellId> covering = new ArrayList<>();
      coverer.getCovering(nestedLoopsPolygon, covering);
    }
    
    @Macrobenchmark void fractalCovering() {
      ArrayList<S2CellId> covering = new ArrayList<>();
      coverer.getCovering(nestedFractalsPolygon, covering);
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
    
    private GeometryTestCase testUtils;
    private LoopGridPolygonFactory factory;
    private S2Polygon polygon;
    private S2Polygon complement;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
      factory = new LoopGridPolygonFactory(numLoops, totalNumVertices);
    }
    
    @BeforeRep void setUpRep() {
      polygon = index(factory.newPolygon(testUtils));
      S2Polygon bound = getBoundingPolygon(polygon);
      complement = new S2Polygon();
      complement.initToDifference(bound,  polygon);
      index(complement);
    }
    
    /**
     * Benchmarks whether a loop grid polygon intersects() its complement (obtained by subtracting
     * the polygon from a simple bound).
     */
    @Macrobenchmark void intersectsCmplLoopGrid() {
      assertFalse(polygon.intersects(complement));
    }
  }

  static class IntersectsCmplNestedFractals {
    @Param({"2"})
    int numLoops;

    @Param({"512", "262144"})
    int totalNumVertices;

    private GeometryTestCase testUtils;
    private NestedFractalsPolygonFactory factory;
    private S2Polygon polygon;
    private S2Polygon complement;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
      factory = new NestedFractalsPolygonFactory(numLoops, totalNumVertices);
    }

    @BeforeRep void setUpRep() {
      polygon = index(factory.newPolygon(testUtils));
      S2Polygon bound = getBoundingPolygon(polygon);
      complement = new S2Polygon();
      complement.initToDifference(bound,  polygon);
      index(complement);
    }
    
    /** Similar to containsSelfNestedFractals, so we just do some spot checks. */
    @Macrobenchmark void intersectsCmplNestedFractals() {
      assertFalse(polygon.intersects(complement));
    }
  }

  /** 
   * The only interesting case is intersection, since for union and difference the output polygon
   * includes almost all of the input polygon.
   */
  static class IntersectFractalWithCovering {
    @Param({"1", "2", "32", "128", "4096"})
    int maxCells;

    @Param({"128", "32768"})
    int totalNumVertices;

    private GeometryTestCase testUtils;
    private NestedFractalsPolygonFactory factory;
    private S2RegionCoverer coverer;
    private S2Polygon polygon, cellPoly, result;
    private ArrayList<S2CellId> covering;
    /** The index of the cell of the current covering. */
    private int icover;
    /** True if we should get a new polygon before the next rep. */
    private boolean getNewPolygon;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
      factory = new NestedFractalsPolygonFactory(1, totalNumVertices);
      getNewPolygon = true;
      icover = 0;
      result = new S2Polygon();
    }

    @BeforeRep void setUpRep() {
      if (getNewPolygon) {
        icover = 0;
        polygon = index(factory.newPolygon(testUtils));
        covering = new ArrayList<>(16);
        coverer = new S2RegionCoverer();
        coverer.setMaxCells(maxCells);
        coverer.getCovering(polygon, covering);
      }
      S2Cell cell = new S2Cell(covering.get(icover++));
      cellPoly = index(new S2Polygon(cell));
      getNewPolygon = (icover == covering.size());
      result = new S2Polygon();
    }
    
    @Macrobenchmark void intersectFractalWithCovering() {
      result.initToIntersection(polygon, cellPoly);
    }
  }

  static class IntersectLoopGridWithCovering {
    @Param({"1", "2", "32", "128", "4096"})
    int maxCells;

    @Param({"128", "32768"})
    int totalNumVertices;

    private static final int NUM_VERTICES_PER_LOOP = 4;
    private GeometryTestCase testUtils;
    private LoopGridPolygonFactory factory;
    private S2RegionCoverer coverer;
    private S2Polygon polygon, cellPoly, result;
    private ArrayList<S2CellId> covering;
    /** The index of the cell of the current covering. */
    private int icover;
    /** True if we should get a new polygon before the next rep. */
    private boolean getNewPolygon;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
      int numLoops = totalNumVertices / NUM_VERTICES_PER_LOOP;
      factory = new LoopGridPolygonFactory(numLoops, totalNumVertices);
      getNewPolygon = true;
      icover = 0;
      result = new S2Polygon();
    }

    @BeforeRep void setUpRep() {
      if (getNewPolygon) {
        icover = 0;
        polygon = index(factory.newPolygon(testUtils));
        covering = new ArrayList<>(16);
        coverer = new S2RegionCoverer();
        coverer.setMaxCells(maxCells);
        coverer.getCovering(polygon, covering);
      }
      S2Cell cell = new S2Cell(covering.get(icover++));
      cellPoly = index(new S2Polygon(cell));
      getNewPolygon = (icover == covering.size());
      result = new S2Polygon();
    }
    
    @Macrobenchmark void intersectLoopGridWithCovering() {
      result.initToIntersection(polygon, cellPoly);
    }
  }

  static class UnionNestedFractalWithSelf {
    @Param({"1", "16"})
    int numLoops;

    @Param({"32768"})
    int totalNumVertices;
    
    private GeometryTestCase testUtils;
    private NestedFractalsPolygonFactory factory;
    private S2Polygon polygon;
    private S2Polygon copy;
    private S2Polygon result;
    
    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
      factory = new NestedFractalsPolygonFactory(numLoops, totalNumVertices);  
    }
    
    @BeforeRep void setUpRep() {
      polygon = index(factory.newPolygon(testUtils));
      copy = new S2Polygon();
      copy.copy(polygon);
      index(copy);
      result = new S2Polygon();
    }
    
    @Macrobenchmark void unionNestedFractalWithSelf() {   
      result.initToUnion(polygon, copy);
    }

    @Macrobenchmark void subtractNestedFractalFromSelf() {   
      result.initToDifference(polygon, copy);
    }
  }
  
  static class UnionLoopGridWithBound {
    @Param({"1", "8", "64", "512", "4096"})
    int numLoops;

    @Param({"32768"})
    int totalNumVertices;

    private GeometryTestCase testUtils;
    private LoopGridPolygonFactory factory;
    private S2Polygon polygon;
    private S2Polygon bound;
    private S2Polygon result;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
      factory = new LoopGridPolygonFactory(numLoops, totalNumVertices);
    }

    @BeforeRep void setUpRep() {
      polygon = index(factory.newPolygon(testUtils));
      bound = index(getBoundingPolygon(polygon));
      result = new S2Polygon();
    }

    @Macrobenchmark void unionLoopGridWithBound() {
      result.initToUnion(polygon, bound);
    }

    @Macrobenchmark void subtractBoundFromLoopGrid() {
      result.initToDifference(polygon, bound);
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
    int numQuerySamples = 10;
    int numPolygons = 5;
    List<S2Polygon> polygons = new ArrayList<S2Polygon>(numPolygons);
    List<List<S2Point>> queries = new ArrayList<List<S2Point>>(numPolygons);
    for (int i = 0; i < numPolygons; ++i) {
      polygons.add(index(factory.newPolygon(testUtils)));
    }
    for (int i = 0; i < numPolygons; ++i) {
      List<S2Point> polygonSpecificQueries = new ArrayList<S2Point>(numQuerySamples);
      S2Polygon polygon = polygons.get(i);
      for (int j = 0; j < numQuerySamples; ++j) {
        polygonSpecificQueries.add(testUtils.samplePoint(polygon.getRectBound()));
      }
      queries.add(polygonSpecificQueries);
    }
    // Prevent loop from being optimized away.
    int count = 0;
    int i = 0;
    S2Polygon polygon = polygons.get(i);
    for (int r = reps; r > 0; --r) {
      count += polygon.contains(queries.get(i).get(r % numQuerySamples)) ? 1 : 0;
      if (r % numQuerySamples == 0) {
        i++;
        if (i == numPolygons) {
          i = 0;
        }
        polygon = polygons.get(i);
      }
    }

    return count;
  }

  /** Helps benchmark whether a polygon generated by 'factory' contains() itself. */
  private static void benchmarkContainsSelf(int reps, PolygonFactory factory,
      GeometryTestCase testUtils) {
    S2Polygon polygon = index(factory.newPolygon(testUtils));
    // Make a separate copy to avoid any special checks for A.contains(A).
    S2Polygon copy = new S2Polygon();
    copy.copy(polygon);
    index(copy);
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
   * Builds the S2ShapeIndex for this polygon, so that the construction doesn't interfere with
   * timings.
   */
  private static S2Polygon index(S2Polygon polygon) {
    // Calling contains(S2Cell) will always call polygon.index.iterator(), which builds the index.
    polygon.contains(S2Cell.fromFace(0));
    Preconditions.checkState(polygon.index.isFresh());
    return polygon;
  }
}
