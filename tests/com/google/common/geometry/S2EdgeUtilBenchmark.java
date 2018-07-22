package com.google.common.geometry;

import com.google.caliper.Benchmark;

import java.util.ArrayList;
import java.util.List;

/** 
 * A benchmark for {@link S2EdgeUtil}. Note that return values of benchmark functions are
 * ignored by Caliper.
 */
public class S2EdgeUtilBenchmark {  
  GeometryTestCase testUtils;
  private static final int NUM_POINTS = 400;

  void setUp() {
    testUtils = new GeometryTestCase();
    testUtils.setUp();
  }

  /**
   * Allows us to use the S2Edgeutil methods edgeOrVertexCrossing() and robustCrossing() as inputs
   * to the benchmarkCrossing method.
   */
  private interface CrossingFunction {
    boolean crossing(S2Point a, S2Point b, S2Point c, S2Point d);
  }
  
  static class EdgeOrVertexCrossing implements CrossingFunction {
    @Override
    public boolean crossing(S2Point a, S2Point b, S2Point c, S2Point d) {
      return S2EdgeUtil.edgeOrVertexCrossing(a, b, c, d);
    }
  }
  
  static class RobustCrossing implements CrossingFunction {
    @Override
    public boolean crossing(S2Point a, S2Point b, S2Point c, S2Point d) {
      return S2EdgeUtil.robustCrossing(a, b, c, d) > 0;
    }
  }
  
  /** Benchmarks the given crossing function on a list of randomly generated points. */
  private static int benchmarkCrossing(CrossingFunction crossing, int reps,
      GeometryTestCase testUtils) {
    // We want to avoid cache effects, so numPoints should be small enough so that the points can
    // be in L1 cache.  The size of an S2Point is at most 40 bytes, so 400 will only take at most
    // ~16 KiB of 64 KiB of L1 cache.
    List<S2Point> p = new ArrayList<>(NUM_POINTS);
    for (int i = 0; i < NUM_POINTS; ++i) {
      p.add(testUtils.randomPoint());
    }
    
    int numCrossings = 0;
    for (int r = reps; r > 0; --r) {
      S2Point a = p.get((r + 0) % NUM_POINTS);
      S2Point b = p.get((r + 1) % NUM_POINTS);
      S2Point c = p.get((r + 2) % NUM_POINTS);
      S2Point d = p.get((r + 3) % NUM_POINTS);
    
      if (crossing.crossing(a, b, c, d)) {
        ++numCrossings;
      }
    }
    
    return numCrossings;
  }
  
  @Benchmark void edgeOrVertexCrossing(int reps) {
    setUp();
    benchmarkCrossing(new EdgeOrVertexCrossing(), reps, testUtils);
  }
  
  @Benchmark void robustCrossing(int reps) {
    setUp();
    benchmarkCrossing(new RobustCrossing(), reps, testUtils);
  }
  
  /** Benchmarks robustCrossing() when approximately 1/4 of of the edges cross. */
  @Benchmark int robustCrosserEdgesCross(int reps) {
    setUp();
    // Copied from benchmarkCrossing() above.
    int numPoints = 400;
    List<S2Point> p = new ArrayList<>(numPoints);
    for (int i = 0; i < numPoints; ++i) {
      p.add(testUtils.randomPoint());
    }
    // Approximately 1/4th of points will cross the edge 'ab'.
    S2Point a = testUtils.randomPoint();
    S2Point b = S2Point.normalize(S2Point.add(S2Point.neg(a), new S2Point(0.1, 0.1, 0.1)));
    
    S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(a, b, p.get(0));
    int crossingSum = 0;
    for (int r = reps; r > 0; --r) {
      S2Point d = p.get(r % numPoints);
      
      if (crosser.robustCrossing(d) > 0) {
        crossingSum += crosser.robustCrossing(d);
      }
    }
    
    return crossingSum;
  }
}