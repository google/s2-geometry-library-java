package com.google.common.geometry;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.caliper.BeforeExperiment;
import com.google.caliper.Benchmark;
import com.google.caliper.Param;
import com.google.caliper.api.BeforeRep;
import com.google.caliper.api.Macrobenchmark;
import com.google.common.base.Preconditions;

import java.util.ArrayList;
import java.util.List;

/** 
 * A benchmark for {@link S2Loop}. Caliper ignores returned values, but these are computed and
 * returned so that corresponding for loops don't get optimized away. We create our own instances
 * of GeometryTestCase because Caliper currently requires that benchmarking classes cannot extend
 * anything other than Object.
 */
public class S2LoopBenchmark {
  /**
   * Number of random loops or loop pairs (in order to average over various loop positions on the
   * sphere).
   */
  private static final int NUM_LOOP_SAMPLES = 16;

  /** Default number of vertices for loops.*/
  private static final int DEFAULT_NUM_VERTICES = 4096;

  /** Desired distance between nested loops, expressed as a multiple of the edge length. */
  private static final double DEFAULT_NESTED_GAP_MULTIPLE = 5.0;

  static class LoopConstructor {
    @Param({"4", "8", "64", "512", "4096", "32768", "262144"})
    int numVertices;

    /** Tests the speed of constructing loops from a given list of vertices. */
    @Benchmark int loopConstructor(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      List<List<S2Point>> loops = new ArrayList<>(NUM_LOOP_SAMPLES);
      for (int i = 0; i < NUM_LOOP_SAMPLES; ++i) {
        S2Loop loop = S2Loop.makeRegularLoop(
            testUtils.randomPoint(), GeometryTestCase.DEFAULT_RADIUS, numVertices);
        List<S2Point> loopVertices = new ArrayList<>(numVertices);
        for (int v = 0; v < numVertices; ++v) {
          loopVertices.add(loop.vertex(v));
        }
        loops.add(loopVertices);
      }
      
      // Cycle through the list of sample loops, constructing them one at a time.
      int vertexSum = 0;
      int index = 0;
      for (int r = reps; r > 0; --r) {
        vertexSum += new S2Loop(loops.get(index)).numVertices();
        if (++index == NUM_LOOP_SAMPLES) {
          index = 0;
        }
      }

      return vertexSum;
    }
  }

  /**
   * Returns a list of 'NUM_LOOP_SAMPLES' regular loops, each of 'numVertices' vertices, with
   * each loop positioned randomly (but repeatably) on the sphere.
   */
  protected static List<S2Loop> getRegularLoops(S1Angle radius, int numVertices,
      GeometryTestCase testUtils) {
    List<S2Loop> loops = new ArrayList<>(NUM_LOOP_SAMPLES);
    for (int i = 0; i < NUM_LOOP_SAMPLES; ++i) {
      loops.add(S2Loop.makeRegularLoop(testUtils.randomPoint(), radius, numVertices));
    }

    return loops;
  }

  static class IsValid  {
    /** Number of vertices in each loop. */
    @Param({"4", "16", "64", "256", "2048", "32768"})
    int numVertices;
    
    private List<S2Loop> loops;

    @BeforeExperiment void setUpExperiment() {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      loops = index(getRegularLoops(GeometryTestCase.DEFAULT_RADIUS, numVertices, testUtils));
    }
    
    /** Tests the speed of isValid() on regular loops. */
    @Benchmark int isValid(int reps) {
      int index = 0;
      int numValid = 0;
      for (int r = reps; r > 0; --r) {
        if (loops.get(index).isValid()) {
          numValid++;
        }
        if (++index == loops.size()) {
          index = 0;
        }
      }

      return numValid;
    }
  }

  static class LoopContains {
    /** Number of vertices in each loop.*/
    @Param({"4", "8", "16", "32", "64", "128", "256", "512", "1024", "2048", "8192", "32768",
        "262144"})
    int numVertices;

    private static final int NUM_QUERIES_PER_LOOP = 100;
    private List<S2Loop> loops;
    private List<List<S2Point>> allQueries;
    
    @BeforeExperiment void setUpExperiment() {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      // Create a list of sample loops.
      loops = index(getRegularLoops(GeometryTestCase.DEFAULT_RADIUS, numVertices, testUtils));

      // For each loop created, create 'NUM_QUERIES_PER_LOOP' query points.  The query points for
      // the i-th loop are in allQueries.get(i).
      allQueries = new ArrayList<>(loops.size()); 
      for (int i = 0; i < loops.size(); i++) {
        S2LatLngRect bound = loops.get(i).getRectBound();

        List<S2Point> loopQueries = new ArrayList<>(NUM_QUERIES_PER_LOOP);
        for (int j = 0; j < NUM_QUERIES_PER_LOOP; ++j) {
          S2Point p = testUtils.samplePoint(bound);
          assertTrue(bound.contains(p));
          loopQueries.add(p);
        }
        allQueries.add(loopQueries);
      }
    }

    /** 
     * Tests the speed of checking whether a regular loop contains a point.  The query points for a
     * loop are chosen so that they all lie in the loop's bounding rectangle (to avoid the 
     * quick-rejection code path.
     */
    @Benchmark int containsPoint (int reps) {
      // Cycle through the loops and queries such that contains() is called 'reps' times.
      int numContainedPoints = 0;
      for (int iloop = 0, iquery = 0, r = reps; r > 0; --r) {
        S2Point query = allQueries.get(iloop).get(iquery);
        if (loops.get(iloop).contains(query)) {
          numContainedPoints++;
        }
        if (++iloop == loops.size()) {
          iloop = 0;
          if (++iquery == NUM_QUERIES_PER_LOOP) {
            iquery = 0;
          }
        }
      }

      return numContainedPoints;
    }
  }

  static class ConstructAndContains {
    /** Number of times to call contains(S2Point) during one iteration of our benchmark tests. */
    @Param({"1", "8", "64", "512", "4096", "32768"})
    int numCalls;

    /** Whether the loop should contain the query point. */
    @Param({"true", "false"})
    boolean shouldContain;

    /** Number of vertices per loop. */
    private static final int NUM_VERTICES_PER_LOOP = 1024;
    S2Loop loop;
    S2Point query;    
    GeometryTestCase testUtils;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /** Constructs a loop of 'numVertices' vertices and calls contains(S2Point) 'numCalls' times. */
    @Macrobenchmark int constructAndContains() {
      query = testUtils.randomPoint();
      // Build a regular loop either around the query point or its antipode, based on whether the
      // loop should contain the query point.
      loop = S2Loop.makeRegularLoop(shouldContain ? query : S2Point.neg(query),
          GeometryTestCase.DEFAULT_RADIUS, NUM_VERTICES_PER_LOOP);
      
      int sumContains = 0;
      for (int n = numCalls; n > 0; --n) {
        sumContains += loop.contains(query) ? 1 : 0;
      }
      
      return sumContains;
    }
  }
  
  /** 
   * The loop relation methods (contains, intersects, compareBoundary) scale similarly with the
   * number of vertices.  So to save time, we only benchmark compareBoundary() on a full range of
   * vertex counts.  For contains() and intersects() we just run a single test using 
   * DEFAULT_NUM_VERTICES.
   */
  static class CompareBoundaryContains {
    @Param({"8", "64", "512", "4096", "32768"})
    int numVertices;
    
    /** Benchmarks compareBoundary() when one loop contains the other. */
    @Benchmark void compareBoundaryContains(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      List<S2Loop> loops = index(GeometryTestCase.makeNestedLoopPair(
          GeometryTestCase.DEFAULT_RADIUS, DEFAULT_NESTED_GAP_MULTIPLE, numVertices, 
          testUtils.randomPoint()));
      for (int r = reps; r > 0; --r) {
        assertEquals(1, loops.get(0).compareBoundary(loops.get(1)));
      }
    }
  }
  
  static class CompareBoundaryCrosses {
    @Param({"8", "64", "512", "4096", "32768"})
    int numVertices;

    /** Benchmarks compareBoundary() when the loops cross each other. */
    @Benchmark void compareBoundaryCrosses(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      List<S2Loop> loops = index(testUtils.makeCrossingLoopPairDefault(
          numVertices, testUtils.randomPoint(), testUtils.randomPoint()));
      for (int r = reps; r > 0; --r) {
        assertEquals(0, loops.get(0).compareBoundary(loops.get(1)));
      }
    }
  }
  
  static class CompareBoundaryDisjoint {
    @Param({"8", "64", "512", "4096", "32768"})
    int numVertices;
    
    /** Benchmarks compareBoundary() when the loops are disjoint. */
    @Benchmark void compareBoundaryDisjoint(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      List<S2Loop> loops = index(GeometryTestCase.makeDisjointLoopPair(
          GeometryTestCase.DEFAULT_RADIUS, DEFAULT_NESTED_GAP_MULTIPLE, numVertices,
          testUtils.randomPoint()));
      for (int r = reps; r > 0; --r) {
        assertEquals(-1, loops.get(0).compareBoundary(loops.get(1)));
      }
    }
  }

  static class ContainsContains {
    /** Benchmarks contains() when one loop contains the other. */
    @Benchmark void containsContains(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      List<S2Loop> loops = index(GeometryTestCase.makeNestedLoopPair(
          GeometryTestCase.DEFAULT_RADIUS, DEFAULT_NESTED_GAP_MULTIPLE, DEFAULT_NUM_VERTICES,
          testUtils.randomPoint()));
      for (int r = reps; r > 0; --r) {
        assertTrue(loops.get(0).contains(loops.get(1)));
      }
    }
  }

  static class ContainsCrosses {
    /** Benchmarks contains() when the loops cross each other. */
    @Benchmark void containsCrosses(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      List<S2Loop> loops = index(testUtils.makeCrossingLoopPairDefault(DEFAULT_NUM_VERTICES,
          testUtils.randomPoint(),  testUtils.randomPoint()));
      for (int r = reps; r > 0; --r) {
        assertFalse(loops.get(0).contains(loops.get(1)));
      }
    }
  }
  
  static class ContainsDisjoint {
    /** Benchmarks contains() when the loops are disjoint. */
    @Benchmark void containsDisjoint(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      List<S2Loop> loops = index(GeometryTestCase.makeDisjointLoopPair(
          GeometryTestCase.DEFAULT_RADIUS, DEFAULT_NESTED_GAP_MULTIPLE, DEFAULT_NUM_VERTICES,
          testUtils.randomPoint()));
      for (int r = reps; r > 0; --r) {
        assertFalse(loops.get(0).contains(loops.get(1)));
      }
    }
  }

  static class IntersectsContains {
    /** Benchmarks intersects() when one loop contains the other. */
    @Benchmark void intersectsContains(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      List<S2Loop> loops = index(GeometryTestCase.makeNestedLoopPair(
          GeometryTestCase.DEFAULT_RADIUS, DEFAULT_NESTED_GAP_MULTIPLE, DEFAULT_NUM_VERTICES,
          testUtils.randomPoint()));
      for (int r = reps; r > 0; --r) {
        assertTrue(loops.get(0).intersects(loops.get(1)));
      }
    }
  }

  /** Benchmarks intersects() when the two loops cross. */
  static class IntersectsCrosses {
    /** Benchmarks intersects() when the loops cross each other. */
    @Benchmark void intersectsCrosses(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      List<S2Loop> loops = index(testUtils.makeCrossingLoopPairDefault(DEFAULT_NUM_VERTICES,
          testUtils.randomPoint(), testUtils.randomPoint()));
      for (int r = reps; r > 0; --r) {
        assertTrue(loops.get(0).intersects(loops.get(1)));
      }
    }
  }
  
  static class IntersectsDisjoint {    
    /** Benchmarks intersects() when the loops are disjoint. */
    @Benchmark void intersectsDisjoint(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      List<S2Loop> loops = index(GeometryTestCase.makeDisjointLoopPair(
          GeometryTestCase.DEFAULT_RADIUS, DEFAULT_NESTED_GAP_MULTIPLE, DEFAULT_NUM_VERTICES,
          testUtils.randomPoint()));
      for (int r = reps; r > 0; --r) {
        assertFalse(loops.get(0).intersects(loops.get(1)));
      }
    }
  }

  static class ContainsVsRadiusMeters {
    /** Radius of the outer loop. */
    @Param({"1", "8", "64", "1024", "65536"})
    int meters;

    private GeometryTestCase testUtils;
    private List<S2Loop> loops;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    @BeforeRep void setUpRep() {
      loops = index(GeometryTestCase.makeNestedLoopPair(
          GeometryTestCase.metersToAngle(meters), DEFAULT_NESTED_GAP_MULTIPLE, DEFAULT_NUM_VERTICES,
          testUtils.randomPoint()));
    }
    
    /**
     * Benchmarks contains() on nested loops as a function of the loop radius.  Performance
     * degrades on very small loops because S2.expensiveCCW() is invoked more often.
     */
    @Macrobenchmark void containsVsRadiusMeters() {
      assertTrue(loops.get(0).contains(loops.get(1)));
    }
  }

  static class ContainsVsGapEdgeMultiple {
    /** See comments for getNestedLoopPair(). */
    @Param({"0", "1", "8", "64"})
    int gapEdgeMultiple;

    private GeometryTestCase testUtils;
    private List<S2Loop> loops;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    @BeforeRep void setUpRep() {
      loops = index(GeometryTestCase.makeNestedLoopPair(GeometryTestCase.DEFAULT_RADIUS,
          gapEdgeMultiple, DEFAULT_NUM_VERTICES, testUtils.randomPoint()));
    }
    
    /**
     * Benchmarks contains() on nested loops as a function of the distance between the loops
     * (expressed as a multiple of the loop edge length).  Performance degrades as the loops get
     * very close together because spatial indexing is not as effective at pruning the
     * intersection candidates.
     */
    @Macrobenchmark void containsVsGapEdgeMultiple () {
      assertTrue(loops.get(0).contains(loops.get(1)));
    }
  }

  static class IntersectsCrossesVsLogRadiusRatio {
    /** 
     * The C++ version also has "10" as a parameter, but it is omitted here because intersects()
     * will return false for evaluating the loops.  The current value of DEFAULT_RADIUS_KM is 10 km.
     * A loop 10 * 2^10 km in radius (along the surface of the sphere) will cover just a bit more
     * than a hemisphere.  So this is indeed a big loop, and perhaps a special case.
     * TODO(eengle) will investigate and enable logRadiusRatio = 10.
     */
    @Param({"-10", "-5", "0", "5"})
    int logRadiusRatio;
    
    /** 
     * Benchmarks intersects() on crossing loops as a function of the relative sizes of the two
     * loops (e.g., one loop radius much larger than the other).  Performance of spatial indexing
     * can degrade when one loop is much larger than the other.
     */
    @Benchmark void intersectsCrossesVsLogRatio(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      S1Angle aRadius = GeometryTestCase.DEFAULT_RADIUS;
      S1Angle bRadius = S1Angle.radians(Math.scalb(aRadius.radians(), logRadiusRatio));
      List<S2Loop> loops = index(testUtils.makeCrossingLoopPair(aRadius, bRadius,
          DEFAULT_NUM_VERTICES, testUtils.randomPoint(), testUtils.randomPoint()));
      for (int r = reps; r > 0; --r) {
        assertTrue(loops.get(0).intersects(loops.get(1)));
      }
    }
  }
  
  /** Ensures the index is built for the given list of loops and returns them. */
  private static List<S2Loop> index(List<S2Loop> loops) {
    for (S2Loop loop : loops) {
      loop.contains(S2Cell.fromFace(0));
      Preconditions.checkState(loop.index.isFresh());
    }
    return loops;
  }
}
