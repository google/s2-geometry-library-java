package com.google.common.geometry;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.caliper.BeforeExperiment;
import com.google.caliper.Benchmark;
import com.google.caliper.Param;
import com.google.caliper.api.BeforeRep;
import com.google.caliper.api.Macrobenchmark;
import com.google.common.collect.Lists;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

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

  /** The default radius for loops. */
  private static final double DEFAULT_RADIUS_KM = 10.0;

  /** Default number of vertices for loops.*/
  private static final int DEFAULT_NUM_VERTICES = 4096;

  /** Desired distance between nested loops, expressed as a multiple of the edge length. */
  private static final double DEFAULT_NESTED_GAP_MULTIPLE = 5.0;

  /** Desired ratio between the radii of crossing loops. */
  private static final double DEFAULT_CROSSING_RADIUS_RATIO = 10.0;

  static class LoopConstructor {
    @Param({"4", "8", "64", "512", "4096", "32768", "262144"})
    int numVertices;

    /** Tests the speed of constructing loops from a given list of vertices. */
    @Benchmark int loopConstructor(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      // Create a list of sample loops.
      List<List<S2Point>> loops = new ArrayList<>(NUM_LOOP_SAMPLES);
      for (int i = 0; i < NUM_LOOP_SAMPLES; ++i) {
        S2Loop loop = S2Loop.makeRegularLoop(testUtils.randomPoint(),
            GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM), numVertices);
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

    /** Tests the speed of isValid() on regular loops. */
    @Benchmark int isValid(int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      List<S2Loop> loops = getRegularLoops(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM),
          numVertices, testUtils);

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

    /** 
     * Tests the speed of checking whether a regular loop contains a point.  The query points for a
     * loop are chosen so that they all lie in the loop's bounding rectangle (to avoid the 
     * quick-rejection code path.
     */
    @Benchmark int containsPoint (int reps) {
      GeometryTestCase testUtils = new GeometryTestCase();
      testUtils.setUp();
      // Create a list of sample loops.
      List<S2Loop> loops = getRegularLoops(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM),
          numVertices, testUtils);

      // For each loop created, create 'NUM_QUERIES_PER_LOOP' query points.  The query points for
      // the i-th loop are in allQueries.get(i).
      List<List<S2Point>> allQueries = new ArrayList<>(loops.size()); 
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
      loop = S2Loop.makeRegularLoop(
          shouldContain ? query : S2Point.neg(query),
          GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM), NUM_VERTICES_PER_LOOP);
      List<S2Point> loopVertices = new ArrayList<>(NUM_VERTICES_PER_LOOP);
      for (int v = 0; v < NUM_VERTICES_PER_LOOP; ++v) {
        loopVertices.add(loop.vertex(v));
      }
      loop = new S2Loop(loopVertices);
      
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

    private GeometryTestCase testUtils;
    private List<S2Loop> loops;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }
    
    @BeforeRep void setUpRep() {
      loops = getNestedLoopPair(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM),
          DEFAULT_NESTED_GAP_MULTIPLE, numVertices, testUtils.randomPoint());
    }

    /** Benchmarks compareBoundary() when one loop contains the other. */
    @Macrobenchmark void compareBoundaryContains() {
      assertEquals(1, loops.get(0).compareBoundary(loops.get(1)));
    }
  }
  
  static class CompareBoundaryCrosses {
    @Param({"8", "64", "512", "4096", "32768"})
    int numVertices;

    private GeometryTestCase testUtils;
    private List<S2Loop> loops;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }
    
    @BeforeRep void setUpRep() {
      loops = getCrossingLoopPairDefault(numVertices, testUtils.randomPoint(),
          testUtils.randomPoint());
    }
    
    /** Benchmarks compareBoundary() when the loops cross each other. */
    @Macrobenchmark void compareBoundaryCrosses() {
      assertEquals(0, loops.get(0).compareBoundary(loops.get(1)));
    }
  }
  
  static class CompareBoundaryDisjoint {
    @Param({"8", "64", "512", "4096", "32768"})
    int numVertices;

    private GeometryTestCase testUtils;
    private List<S2Loop> loops;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }
    
    @BeforeRep void setUpRep() {
      loops = getDisjointLoopPair(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM),
          DEFAULT_NESTED_GAP_MULTIPLE, numVertices, testUtils.randomPoint());
     }
    
    /** Benchmarks compareBoundary() when the loops are disjoint. */
    @Macrobenchmark void compareBoundaryDisjoint() {
      assertEquals(1, loops.get(0).compareBoundary(loops.get(1)));
    }
  }

  static class ContainsContains {
    private GeometryTestCase testUtils;
    private List<S2Loop> loops;
    
    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }
    
    @BeforeRep void setUpRep() {
      loops = getNestedLoopPair(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM),
          DEFAULT_NESTED_GAP_MULTIPLE, DEFAULT_NUM_VERTICES, testUtils.randomPoint());
    }
    
    /** Benchmarks contains() when one loop contains the other. */
    @Macrobenchmark void containsContains() {
      assertTrue(loops.get(0).contains(loops.get(1)));
    }
  }

  static class ContainsCrosses {
    private GeometryTestCase testUtils;
    private List<S2Loop> loops;
    
    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }
    
    @BeforeRep void setUpRep() {
      loops = getCrossingLoopPairDefault(DEFAULT_NUM_VERTICES, testUtils.randomPoint(), 
          testUtils.randomPoint());
    }
      
    /** Benchmarks contains() when the loops cross each other. */
    @Macrobenchmark void containsCrosses() {
      assertFalse(loops.get(0).contains(loops.get(1)));
    }
  }
  
  static class ContainsDisjoint {
    private GeometryTestCase testUtils;
    private List<S2Loop> loops;
    
    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }
    
    @BeforeRep void setUpRep() {
      loops = getDisjointLoopPair(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM),
          DEFAULT_NESTED_GAP_MULTIPLE, DEFAULT_NUM_VERTICES, testUtils.randomPoint());
    }
    
    /** Benchmarks contains() when the loops are disjoint. */
    @Macrobenchmark void containsDisjoint() {
      assertTrue(loops.get(0).contains(loops.get(1)));
    }
  }

  static class IntersectsContains {
    private GeometryTestCase testUtils;
    private List<S2Loop> loops;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    @BeforeRep void setUpRep() {
      loops = getNestedLoopPair(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM),
          DEFAULT_NESTED_GAP_MULTIPLE, DEFAULT_NUM_VERTICES, testUtils.randomPoint());
    }

    /** Benchmarks intersects() when one loop contains the other. */
    @Macrobenchmark void intersectsContains() {
      assertTrue(loops.get(0).intersects(loops.get(1)));
    }
  }

  /** Benchmarks intersects() when the two loops cross. */
  static class IntersectsCrosses {
    private GeometryTestCase testUtils;
    private List<S2Loop> loops;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    @BeforeRep void setUpRep() {
      loops = getCrossingLoopPairDefault(DEFAULT_NUM_VERTICES, testUtils.randomPoint(), 
          testUtils.randomPoint());
    }

    /** Benchmarks intersects() when the loops cross each other. */
    @Macrobenchmark void intersectsCrosses() {
      assertTrue(loops.get(0).intersects(loops.get(1)));
    }
  }
  
  static class IntersectsDisjoint {
    private GeometryTestCase testUtils;
    private List<S2Loop> loops;
    
    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }
    
    @BeforeRep void setUpRep() {
      loops = getDisjointLoopPair(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM),
          DEFAULT_NESTED_GAP_MULTIPLE, DEFAULT_NUM_VERTICES, testUtils.randomPoint());
    }
    
    /** Benchmarks intersects() when the loops are disjoint. */
    @Macrobenchmark void intersectsDisjoint() {
      assertTrue(loops.get(0).intersects(loops.get(1)));
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
      loops = getNestedLoopPair(GeometryTestCase.metersToAngle(meters), DEFAULT_NESTED_GAP_MULTIPLE,
          DEFAULT_NUM_VERTICES, testUtils.randomPoint());
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
      loops = getNestedLoopPair(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM), gapEdgeMultiple,
          DEFAULT_NUM_VERTICES, testUtils.randomPoint());
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

    private GeometryTestCase testUtils;
    private List<S2Loop> loops;
    private S1Angle aRadius, bRadius;

    @BeforeExperiment void setUpExperiment() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
      aRadius = GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM);
      bRadius = S1Angle.radians(Math.scalb(aRadius.radians(), logRadiusRatio));
    }

    @BeforeRep void setUpRep() {
      loops = getCrossingLoopPair(aRadius, bRadius, DEFAULT_NUM_VERTICES, testUtils.randomPoint(),
          testUtils.randomPoint());
    }
    
    /** 
     * Benchmarks intersects() on crossing loops as a function of the relative sizes of the two
     * loops (e.g., one loop radius much larger than the other).  Performance of spatial indexing
     * can degrade when one loop is much larger than the other.
     */
    @Macrobenchmark void intersectsCrossesVsLogRatio() {
      assertTrue(loops.get(0).intersects(loops.get(1)));
    }
  }

  /**
   * Returns a pair of nested loops, such that the second loop in the returned list is nested within
   * the first.  Both loops are centered at 'p', and each has 'numVertices' vertices.  The outer 
   * loop has the given radius 'outerRadius'.  The inner loop is inset by a small distance ('gap') 
   * from the outer loop which is approximately equal to 'gapEdgeMultiple' times the edge length of
   * the outer loop.  (This allows better spatial indexing, which becomes less effective at pruning 
   * intersection candidates as the loops get closer together.)
   * 
   * <p>Caveats: The gap is actually measured to the incircle of the outer loop, and the gap is
   * clamped if necessary to prevent the inner loop from becoming vanishingly small.  (Rule of
   * thumb: to obtain a 'gapEdgeMultiple' of 'm', the loops must have approximately 7 * m
   * vertices or more.
   */
  private static List<S2Loop> getNestedLoopPair(S1Angle outerRadius, double gapEdgeMultiple, 
      int numVertices, S2Point p) {
    // The inner loop is inscribed within the incircle (maximum inscribed circle) of the outer
    // loop.
    S1Angle incircleRadius =
        S1Angle.radians(outerRadius.radians() * Math.cos(Math.PI / numVertices));
    S1Angle edgeLen = S1Angle.radians(outerRadius.radians() * (2 * Math.PI / numVertices));

    // If the edge count is too small, it may not be possible to inset the inner loop by the given
    // multiple of the edge length.  We handle this by clamping 'innerRadius' to be at least 1% of
    // 'outerRadius'.
    S1Angle innerRadius = 
        S1Angle.radians(Math.max(incircleRadius.radians() - gapEdgeMultiple * edgeLen.radians(),
            0.01 * incircleRadius.radians()));

    // Generate two loops with the same center.
    List<S2Loop> loops = new ArrayList<>(2);
    loops.add(S2Loop.makeRegularLoop(p, outerRadius, numVertices));
    loops.add(S2Loop.makeRegularLoop(p, innerRadius, numVertices));
    
    return loops;
  }

  /**
   * Returns a pair of crossing loops. The first loop in the returned list will have center 'aPoint'
   * and radius 'aRadius'.  The second will have its center along the arc containing aPoint and
   * bPoint, and it will have radius 'bRadius'.  Both loops have 'numVertices' vertices.
   */
  private static List<S2Loop> getCrossingLoopPair(S1Angle aRadius, S1Angle bRadius,
      int numVertices, S2Point aPoint, S2Point bPoint) {
    // The edges of each loop are bounded by two circles, one circumscribed around the loop (the
    // circumcircle), and the other inscribed within the loop (the incircle).  Our strategy is to
    // place the smaller loop such that its incircle crosses both circles of the larger loop.
    S1Angle maxRadius = S1Angle.radians(Math.max(aRadius.radians(), bRadius.radians()));
    S1Angle minRadius = S1Angle.radians(Math.min(aRadius.radians(), bRadius.radians()));

    // Check that the smaller loop is big enough that its incircle can span the gap between the
    // incircle and the circumcircle of the larger loop.
    double incircleFactor = Math.cos(Math.PI / numVertices);
    assertTrue(minRadius.radians() * incircleFactor > maxRadius.radians() * (1 - incircleFactor));

    // Compute the range of distances between the two loop centers such that the incircle of the
    // smaller loop crosses both circles of the larger loop.
    S1Angle minDist = S1Angle.radians(maxRadius.radians() - incircleFactor * minRadius.radians());
    S1Angle maxDist = S1Angle.radians(incircleFactor * (minRadius.radians() + maxRadius.radians()));

    // Now generate a pair of loops whose centers are separated by distances in the given range.
    // Loop orientations are chosen randomly.
    Random rand = new Random(0);
    List<S2Loop> loops = new ArrayList<>(2);
    loops.add(S2Loop.makeRegularLoop(aPoint, aRadius, numVertices));
    S1Angle angle = S1Angle.radians(
        minDist.radians() + rand.nextDouble() * (maxDist.radians() - minDist.radians()));
    S2Point bCenter = S2EdgeUtil.interpolateAtDistance(angle, aPoint, bPoint);
    loops.add(S2Loop.makeRegularLoop(bCenter, bRadius, numVertices));
    
    return loops;
  }

  /**
   * Returns a pair of disjoint loops . The loops are constructed so that it is impossible to
   * determine the relationship between them based solely on their bounds (they could be nested,
   * crossing, or disjoint).  The outer loop (1st loop in returned list) will look somewhat
   * like the outline of the  letter "C": it consists of two nested loops (the "outside shell" and
   * the "inside shell"), which each have a single edge removed and are then joined together to form
   * a single loop.  The inner loop (2nd loop in returned list) is then nested within the inside 
   * shell of the outer loop.
   * 
   * <p>The outer loop has 'numVertices' vertices split between its outside and inside shells.  The
   * radius of the outside shell is 'outerRadius', while the radius of the inside shell is
   * (0.9 * outerRadius).
   * 
   * <p>The inner loop has 'numVertices' vertices, and is separated from the inside shell of the
   * outer loop by a small distance ("gap") which is approximately equal to 'gapMultipleEdges'
   * times the edge length of the inside shell.  (See getNestedLoopPair for details.)
   */
  private static List<S2Loop> getDisjointLoopPair(S1Angle outerRadius, double gapEdgeMultiple, 
      int numVertices, S2Point p) {
    // Compute the radius of the inside shell of the outer loop, the edge length of the outer
    // shell, and finally the incircle radius of the inside shell (this is the maximum possible
    // radius of the inner loop).
    S1Angle outerInsideRadius = 
        S1Angle.radians(0.9 * outerRadius.radians() * Math.cos(2 * Math.PI / numVertices));
    S1Angle edgeLen = S1Angle.radians(outerInsideRadius.radians() * (Math.PI / numVertices));
    S1Angle incircleRadius =
        S1Angle.radians(outerInsideRadius.radians() * Math.cos(2 * Math.PI / numVertices));

    // See comments in getNestedLoopPair().
    S1Angle innerRadius = 
        S1Angle.radians(Math.max(incircleRadius.radians() - gapEdgeMultiple * edgeLen.radians(),
            0.01 * incircleRadius.radians()));

    S2Loop outerOutside = S2Loop.makeRegularLoop(p, outerRadius, Math.max(4, numVertices / 2));
    S2Loop outerInside = 
        S2Loop.makeRegularLoop(p, outerInsideRadius, Math.max(4, numVertices / 2));
    List<S2Point> vertices = 
        new ArrayList<>(outerInside.numVertices() + outerOutside.numVertices());

    // Join together the outside and inside shells to form the outer loop.
    for (int j = 0; j < outerInside.numVertices(); ++j) {
      vertices.add(outerInside.vertex(j));
    }
    Lists.reverse(vertices);
    for (int j = 0; j < outerOutside.numVertices(); ++j) {
      vertices.add(outerOutside.vertex(j));
    }
    
    return Arrays.asList(new S2Loop(vertices), S2Loop.makeRegularLoop(p, innerRadius, numVertices));
  }

  /**
   * Returns the pair of crossing loops given by getCrossingLoopPair() when using default values for
   * the radii of the two loops.  The loops will have centers 'aPoint' and 'bPoint', and
   * will have 'numVertices' vertices each.
   */
  private static List<S2Loop> getCrossingLoopPairDefault(int numVertices, S2Point aPoint,
      S2Point bPoint) {
    S1Angle aRadius = GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM);
    S1Angle bRadius = S1Angle.radians(aRadius.radians() * DEFAULT_CROSSING_RADIUS_RATIO);
    return getCrossingLoopPair(aRadius, bRadius, numVertices, aPoint, bPoint);
  }
}