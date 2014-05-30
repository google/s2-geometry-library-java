package com.google.common.geometry;

import static org.junit.Assert.assertTrue;

import com.google.caliper.Benchmark;
import com.google.caliper.Param;
import com.google.common.collect.Lists;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * A benchmark for {@link S2Loop}. Caliper ignores returned values, but these are computed and
 * returned so that corresponding for loops don't get optimized away.
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

  /** 
   * Ability to use methods / fields from GeometryTestCase.  Caliper insists that S2LoopBenchmark
   * not extend anything other than Object.
   */
  GeometryTestCase testUtils;

  /** 
   * Initializes random seed and lets us use GeometryTestCase methods.  We call this once at the
   * beginning of each benchmark method that needs it in order to avoid non-deterministic
   * benchmarks.  We don't use @BeforeExperiment here, because if we did, setUp() would be called
   * only before all benchmark methods of the same kind are run, which can lead to non-deterministic
   * results.
   */
  private void setUp() {
    testUtils = new GeometryTestCase();
    testUtils.setUp();
  }

  static class LoopConstructor {
    @Param({"4", "8", "64", "512", "4096", "32768", "262144"})
    int numVertices;

    GeometryTestCase testUtils;

    private void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /** Tests the speed of constructing loops from a given list of vertices. */
    @Benchmark int loopConstructor(int reps) {
      setUp();
      // Create a list of sample loops.
      List<List<S2Point>> loops = new ArrayList<>(NUM_LOOP_SAMPLES);
      for (int i = 0; i < NUM_LOOP_SAMPLES; ++i) {
        S2Loop loop = S2Loop.makeRegularLoop(testUtils.randomPoint(),
            GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM),
            numVertices);
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
   * Return a list of 'NUM_LOOP_SAMPLES' regular loops, each of 'numVertices' vertices, with
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

    GeometryTestCase testUtils;

    private void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /** Tests the speed of isValid() on regular loops. */
    @Benchmark int isValid(int reps) {
      setUp();
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
    @Param({"4", "8", "16", "32", "64", "128", "256", "512", "1024", "2048", "8192", 
      "32768", "262144"})
    int numVertices;

    static final int NUM_QUERIES_PER_LOOP = 100;

    GeometryTestCase testUtils;

    private void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /** 
     * Tests the speed of checking whether a regular loop contains a point.  The query points for a
     * loop are chosen so that they all lie in the loop's bounding rectangle (to avoid the 
     * quick-rejection code path.
     */
    @Benchmark int containsPoint (int reps) {
      setUp();
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
    private int numVertices = 1024;

    GeometryTestCase testUtils;

    private void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /**
     * Constructs a loop of 'numVertices' vertices and calls contains(S2Point) 'numCalls' times.
     * The whole process is repeated 'reps' times.
     */
    @Benchmark int constructAndContains(int reps) {
      setUp();
      S2Point query = testUtils.randomPoint();
      // Build a regular loop either around the query point or its antipode, based on whether the
      // loop should contain the query point.
      S2Loop loop = S2Loop.makeRegularLoop(
          shouldContain ? query : S2Point.neg(query),
              GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM), numVertices);
      List<S2Point> loopVertices = new ArrayList<>(numVertices);
      for (int v = 0; v < numVertices; ++v) {
        loopVertices.add(loop.vertex(v));
      }

      int numContains = 0;
      for (int r = reps; r > 0; --r) {
        S2Loop tempLoop = new S2Loop(loopVertices);
        for (int n = numCalls; n > 0; --n) {
          if (tempLoop.contains(query)) {
            numContains++;
          }
        }
      }

      return numContains;
    }
  }

  /**
   * Use an interface so that compareBoundary(), intersects(), and contains() can be used as 
   * arguments.
   */
  private interface SpatialRelation {
    int testRelation(S2Loop a, S2Loop b);
  }

  private static class CompareBoundary implements SpatialRelation {
    @Override public int testRelation(S2Loop a, S2Loop b) {
      return a.compareBoundary(b);
    }
  }

  private static class Intersects implements SpatialRelation {
    @Override public int testRelation(S2Loop a, S2Loop b) {
      return a.intersects(b) ? 1 : 0;
    }
  }

  private static class Contains implements SpatialRelation {
    @Override public int testRelation(S2Loop a, S2Loop b) {
      return a.contains(b) ? 1 : 0;
    }
  }

  /** 
   * The loop relation methods (contains, intersects, compareBoundary) scale similarly with the
   * number of vertices.  So to save time, we only benchmark compareBoundary() on a full range of
   * vertex counts.  For contains() and intersects() we just run a single test using 
   * DEFAULT_NUM_VERTICES.
   */
  static class CompareBoundaryBenchmarks {
    @Param({"8", "64", "512", "4096", "32768"})
    int numVertices;

    GeometryTestCase testUtils;

    private void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /** Benchmarks compareBoundary() when one loop contains the other. */
    @Benchmark void compareBoundaryContains(int reps) {
      setUp();
      relationContains(reps, new CompareBoundary(), 1, numVertices, randomSamples(testUtils));
    }

    /** Benchmarks compareBoundary() when the loops cross each other. */
    @Benchmark void compareBoundaryCrosses(int reps) {
      setUp();
      relationCrosses(reps, new CompareBoundary(), 0, numVertices, randomSamples(testUtils),
          randomSamples(testUtils));
    }

    /** Benchmarks compareBoundary() when the loops are disjoint. */
    @Benchmark void compareBoundaryDisjoint(int reps) {
      setUp();
      relationDisjoint(reps, new CompareBoundary(), -1, numVertices, randomSamples(testUtils));
    }
  }

  /** Benchmarks contains() when one loop contains the other. */
  @Benchmark void containsContains(int reps) {
    setUp();
    relationContains(reps, new Contains(), 1, DEFAULT_NUM_VERTICES, randomSamples(testUtils));
  }

  /** Benchmarks contains() when the loops cross each other. */
  @Benchmark void containsCrosses(int reps) {
    setUp();
    relationCrosses(reps, new Contains(), 0, DEFAULT_NUM_VERTICES, randomSamples(testUtils),
        randomSamples(testUtils));
  }

  /** Benchmarks contains() when the loops are disjoint. */
  @Benchmark void containsDisjoint (int reps) {
    setUp();
    relationDisjoint(reps, new Contains(), 1, DEFAULT_NUM_VERTICES, randomSamples(testUtils));
  }  

  /** Benchmarks intersects() when one loop contains the other. */
  @Benchmark void intersectsContains(int reps) {
    setUp();
    List<S2Point> randPoints = new ArrayList<>(NUM_LOOP_SAMPLES);
    for (int i = 0; i < NUM_LOOP_SAMPLES; ++i) {
      randPoints.add(testUtils.randomPoint());
    }
    relationContains(reps, new Intersects(), 1, DEFAULT_NUM_VERTICES, randomSamples(testUtils));
  }

  /** Benchmarks intersects() when the two loops cross. */
  @Benchmark void intersectsCrosses(int reps) {
    setUp();
    relationCrosses(reps, new Intersects(), 1, DEFAULT_NUM_VERTICES, randomSamples(testUtils),
        randomSamples(testUtils));
  }

  /** Benchmarks intersects() when the loops are disjoint. */
  @Benchmark void intersectsDisjoint(int reps) {
    setUp();
    relationDisjoint(reps, new Intersects(), 1, DEFAULT_NUM_VERTICES, randomSamples(testUtils));
  }

  static class ContainsVsRadiusMeters {
    /** Radius of the outer loop. */
    @Param({"1", "8", "64", "1024", "65536"})
    int meters;

    GeometryTestCase testUtils;

    private void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /**
     * Benchmark contains() on nested loops as a function of the loop radius.  Performance
     * degrades on very small loops because S2.expensiveCCW() is invoked more often.  'meters' is
     * the radius of the outer loop.
     */
    @Benchmark void containsVsRadiusMeters(int reps) {
      setUp();
      List<S2Loop> aLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
      List<S2Loop> bLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
      getNestedLoopPairs(GeometryTestCase.metersToAngle(meters), DEFAULT_NESTED_GAP_MULTIPLE,
          DEFAULT_NUM_VERTICES, aLoops, bLoops, randomSamples(testUtils));
      loopRelation(reps, new Contains(), 1, aLoops, bLoops);
    }
  }

  static class ContainsVsGapEdgeMultiple {
    /** See comments for getNestedLoopPairs(). */
    @Param({"0", "1", "8", "64"})
    int gapEdgeMultiple;

    GeometryTestCase testUtils;

    private void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /**
     * Benchmark contains() on nested loops as a function of the distance between the loops
     * (expressed as a multiple of the loop edge length).  Performance degrades as the loops get
     * very close together because spatial indexing is not as effective at pruning the
     * intersection candidates.
     */
    @Benchmark void containsVsGapEdgeMultiple (int reps) {
      setUp();
      List<S2Loop> aLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
      List<S2Loop> bLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
      getNestedLoopPairs(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM), gapEdgeMultiple,
          DEFAULT_NUM_VERTICES, aLoops, bLoops, randomSamples(testUtils));
      loopRelation(reps, new Contains(), 1, aLoops, bLoops);
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

    GeometryTestCase testUtils;

    private void setUp() {
      testUtils = new GeometryTestCase();
      testUtils.setUp();
    }

    /** 
     * Benchmark intersects() on crossing loops as a function of the relative sizes of the two
     * loops (e.g., one looop radius much larger tha nthe other).  Performance of spatial indexing
     * can degrade when one loop is much larger than the other.
     */
    @Benchmark void intersectsCrossesVsLogRatio(int reps) {
      setUp();
      List<S2Loop> aLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
      List<S2Loop> bLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
      S1Angle aRadius = GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM);
      S1Angle bRadius = S1Angle.radians(aRadius.radians() * Math.pow(2.0, logRadiusRatio));

      getCrossingLoopPairs(aRadius, bRadius, DEFAULT_NUM_VERTICES, aLoops, bLoops,
          randomSamples(testUtils), randomSamples(testUtils));
      loopRelation(reps, new Intersects(), 1, aLoops, bLoops);
    }
  }

  /**
   * Adds an outer and inner loop to 'outerLoops' and 'innerLoops' respectively, for each central
   * vertex given in 'randPoints'.  All loops have 'numVertices' vertices.  Each outer loop has the
   * given 'outerRadius'.  Each inner loop is inset by a small distance ('gap') from the outer loop
   * which is approximately equal to 'gapEdgeMultiple' times the edge length of the outer loop.
   * (This allows better spatial indexing, which becomes less effective at pruning intersection
   * candidates as the loops get closer together.)  Loops are otherwise positioned randomly
   * (but repeatably) on the sphere.
   * 
   * <p>Caveats: The gap is actually measured to the incircle of the outer loop, and the gap is
   * clamped if necessary to prevent the inner loop from becoming vanishingly small.  (Rule of
   * thumb: to obtain a 'gapEdgeMultiple' of 'm', the loops must have approximately 7 * m
   * vertices or more.
   */
  private static void getNestedLoopPairs(S1Angle outerRadius, double gapEdgeMultiple, 
      int numVertices, List<S2Loop> outerLoops, List<S2Loop> innerLoops, List<S2Point> randPoints) {
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

    for (int i = 0; i < randPoints.size(); ++i) {
      // Generate two loops with the same center but with random different orientations.
      S2Point p = randPoints.get(i);
      outerLoops.add(S2Loop.makeRegularLoop(p, outerRadius, numVertices));
      innerLoops.add(S2Loop.makeRegularLoop(p, innerRadius, numVertices));
    }
  }

  /**
   * Adds a pair of loops to aLoops and bLoops whose boundaries cross each other, for each point
   * in aPoints. For the i-th pair of loops, the center of one loop will be at aPoints.get(i), and
   * the center of the other loop will be lie along the arc containing aPoints.get(i) and 
   * bPoints.get(i).  All loops have numVertices vertices.  The aLoops must all have radius aRadius,
   * and similarly for bLoops. Loops are otherwise positioned randomly (but repeatably) on the
   * sphere. Requires aPoints and bPoints lists to be equal length.
   */
  private static void getCrossingLoopPairs (S1Angle aRadius, S1Angle bRadius, int numVertices,
      List<S2Loop> aLoops, List<S2Loop> bLoops, List<S2Point> aPoints, List<S2Point> bPoints) {
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

    // Now generate pairs of loops whose centers are separated by distances in the given range.
    // Loop orientations are chosen randomly.
    Random rand = new Random(0);
    for (int i = 0; i < aPoints.size(); ++i) {
      S2Point p = aPoints.get(i);
      aLoops.add(S2Loop.makeRegularLoop(p, aRadius, numVertices));
      S1Angle angle = S1Angle.radians(
          minDist.radians() + rand.nextDouble() * (maxDist.radians() - minDist.radians()));
      S2Point bCenter = S2EdgeUtil.interpolateAtDistance(angle, p, bPoints.get(i));
      bLoops.add(S2Loop.makeRegularLoop(bCenter, bRadius, numVertices));
    }
  }

  /**
   * Adds an outer and inner loop to 'outerLoops' and 'innerLoops' respectively, for each central
   * vertex given in 'randPoints'.  The disjoint pairs are constructed so that it is impossible to
   * determine the relationship between the two loops based solely on their bounds (they could be 
   * nested, crossing, or disjoint).  Each 'outer loop' looks somewhat like the outline of the 
   * letter "C": it consists of two nested loops (the "outside shell" and the "inside shell"), 
   * which each have a single edge removed and are then joined together to form a single loop.  The
   * 'inner loop' is then nested within the inside shell of the outer loop.
   * 
   * <p>The outer loop has 'numVertices' vertices split between its outside and inside shells.  The
   * radius of the outside shell is 'outerRadius', while the radius of the inside shell is
   * (0.9 * outerRadius).
   * 
   * <p>The inner loop has 'numVertices' vertices, and is separated from the inside shell of the
   * outer loop by a small distance ("gap") which is approximately equal to 'gapMultipleEdges'
   * times the edge length of the inside shell.  (See getNestedLoopPairs for details.)
   */
  private static void getDisjointLoopPairs(S1Angle outerRadius, double gapEdgeMultiple, 
      int numVertices, List<S2Loop> outerLoops, List<S2Loop> innerLoops, List<S2Point> randPoints) {
    // Compute the radius of the inside shell of the outer loop, the edge length of the outer
    // shell, and finally the incircle radius of the inside shell (this is the maximum possible
    // radius of the inner loop).
    S1Angle outerInsideRadius = 
        S1Angle.radians(0.9 * outerRadius.radians() * Math.cos(2 * Math.PI / numVertices));
    S1Angle edgeLen = S1Angle.radians(outerInsideRadius.radians() * (Math.PI / numVertices));
    S1Angle incircleRadius =
        S1Angle.radians(outerInsideRadius.radians() * Math.cos(2 * Math.PI / numVertices));

    // See comments in getNestedLoopPairs().
    S1Angle innerRadius = 
        S1Angle.radians(Math.max(incircleRadius.radians() - gapEdgeMultiple * edgeLen.radians(),
            0.01 * incircleRadius.radians()));
    for (int i = 0; i < randPoints.size(); ++i) {
      S2Point p = randPoints.get(i);
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
      outerLoops.add(new S2Loop(vertices));

      // Now construct the inner loop, which has the same center.
      innerLoops.add(S2Loop.makeRegularLoop(p, innerRadius, numVertices));
    }
  }

  /**
   * Takes in two lists of loops of equal size and tests the speed of the relation method
   * (e.g., compareBoundary(), intersects(), or contains()).
   */
  private static void loopRelation(int reps, SpatialRelation method,
      @SuppressWarnings("unused") int expectedResult, List<S2Loop> aLoops, List<S2Loop> bLoops) {
    for (int i = 0, r = reps; r > 0; --r) {
      // TODO(eengle) Get assertEquals working for loopRelation.
      // assertEquals(expectedResult, method.testRelation(aLoops.get(i), bLoops.get(i)));
      method.testRelation(aLoops.get(i), bLoops.get(i));
      if (++i == aLoops.size()) {
        i = 0;
      }
    }
  }

  /**
   * Tests the speed of the given S2Loop loop relation method (e.g., intersects) on a collection of
   * loop pairs where one loop contains the other loop (i.e., nested loops).  Verifies that the
   * method always returns 'expectedResult'.
   */
  private static void relationContains(int reps, SpatialRelation method, int expectedResult,
      int numVertices, List<S2Point> randomPoints) {
    List<S2Loop> aLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
    List<S2Loop> bLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
    getNestedLoopPairs(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM), DEFAULT_NESTED_GAP_MULTIPLE, 
        numVertices, aLoops, bLoops, randomPoints);
    loopRelation(reps, method, expectedResult, aLoops, bLoops);
  }

  /**
   * Tests the speed of the given S2Loop loop relation method (e.g., intersects) on a collection of
   * loop pairs where the loop boundaries cross each other.  Verifies that the method always returns
   * 'expectedResult'.
   */
  private static void relationCrosses(int reps, SpatialRelation method, int expectedResult, 
      int numVertices, List<S2Point> aPoints, List<S2Point> bPoints) {
    List<S2Loop> aLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
    List<S2Loop> bLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
    S1Angle aRadius = GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM);
    S1Angle bRadius = S1Angle.radians(aRadius.radians() * DEFAULT_CROSSING_RADIUS_RATIO);
    getCrossingLoopPairs(aRadius, bRadius, numVertices, aLoops, bLoops, aPoints, bPoints);
    loopRelation(reps, method, expectedResult, aLoops, bLoops);
  }

  /**
   * Tests the speed of the given S2Loop loop relation method (e.g., intersects) on a collection of
   * disjoint loop pairs.  Verifies that the method always returns 'expectedResult'.
   */
  private static void relationDisjoint(int reps, SpatialRelation method, int expectedResult,
      int numVertices, List<S2Point> randPoints) {
    List<S2Loop> aLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
    List<S2Loop> bLoops = new ArrayList<>(NUM_LOOP_SAMPLES);
    getDisjointLoopPairs(GeometryTestCase.kmToAngle(DEFAULT_RADIUS_KM), DEFAULT_NESTED_GAP_MULTIPLE,
        numVertices, aLoops, bLoops, randPoints);
    loopRelation(reps, method, expectedResult, aLoops, bLoops);
  }

  private static final List<S2Point> randomSamples(GeometryTestCase testUtils) {
    List<S2Point> points = new ArrayList<>(NUM_LOOP_SAMPLES);
    for (int i = 0; i < NUM_LOOP_SAMPLES; ++i) {
      points.add(testUtils.randomPoint());
    }
    return points;
  }
}