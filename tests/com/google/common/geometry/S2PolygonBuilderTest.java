/*
 * Copyright 2006 Google Inc.
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

package com.google.common.geometry;

import com.google.common.annotations.GwtCompatible;
import com.google.common.base.Joiner;
import com.google.common.collect.Lists;

import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Tests for {@link S2PolygonBuilder}.
 *
 */
@GwtCompatible
public strictfp class S2PolygonBuilderTest extends GeometryTestCase {
  /** Holds the original S2Loop loglevel. */
  private Level oldLoopLevel;
  /** Holds the original S2Polygon loglevel. */
  private Level oldPolygonLevel;

  @Override
  protected void setUp() {
    super.setUp();

    // Disable the logging done by S2Loop and S2Polygon, since it will get in
    // the way of the logging we need to see from S2PolygonBuilder*java, and
    // save the levels to restore when we teardown the test.
    Logger loopLogger = Platform.getLoggerForClass(S2Loop.class);
    oldLoopLevel = loopLogger.getLevel();
    loopLogger.setLevel(Level.OFF);
    Logger polyLogger = Platform.getLoggerForClass(S2Polygon.class);
    oldPolygonLevel = polyLogger.getLevel();
    polyLogger.setLevel(Level.OFF);
  }

  @Override
  protected void tearDown() {
    // Restore the logging done by S2Loop and S2Polygon to the prior settings.
    Platform.getLoggerForClass(S2Loop.class).setLevel(oldLoopLevel);
    Platform.getLoggerForClass(S2Polygon.class).setLevel(oldPolygonLevel);
  }

  /**
   * A chain represents either a polyline or a loop, depending on whether
   * "closed" is true.
   */
  private static final class Chain {
    public final String str;
    public final boolean closed;
    public Chain(String str, boolean closed) {
      this.str = str;
      this.closed = closed;
    }
  }

  /** Test 0: No loops. */
  public void testNoLoops() {
    Chain[] chainsIn = {new Chain(null, false)};
    String[] loopsOut = {};
    runTest(0, 0, true, 0.0, 10.0, 90.0, chainsIn, loopsOut, 0, true);
  }

  /** Test 1: One loop with some extra edges. */
  public void testLoopWithExtras() {
    Chain[] chainsIn = {
        new Chain("0:0, 0:10, 10:5", true),
        new Chain("0:0, 5:5", false),
        new Chain("10:5, 20:7, 30:10, 40:15, 50:3, 60:-20", false)};
    String[] loopsOut = {"0:0, 0:10, 10:5"};
    runTest(0, 0, true, 0.0, 4.0, 15.0, chainsIn, loopsOut, 6, true);
  }

  /**
   * Test 2: One loop that has an edge removed by XORing, plus lots of extra
   * edges.
   */
  public void testLoopEdgeXorWithExtras() {
    Chain[] chainsIn = {
        new Chain("0:0, 0:10, 5:15, 10:10, 10:0", true),
        new Chain("10:10, 12:12, 14:14, 16:16, 18:18", false),
        new Chain("14:14, 14:16, 14:18, 14:20", false),
        new Chain("14:18, 16:20, 18:22", false),
        new Chain("18:12, 16:12, 14:12, 12:12", false),
        new Chain("20:18, 18:16, 16:14, 14:12", false),
        new Chain("20:14, 18:14, 16:14", false),
        new Chain("5:15, 0:10", false)};
    String[] loopsOut = {};
    runTest(0, 1, true, 0.0, 1.0, 45.0, chainsIn, loopsOut, 21, true);
  }

  /** Test 3: Three loops (two shells and one hole) that combine into one. */
  public void testOneShellTwoHoles() {
    Chain[] chainsIn = {
        new Chain("0:0, 0:10, 5:10, 10:10, 10:5, 10:0", true),
        new Chain("0:10, 0:15, 5:15, 5:10", true),
        new Chain("10:10, 5:10, 5:5, 10:5", true)};
    String[] loopsOut = {"0:0, 0:10, 0:15, 5:15, 5:10, 5:5, 10:5, 10:0"};
    // XOR.
    runTest(0, 1, true, 0.0, 4.0, 90.0, chainsIn, loopsOut, 0, true);
  }

  /**
   * Test 4: A big CCW triangle contain 3 CW triangular holes. The whole thing
   * looks like a pyramid of nine small triangles (with two extra edges).
   */
  public void testTriangleWithHoles() {
    Chain[] chainsIn = {
        new Chain("0:0, 0:2, 0:4, 0:6, 1:5, 2:4, 3:3, 2:2, 1:1", true),
        new Chain("0:2, 1:1, 1:3", true),
        new Chain("0:4, 1:3, 1:5", true),
        new Chain("1:3, 2:2, 2:4", true),
        new Chain("0:0, -1:1", false),
        new Chain("3:3, 5:5", false)};
    String[] loopsOut = {
        "0:0, 0:2, 1:1",
        "0:2, 0:4, 1:3",
        "0:4, 0:6, 1:5",
        "1:1, 1:3, 2:2",
        "1:3, 1:5, 2:4",
        "2:2, 2:4, 3:3"};
    // Directed edges required for unique result.
    runTest(-1, 0, true, 0.0, 0.9, 30.0, chainsIn, loopsOut, 2, true);
  }

  /**
   * Test 5: A square divided into four subsquares. In this case we want to
   * extract the four loops rather than taking their union. There are four extra
   * edges as well.
   */
  public void testFourSubsquaresLoops() {
    Chain[] chainsIn = {
        new Chain("0:0, 0:5, 5:5, 5:0", true),
        new Chain("0:5, 0:10, 5:10, 5:5", true),
        new Chain("5:0, 5:5, 10:5, 10:0", true),
        new Chain("5:5, 5:10, 10:10, 10:5", true),
        new Chain("0:10, 0:15, 0:20", false),
        new Chain("20:0, 15:0, 10:0", false)};
    String[] loopsOut = {
        "0:0, 0:5, 5:5, 5:0",
        "0:5, 0:10, 5:10, 5:5",
        "5:0, 5:5, 10:5, 10:0",
        "5:5, 5:10, 10:10, 10:5"};
    // Don't XOR.
    runTest(0, -1, true, 0.0, 4.0, 90.0, chainsIn, loopsOut, 4, false);
  }

  /** Test 6: Five nested loops that touch at a point. */
  public void testFiveNestedLoops() {
    Chain[] chainsIn = {
        new Chain("0:0, 0:10, 10:10, 10:0", true),
        new Chain("0:0, 1:9, 9:9, 9:1", true),
        new Chain("0:0, 2:8, 8:8, 8:2", true),
        new Chain("0:0, 3:7, 7:7, 7:3", true),
        new Chain("0:0, 4:6, 6:6, 6:4", true)};
    String[] loopsOut = {
        "0:0, 0:10, 10:10, 10:0",
        "0:0, 1:9, 9:9, 9:1",
        "0:0, 2:8, 8:8, 8:2",
        "0:0, 3:7, 7:7, 7:3",
        "0:0, 4:6, 6:6, 6:4"};
    runTest(0, 0, true, 0.0, 0.8, 5.0, chainsIn, loopsOut, 0, true);
  }

  /** Test 7: Four diamonds nested within each other touching at two points. */
  public void testFourNestedDiamonds() {
    Chain[] chainsIn = {
        new Chain("0:-20, -10:0, 0:20, 10:0", true),
        new Chain("0:10, -10:0, 0:-10, 10:0", true),
        new Chain("0:-10, -5:0, 0:10, 5:0", true),
        new Chain("0:5, -5:0, 0:-5, 5:0", true)};
    String[] loopsOut = {
        "0:-20, -10:0, 0:-10, 10:0",
        "0:-10, -5:0, 0:-5, 5:0",
        "0:5, -5:0, 0:10, 5:0",
        "0:10, -10:0, 0:20, 10:0"};
    // Directed edges required for unique result.
    runTest(-1, 0, true, 0.0, 4.0, 15.0, chainsIn, loopsOut, 0, true);
  }

  /**
   * Test 8: Seven diamonds nested within each other touching at one point
   * between each nested pair.
   */
  public void testSevenDiamonds() {
    Chain[] chainsIn = {
        new Chain("0:-70, -70:0, 0:70, 70:0", true),
        new Chain("0:-70, -60:0, 0:60, 60:0", true),
        new Chain("0:-50, -60:0, 0:50, 50:0", true),
        new Chain("0:-40, -40:0, 0:50, 40:0", true),
        new Chain("0:-30, -30:0, 0:30, 40:0", true),
        new Chain("0:-20, -20:0, 0:30, 20:0", true),
        new Chain("0:-10, -20:0, 0:10, 10:0", true)};
    String[] loopsOut = {
        "0:-70, -70:0, 0:70, 70:0",
        "0:-70, -60:0, 0:60, 60:0",
        "0:-50, -60:0, 0:50, 50:0",
        "0:-40, -40:0, 0:50, 40:0",
        "0:-30, -30:0, 0:30, 40:0",
        "0:-20, -20:0, 0:30, 20:0",
        "0:-10, -20:0, 0:10, 10:0"};
    runTest(0, 0, true, 0.0, 9.0, 4.0, chainsIn, loopsOut, 0, true);
  }

  /** Test 9: A triangle and a self-intersecting bowtie. */
  public void testTriangleAndBowtie() {
    Chain[] chainsIn = {
        new Chain("0:0, 0:10, 5:5", true),
        new Chain("0:20, 0:30, 10:20", false),
        new Chain("10:20, 10:30, 0:20", false)};
    String[] loopsOut = {"0:0, 0:10, 5:5"};
    runTest(0, 0, false, 0.0, 4.0, 45.0, chainsIn, loopsOut, 4, true);
  }

  /** Test 10: Two triangles that intersect each other. */
  public void testTwoIntersectingTriangles() {
    Chain[] chainsIn = {
        new Chain("0:0, 0:12, 6:6", true),
        new Chain("3:6, 3:18, 9:12", true)};
    String[] loopsOut = {};
    runTest(0, 0, false, 0.0, 2.0, 45.0, chainsIn, loopsOut, 6, true);
  }

  /**
   * Test 11: Four squares that combine to make a big square. The nominal edges
   * of the square are at +/-8.5 degrees in latitude and longitude. All vertices
   * except the center vertex are perturbed by up to 0.5 degrees in latitude
   * and/or longitude.
   *
   * <p>The various copies of the center vertex are misaligned by more than this
   * (i.e. they are structured as a tree where adjacent vertices are separated
   * by at most 1 degree in latitude and/or longitude) so that the clustering
   * algorithm needs more than one iteration to find them all. Note that the
   * merged position of this vertex doesn't matter because it is XORed away in
   * the output. However, it's important that all edge pairs that need to be
   * XORed are separated by no more than 'min_merge' below.
   */
  public void testFourSubsquaresUnion() {
    Chain[] chainsIn = {
        new Chain("-8:-8, -8:0", false),
        new Chain("-8:1, -8:8", false),
        new Chain("0:-9, 1:-1", false),
        new Chain("1:2, 1:9", false),
        new Chain("0:8, 2:2", false),
        new Chain("0:-2, 1:-8", false),
        new Chain("8:9, 9:1", false),
        new Chain("9:0, 8:-9", false),
        new Chain("9:-9, 0:-8", false),
        new Chain("1:-9, -9:-9", false),
        new Chain("8:0, 1:0", false),
        new Chain("-1:1, -8:0", false),
        new Chain("-8:1, -2:0", false),
        new Chain("0:1, 8:1", false),
        new Chain("-9:8, 1:8", false),
        new Chain("0:9, 8:8", false)};
    String[] loopsOut = {
        "8.5:8.5, 8.5:0.5, 8.5:-8.5, 0.5:-8.5, " +
        "-8.5:-8.5, -8.5:0.5, -8.5:8.5, 0.5:8.5"};
    // XOR, min_merge > sqrt(2), max_merge < 6.
    runTest(0, 1, true, 1.7, 5.8, 70.0, chainsIn, loopsOut, 0, true);
  }

  /**
   * Parse the vertices in {@code str}, transforms them into the given
   * {@code basis}, and adds the resulting points to {@code vertices}.
   */
  private static void getVertices(String str, Matrix3x3 basis, List<S2Point> vertices) {
    S2Polyline line = makePolyline(str);
    for (int i = 0; i < line.numVertices(); ++i) {
      vertices.add(S2Point.normalize(rotate(line.vertex(i), basis)));
    }
  }

  /**
   * Returns a fraction between 0 and 1 where small values are more likely. In
   * particular it often returns exactly 0, and often returns a fraction whose
   * logarithm is uniformly distributed over some interval.
   */
  private double smallFraction() {
    double r = rand.nextDouble();
    double u = rand.nextDouble();
    if (r < 0.3) {
      return 0.0;
    }
    if (r < 0.6) {
      return u;
    }
    return Math.pow(1e-10, u);
  }

  /**
   * Returns a copy of the given point {@code p} after rotating it by the
   * rotation matrix {@code r}.
   */
  private static S2Point rotate(S2Point p, Matrix3x3 r) {
    Matrix3x3 rotated = r.mult(new Matrix3x3(1, p.x, p.y, p.z));
    return new S2Point(rotated.get(0, 0), rotated.get(1, 0), rotated.get(2, 0));
  }

  /** Returns a random rotation matrix. */
  private Matrix3x3 getRotationMatrix() {
    List<S2Point> points = getRandomFrame();
    S2Point a = points.get(0);
    S2Point b = points.get(1);
    S2Point c = points.get(2);
    return new Matrix3x3(3,
        a.getX(), b.getX(), c.getX(),
        a.getY(), b.getY(), c.getY(),
        a.getZ(), b.getZ(), c.getZ());
  }

  /** Returns the point "x" randomly perturbed within a radius of maxPerturb. */
  private S2Point perturb(S2Point x, double maxPerturb) {
    if (maxPerturb == 0) {
      return x;
    } else {
      return samplePoint(S2Cap.fromAxisAngle(S2Point.normalize(x), S1Angle.radians(maxPerturb)));
    }
  }

  /**
   * Adds an edge from {@code v0} to {@code v1}, possibly splitting it
   * recursively up to {@code maxSplits} times, and perturbing each vertex up to
   * a distance of {@code maxPerturb}. No edge shorter than {@code minEdge} will
   * be created due to splitting.
   */
  private void addEdge(S2Point v0, S2Point v1, int maxSplits, double maxPerturb,
      double minEdge, S2PolygonBuilder builder) {
    double length = v0.angle(v1);
    if (maxSplits > 0 && rand.nextInt(2) == 0 && length >= 2 * minEdge) {
      // Choose an interpolation parameter such that the length of each
      // piece is at least minEdge.
      double f = minEdge / length;
      double t = f + (1 - 2 * f) * rand.nextDouble();

      // Now add the two sub-edges recursively.
      S2Point vmid = S2EdgeUtil.interpolate(t, v0, v1);
      addEdge(v0, vmid, maxSplits - 1, maxPerturb, minEdge, builder);
      addEdge(vmid, v1, maxSplits - 1, maxPerturb, minEdge, builder);
    } else {
      builder.addEdge(perturb(v0, maxPerturb), perturb(v1, maxPerturb));
    }
  }

  /**
   * Return true if "loop" matches any of the given candidates. The type of
   * matching depends on whether any edge splitting was done.
   */
  private static boolean findLoop(S2Loop loop, List<S2Loop> candidates, int maxSplits,
      double maxError) {
    for (S2Loop candidate : candidates) {
      if (maxSplits == 0) {
        // The two loops should match except for vertex perturbations.
        if (loop.boundaryApproxEquals(candidate, maxError)) {
          return true;
        }
      } else {
        // The two loops may have different numbers of vertices.
        if (loop.boundaryNear(candidate, maxError)) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Returns true if any loops in {@code actual} are missing from
   * {@code expected}, and log any loops from "actual" that are not present in
   * "expected".
   */
  private static boolean findMissingLoops(List<S2Loop> actual, List<S2Loop> expected,
      Matrix3x3 m, int maxSplits, double maxError, String label) {
    for (int i = 0; i < actual.size(); ++i) {
      if (!findLoop(actual.get(i), expected, maxSplits, maxError)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Transform the given edge chain to the reference frame given by
   * {@code basis}, optionally split each edge into pieces and/or perturb the
   * vertices up to the given radius, and add them to the builder.
   */
  private void addChain(Chain chain, Matrix3x3 basis, int maxSplits, double maxPerturb,
      double minEdge, S2PolygonBuilder builder) {
    List<S2Point> vertices = Lists.newArrayList();
    getVertices(chain.str, basis, vertices);
    if (chain.closed) {
      vertices.add(vertices.get(0));
    }
    for (int i = 1; i < vertices.size(); ++i) {
      addEdge(vertices.get(i - 1), vertices.get(i), maxSplits, maxPerturb, minEdge, builder);
    }
  }

  private boolean evalTristate(int state) {
    return (state > 0) ? true : (state < 0) ? false : (rand.nextDouble() > 0.5);
  }

  /**
   * Returns true if the actual number of unused edges is inconsistent with the
   * expected number of unused edges. If there are no splits, the number of
   * unused edges should match exactly. Otherwise, both values should be zero or
   * both should be non-zero.
   */
  private static boolean unexpectedUnusedEdgeCount(int numActual, int numExpected, int maxSplits) {
    if (maxSplits == 0) {
      return numActual != numExpected;
    } else {
      return (numActual > 0) != (numExpected > 0);
    }
  }

  /**
   * Print the unused edges, transformed back into their original
   * latitude-longitude space in degrees.
   */
  private static String dumpUnusedEdges(List<S2Edge> unusedEdges, Matrix3x3 m, int numExpected) {
    if (unusedEdges.size() == numExpected) {
      return "";
    }
    List<String> lines = Lists.newArrayList();
    // Get inverse to rotate coordinates back to starting location.
    Matrix3x3 inverse = m.transpose();
    lines.add(Platform.formatString("Wrong number of unused edges (%d expected, %d actual):",
        numExpected, unusedEdges.size()));
    for (int i = 0; i < unusedEdges.size(); ++i) {
      S2LatLng p0 = new S2LatLng(rotate(unusedEdges.get(i).getStart(), inverse));
      S2LatLng p1 = new S2LatLng(rotate(unusedEdges.get(i).getEnd(), inverse));
      lines.add(Platform.formatString("  [%.6f, %.6f] -> [%.6f, %.5f]",
          p0.lat().degrees(), p0.lng().degrees(),
          p1.lat().degrees(), p1.lng().degrees()));
    }
    return Joiner.on("\n").join(lines) + "\n";
  }

  /**
   * Runs a test with the given arguments.
   *
   * @param undirectedEdges +1 = undirected, -1 = directed, 0 = either one.
   * @param xorEdges +1 = XOR, -1 = don't XOR, 0 = either one.
   * @param canSplit Can edges be split for this test case.
   * @param minMerge Minimum merge distance for this test case in degrees.
   * @param maxMerge Maximum merge distance for this test case in degrees.
   * @param minVertexAngle Minimum angle in degrees between any two edges
   * <strong>after</strong> vertex merging.
   * @param chainsIn Each test case consists of a set of input loops and
   * polylines.
   * @param loopsOut The expected set of output loops, directed appropriately.
   * @param numUnusedEdges The expected number of unused edges.
   * @param makesPolygon If true, we expect a valid polygon from the loops.
   */
  private void runTest(
      final int undirectedEdges,
      final int xorEdges,
      final boolean canSplit,
      final double minMerge,
      final double maxMerge,
      final double minVertexAngle,
      final Chain[] chainsIn,
      final String[] loopsOut,
      final int numUnusedEdges,
      final boolean makesPolygon) {
    for (int iter = 0; iter < 500; ++iter) {
      // Initialize to the default options, which are changed below
      S2PolygonBuilder.Options options = S2PolygonBuilder.Options.builder()
          .setUndirectedEdges(evalTristate(undirectedEdges))
          .setXorEdges(evalTristate(xorEdges))
          .build();

      // Each test has a minimum and a maximum merge radius.  The merge
      // radius must be at least the given minimum to ensure that all expected
      // merging will take place, and it must be at most the given maximum to
      // ensure that no unexpected merging takes place.
      //
      // If the minimum and maximum values are different, we have some latitude
      // to perturb the vertices as long as the merge radius is adjusted
      // appropriately.  If "p" is the maximum perturbation radius, "m" and
      // "M" are the min/max merge radii, and "v" is the vertex merge radius
      // for this test, we require that
      //
      //       v >= m + 2*p    and    v <= M - 2*p .
      //
      // This implies that we can choose "v" in the range [m,M], and then choose
      //
      //       p <= 0.5 * min(v - m, M - v) .
      //
      // Things get more complicated when we turn on edge splicing.  Since the
      // min/max merge radii apply to vertices, we need to adjust them to ensure
      // that vertices are not accidentally spliced into nearby edges.  Recall
      // that the edge splice radius is defined as (e = v * f) where "f" is the
      // edge splice fraction.  Letting "a" be the minimum angle between two
      // edges at a vertex, we need to ensure that
      //
      //     e <= M * sin(a) - 2*p .
      //
      // The right-hand side is a lower bound on the distance from a vertex to a
      // non-incident edge.  (To simplify things, we ignore this case and fold
      // it into the case below.)
      //
      // If we also split edges by introducing new vertices, things get even
      // more complicated.  First, the vertex merge radius "v" must be chosen
      // such that
      //
      //      e >= m + 2*p    and  v <= M * sin(a) - 2*p .
      //
      // Note that the right-hand inequality now applies to "v" rather than "e",
      // since a new vertex can be introduced anywhere along a split edge.
      //
      // Finally, we need to ensure that the new edges created by splitting an
      // edge are not too short, otherwise unbounded vertex merging and/or edge
      // splicing can occur.  Letting "g" be the minimum distance (gap) between
      // vertices along a split edge, we require that
      //
      //      2 * sin(a/2) * (g - m) - 2*p >= v
      //
      // which is satisfied whenever
      //
      //      g >= m + (v + 2*p) / sin(a)
      //
      // This inequality is derived by considering two edges of length "g"
      // meeting at an angle "a", where both vertices are perturbed by distance
      // "p" toward each other, and the shared vertex is perturbed by the
      // minimum merge radius "m" along one of the two edges.

      double minMergeRadians = S1Angle.degrees(minMerge).radians();
      double maxMergeRadians = S1Angle.degrees(maxMerge).radians();
      double minSin = Math.sin(S1Angle.degrees(minVertexAngle).radians());

      // Half of the time we allow edges to be split into smaller pieces
      // (up to 5 levels, i.e. up to 32 pieces).
      int maxSplits = Math.max(0, rand.nextInt(10) - 4);
      if (!canSplit) {
        maxSplits = 0;
      }

      // We choose randomly among two different values for the edge fraction,
      // just to exercise that code.
      double edgeFraction = options.getEdgeSpliceFraction();
      double vertexMerge, maxPerturb;
      if (minSin < edgeFraction && rand.nextInt(2) == 0) {
        edgeFraction = minSin;
      }
      if (maxSplits == 0 && rand.nextInt(2) == 0) {
        // Turn off edge splicing completely.
        edgeFraction = 0;
        vertexMerge = minMergeRadians + smallFraction() * (maxMergeRadians - minMergeRadians);
        maxPerturb = 0.5 * Math.min(vertexMerge - minMergeRadians,
                                maxMergeRadians - vertexMerge);
      } else {
        // Splice edges.  These bounds also assume that edges may be split
        // (see detailed comments above).
        //
        // If edges are actually split, need to bump up the minimum merge radius
        // to ensure that split edges in opposite directions are unified.
        // Otherwise there will be tiny degenerate loops created.
        if (maxSplits > 0) {
          minMergeRadians += 1e-15;
        }
        minMergeRadians /= edgeFraction;
        maxMergeRadians *= minSin;
        assertTrue(maxMergeRadians >= minMergeRadians);

        vertexMerge = minMergeRadians + smallFraction() * (maxMergeRadians - minMergeRadians);
        maxPerturb = 0.5 * Math.min(edgeFraction * (vertexMerge - minMergeRadians),
                                maxMergeRadians - vertexMerge);
      }

      // We can perturb by any amount up to the maximum, but choosing a
      // lower maximum decreases the error bounds when checking the output.
      maxPerturb *= smallFraction();

      // This is the minimum length of a split edge to prevent unexpected
      // merging and/or splicing (the "g" value mentioned above).
      double minEdge = minMergeRadians + (vertexMerge + 2 * maxPerturb) / minSin;

      S2PolygonBuilder.Options newOptions = options.toBuilder()
          .setMergeDistance(S1Angle.radians(vertexMerge))
          .setEdgeSpliceFraction(edgeFraction)
          .setValidate(true)
          .build();
      S2PolygonBuilder builder = new S2PolygonBuilder(newOptions);

      // On each iteration we randomly rotate the test case around the sphere.
      // This causes the S2PolygonBuilder to choose different first edges when
      // trying to build loops.
      Matrix3x3 m = getRotationMatrix();
      for (Chain chain : chainsIn) {
        addChain(chain, m, maxSplits, maxPerturb, minEdge, builder);
      }
      List<S2Loop> loops = Lists.newArrayList();
      List<S2Edge> unusedEdges = Lists.newArrayList();
      if (xorEdges < 0) {
        builder.assembleLoops(loops, unusedEdges);
      } else {
        S2Polygon polygon = new S2Polygon();
        builder.assemblePolygon(polygon, unusedEdges);
        polygon.release(loops);
      }
      List<S2Loop> expected = Lists.newArrayList();
      for (String loop : loopsOut) {
        List<S2Point> vertices = Lists.newArrayList();
        getVertices(loop, m, vertices);
        expected.add(new S2Loop(vertices));
      }

      // We assume that the vertex locations in the expected output polygon
      // are separated from the corresponding vertex locations in the input
      // edges by at most half of the minimum merge radius.  Essentially
      // this means that the expected output vertices should be near the
      // centroid of the various input vertices.
      //
      // If any edges were split, we need to allow a bit more error due to
      // inaccuracies in the interpolated positions.  Similarly, if any vertices
      // were perturbed, we need to bump up the error to allow for numerical
      // errors in the actual perturbation.
      double maxError = 0.5 * minMergeRadians + maxPerturb;
      if (maxSplits > 0 || maxPerturb > 0) {
        maxError += 1e-15;
      }

      assertTrue(
          Platform.formatString(
              "During iteration %d:\n  undirected: %b\n  xor: %b\n" +
                  "  maxSplits: %d\n  maxPerturb: %.6g\n" +
                  "  vertexMerge_radius: %.6g\n  edgeSpliceFraction: %.6g\n" +
                  "  minEdge: %.6g\n  maxError: %.6g\n%s%s%s",
              iter, options.getUndirectedEdges(), options.getXorEdges(),
              maxSplits, S1Angle.radians(maxPerturb).degrees(),
              options.getMergeDistance().degrees(),
              options.getEdgeSpliceFraction(),
              S1Angle.radians(minEdge).degrees(),
              S1Angle.radians(maxError).degrees(),
              loopsToString("Expected", m, expected),
              loopsToString("Actual", m, loops),
              dumpUnusedEdges(unusedEdges, m, numUnusedEdges)),
          !findMissingLoops(loops, expected, m, maxSplits, maxError, "Actual") &&
              !findMissingLoops(expected, loops, m, maxSplits, maxError, "Expected") &&
              !unexpectedUnusedEdgeCount(unusedEdges.size(), numUnusedEdges, maxSplits) &&
              makesPolygon == new S2Polygon(loops).isValid());
    }
  }

  private static String loopsToString(String label, Matrix3x3 m, List<S2Loop> actual) {
    List<String> lines = Lists.newArrayList();
    for (int i = 0; i < actual.size(); i++) {
      lines.add(Platform.formatString("%s loop %d:", label, i));
      S2Loop loop = actual.get(i);
      Matrix3x3 inverse = m.transpose();
      for (int j = 0; j < loop.numVertices(); ++j) {
        S2LatLng ll = new S2LatLng(rotate(loop.vertex(j), inverse));
        lines.add(Platform.formatString("  [%.6f, %.6f]",
            ll.lat().degrees(), ll.lng().degrees()));
      }
    }
    return Joiner.on("\n").join(lines) + "\n";
  }

  /**
   * Verifies that {@link S2PolygonBuilder#assemblePolygon(S2Polygon,List)}
   * produces an {@code S2Polygon} with consistently merged points as long as
   * the edges are added in the same order and have the same value according to
   * {@link S2Point#equals(Object)}, even if the {@code S2Point}s are different
   * instances each time.  For example, when building an {@code S2Polygon} from
   * coordinates read from a file, most applications desire the builder to
   * always merge the same vertices, and this test verifies that it will.
   */
  public void testVertexMergeDeterminism() {
    final int NUM_TESTS = 100;

    // Make sure the points vary by reference.
    S2Point[][] points = new S2Point[NUM_TESTS][5];
    for (int i = 0; i < NUM_TESTS; i++) {
      points[i][0] = makePoint("0.4330415572436411:0.25");
      points[i][1] = makePoint("0.4330415572436411:0.75");
      points[i][2] = makePoint("0:0.4999999999999999");
      points[i][3] = makePoint("0:0.5");
      points[i][4] = makePoint("0.8660254037844386:0.5");
    }

    // Assemble NUM_TESTS polygons, making sure each result is the same as the
    // one before it.
    S2Polygon last = null;
    for (int i = 0; i < NUM_TESTS; i++) {
      S2PolygonBuilder.Options options = S2PolygonBuilder.Options.DIRECTED_XOR.toBuilder()
          .setMergeDistance(S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE)
          .build();
      S2PolygonBuilder builder = new S2PolygonBuilder(options);
      S2Polygon result = new S2Polygon();

      builder.addEdge(points[i][0], points[i][2]);
      builder.addEdge(points[i][3], points[i][1]);
      builder.addEdge(points[i][3], points[i][2]);
      builder.addEdge(points[i][1], points[i][4]);
      builder.addEdge(points[i][4], points[i][0]);

      builder.assemblePolygon(result, null);
      if (last != null) {
        assertPolygonsExactlyEqual(last, result);
      }
      last = result;
    }
  }

  private void assertPolygonsExactlyEqual(S2Polygon a, S2Polygon b) {
    assertEquals(a.numLoops(), b.numLoops());
    for (int i = 0; i < a.numLoops(); i++) {
      assertEquals(a.loop(i).numVertices(), b.loop(i).numVertices());
      for (int j = 0; j < a.loop(i).numVertices(); j++) {
        assertEquals(a.loop(i).vertex(j), b.loop(i).vertex(j));
      }
    }
  }
}
