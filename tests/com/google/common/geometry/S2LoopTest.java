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
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;

import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.logging.Logger;

/**
 * Tests for {@link S2Loop}.
 *
 *  Note that testLoopRelations2() is suppressed because it fails in corner
 * cases due to a problem with S2.robustCCW().
 *
 */
@GwtCompatible
public strictfp class S2LoopTest extends GeometryTestCase {
  private static final Logger log = Platform.getLoggerForClass(S2LoopTest.class);

  // A stripe that slightly over-wraps the equator.
  private S2Loop candyCane = makeLoop("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70");

  // A small clockwise loop in the northern & eastern hemisperes.
  private S2Loop smallNeCw = makeLoop("35:20, 45:20, 40:25");

  // Loop around the north pole at 80 degrees.
  private S2Loop arctic80 = makeLoop("80:-150, 80:-30, 80:90");

  // Loop around the south pole at 80 degrees.
  private S2Loop antarctic80 = makeLoop("-80:120, -80:0, -80:-120");

  // A completely degenerate triangle along the equator that RobustCCW()
  // considers to be CCW.
  private S2Loop lineTriangle = makeLoop("0:1, 0:2, 0:3");

  // A nearly-degenerate CCW chevron near the equator with very long sides
  // (about 80 degrees).  Note that the precision is less than the C++
  // equivalent due to Java strictfp precision issues, but it is precise
  // enough to still have incredibly small area.
  private S2Loop skinnyChevron = makeLoop("0:0, -1e-80:80, 0:1e-80, 1e-80:80");

  // The northern hemisphere, defined using two pairs of antipodal points.
  private S2Loop northHemi = makeLoop("0:-180, 0:-90, 0:0, 0:90");

  // The northern hemisphere, defined using three points 120 degrees apart.
  private S2Loop northHemi3 = makeLoop("0:-180, 0:-60, 0:60");

  // The western hemisphere, defined using two pairs of antipodal points.
  private S2Loop westHemi = makeLoop("0:-180, -90:0, 0:0, 90:0");

  // The "near" hemisphere, defined using two pairs of antipodal points.
  private S2Loop nearHemi = makeLoop("0:-90, -90:0, 0:90, 90:0");

  // A diamond-shaped loop around the point 0:180.
  private S2Loop loopA = makeLoop("0:178, -1:180, 0:-179, 1:-180");

  // Another diamond-shaped loop around the point 0:180.
  private S2Loop loopB = makeLoop("0:179, -1:180, 0:-178, 1:-180");

  // A shape gotten from A by adding a triangle to one edge, and
  // subtracting a triangle from the opposite edge.
  private S2Loop loopC = makeLoop("0:178, 0:180, -1:180, 0:-179, 1:-179, 1:-180");

  // A shape gotten from A by adding a triangle to one edge, and
  // adding another triangle to the opposite edge.
  private S2Loop loopD = makeLoop("0:178, -1:178, -1:180, 0:-179, 1:-179, 1:-180");

  //   3------------2
  //   |            |               ^
  //   |  7-8  b-c  |               |
  //   |  | |  | |  |      Latitude |
  //   0--6-9--a-d--1               |
  //   |  | |       |               |
  //   |  f-e       |               +----------->
  //   |            |                 Longitude
  //   4------------5
  //
  // Important: It is not okay to skip over collinear vertices when
  // defining these loops (e.g. to define loop E as "0,1,2,3") because S2
  // uses symbolic perturbations to ensure that no three vertices are
  // *ever* considered collinear (e.g., vertices 0, 6, 9 are not
  // collinear).  In other words, it is unpredictable (modulo knowing the
  // details of the symbolic perturbations) whether 0123 contains 06123,
  // for example.
  //
  // Loop E:  0,6,9,a,d,1,2,3
  // Loop F:  0,4,5,1,d,a,9,6
  // Loop G:  0,6,7,8,9,a,b,c,d,1,2,3
  // Loop H:  0,6,f,e,9,a,b,c,d,1,2,3
  // Loop I:  7,6,f,e,9,8
  private S2Loop loopE = makeLoop("0:30, 0:34, 0:36, 0:39, 0:41, 0:44, 30:44, 30:30");
  private S2Loop loopF = makeLoop("0:30, -30:30, -30:44, 0:44, 0:41, 0:39, 0:36, 0:34");
  private S2Loop loopG = makeLoop(
      "0:30, 0:34, 10:34, 10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30");
  private S2Loop loopH = makeLoop(
      "0:30, 0:34, -10:34, -10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30");
  private S2Loop loopI = makeLoop("10:34, 0:34, -10:34, -10:36, 0:36, 10:36");

  // The intersection of A and B.
  private S2Loop aIntersectB = makeLoop("0:179, -1:180, 0:-179, 1:-180");

  // The union of A and B.
  private S2Loop aUnionB = makeLoop("0:178, -1:180, 0:-178, 1:-180");

  // A minus B (concave)
  private S2Loop aMinusB = makeLoop("0:178, -1:180, 0:179, 1:-180");

  // B minus A (concave)
  private S2Loop bMinusA = makeLoop("0:-179, -1:180, 0:-178, 1:-180");

  // A self-crossing loop with a duplicated vertex
  private S2Loop bowtie = makeLoop("0:0, 2:0, 1:1, 0:2, 2:2, 1:1");

  private S2Loop southHemi = invert(northHemi);
  private S2Loop eastHemi = invert(westHemi);
  private S2Loop farHemi = invert(nearHemi);

  private static S2Loop invert(S2Loop loop) {
    S2Loop copy = new S2Loop(loop);
    copy.invert();
    return copy;
  }

  /** Verifies that a is a loose bound of b, with boundaries at most maxError apart. */
  private static void assertBoundsEqual(S2LatLngRect a, S2LatLngRect b, S2LatLng maxError) {
    assertTrue(a.approxEquals(b, maxError));
    assertTrue(a.contains(b));
  }

  /** Verifies that a is a loose bound of b, with boundaries at most maxError apart. */
  private static void assertBoundsEqual(R1Interval a, R1Interval b, double radianError) {
    assertTrue(a.approxEquals(b, radianError));
    assertTrue(a.contains(b));
  }

  public void testBounds() {
    assertTrue(S2Loop.empty().getRectBound().isEmpty());
    assertTrue(S2Loop.full().getRectBound().isFull());
    assertTrue(candyCane.getRectBound().lng().isFull());
    assertTrue(candyCane.getRectBound().latLo().degrees() < -20);
    assertTrue(candyCane.getRectBound().latHi().degrees() > 10);
    assertTrue(smallNeCw.getRectBound().isFull());
    // Permit small latitude error since we pad the loop bounds a bit, but longitude must be exact.
    double latError = 2E-15;
    assertBoundsEqual(arctic80.getRectBound(), new S2LatLngRect(S2LatLng.fromDegrees(80, -180),
        S2LatLng.fromDegrees(90, 180)), S2LatLng.fromRadians(latError, 0));
    assertBoundsEqual(antarctic80.getRectBound(), new S2LatLngRect(S2LatLng.fromDegrees(-90, -180),
        S2LatLng.fromDegrees(-80, 180)), S2LatLng.fromRadians(latError, 0));

    arctic80.invert();
    // The highest latitude of each edge is attained at its midpoint.
    S2Point mid = S2Point.mul(S2Point.add(arctic80.vertex(0), arctic80.vertex(1)), 0.5);
    assertDoubleNear(arctic80.getRectBound().latHi().radians(), new S2LatLng(mid).lat().radians());
    arctic80.invert();

    assertTrue(southHemi.getRectBound().lng().isFull());
    assertBoundsEqual(southHemi.getRectBound().lat(), new R1Interval(-S2.M_PI_2, 0), latError);
  }

  public void testFastConstructor() {
    List<S2Point> vertices = Lists.newArrayList();
    S2LatLngRect bound = parseVertices("-80:120, -80:0, -80:-120", vertices);
    S2Loop loop = S2Loop.newLoopWithTrustedDetails(vertices, false, bound);
    assertTrue(loop.isValid());
    assertEquals(0, loop.depth());
  }

  public void testAreaCentroid() {
    assertEquals(0.0, S2Loop.empty().getArea());
    assertEquals(4 * S2.M_PI, S2Loop.full().getArea());
    assertEquals(null, S2Loop.empty().getCentroid());
    assertEquals(null, S2Loop.full().getCentroid());

    assertDoubleNear(northHemi.getArea(), 2 * S2.M_PI);
    assertDoubleNear(eastHemi.getArea(), 2 * S2.M_PI);

    // Construct spherical caps of random height, and approximate their boundary
    // with closely spaces vertices. Then check that the area and centroid are
    // correct.

    for (int i = 0; i < 100; ++i) {
      // Choose a coordinate frame for the spherical cap.
      S2Point x = randomPoint();
      S2Point y = S2Point.normalize(S2Point.crossProd(x, randomPoint()));
      S2Point z = S2Point.normalize(S2Point.crossProd(x, y));

      // Given two points at latitude phi and whose longitudes differ by dtheta,
      // the geodesic between the two points has a maximum latitude of
      // atan(tan(phi) / cos(dtheta/2)). This can be derived by positioning
      // the two points at (-dtheta/2, phi) and (dtheta/2, phi).
      //
      // We want to position the vertices close enough together so that their
      // maximum distance from the boundary of the spherical cap is kMaxDist.
      // Thus we want fabs(atan(tan(phi) / cos(dtheta/2)) - phi) <= kMaxDist.
      double kMaxDist = TestPlatform.S2_LOOP_TEST_MAX_ANGLE_DELTA;
      double height = 2 * rand.nextDouble();
      double phi = Math.asin(1 - height);
      double maxDtheta =
          2 * Math.acos(Math.tan(Math.abs(phi)) / Math.tan(Math.abs(phi) + kMaxDist));
      maxDtheta = Math.min(S2.M_PI, maxDtheta); // At least 3 vertices.

      List<S2Point> vertices = Lists.newArrayList();
      for (double theta = 0; theta < 2 * S2.M_PI; theta += rand.nextDouble() * maxDtheta) {

        S2Point xCosThetaCosPhi = S2Point.mul(x, (Math.cos(theta) * Math.cos(phi)));
        S2Point ySinThetaCosPhi = S2Point.mul(y, (Math.sin(theta) * Math.cos(phi)));
        S2Point zSinPhi = S2Point.mul(z, Math.sin(phi));

        S2Point sum = S2Point.add(S2Point.add(xCosThetaCosPhi, ySinThetaCosPhi), zSinPhi);

        vertices.add(sum);
      }

      S2Loop loop = new S2Loop(vertices);
      S2AreaCentroid areaCentroid = loop.getAreaAndCentroid();

      double area = loop.getArea();
      S2Point centroid = loop.getCentroid();
      double expectedArea = 2 * S2.M_PI * height;
      assertTrue(areaCentroid.getArea() == area);
      assertTrue(centroid.equals(areaCentroid.getCentroid()));
      assertTrue(Math.abs(area - expectedArea) <= 2 * S2.M_PI * kMaxDist);

      // high probability
      assertTrue(Math.abs(area - expectedArea) >= 0.01 * kMaxDist);

      S2Point expectedCentroid = S2Point.mul(z, expectedArea * (1 - 0.5 * height));

      assertTrue(S2Point.sub(centroid, expectedCentroid).norm() <= 2 * kMaxDist);
    }
  }

  private static S2Loop rotate(S2Loop loop) {
    List<S2Point> vertices = Lists.newArrayList();
    for (int i = 1; i <= loop.numVertices(); ++i) {
      vertices.add(loop.vertex(i));
    }
    return new S2Loop(vertices);
  }

  public void testContains() {
    // Check the full and empty loops have the correct containment relationship
    // with the special "vertex" that defines them.
    assertFalse(S2Loop.empty().contains(S2Loop.EMPTY_VERTEX));
    assertTrue(S2Loop.full().contains(S2Loop.FULL_VERTEX));

    assertTrue(candyCane.contains(S2LatLng.fromDegrees(5, 71).toPoint()));
    for (int i = 0; i < 4; ++i) {
      assertTrue(northHemi.contains(new S2Point(0, 0, 1)));
      assertTrue(!northHemi.contains(new S2Point(0, 0, -1)));
      assertTrue(!southHemi.contains(new S2Point(0, 0, 1)));
      assertTrue(southHemi.contains(new S2Point(0, 0, -1)));
      assertTrue(!westHemi.contains(new S2Point(0, 1, 0)));
      assertTrue(westHemi.contains(new S2Point(0, -1, 0)));
      assertTrue(eastHemi.contains(new S2Point(0, 1, 0)));
      assertTrue(!eastHemi.contains(new S2Point(0, -1, 0)));
      northHemi = rotate(northHemi);
      southHemi = rotate(southHemi);
      eastHemi = rotate(eastHemi);
      westHemi = rotate(westHemi);
    }

    // This code checks each cell vertex is contained by exactly one of
    // the adjacent cells.
    for (int level = 0; level < 3; ++level) {
      List<S2Loop> loops = Lists.newArrayList();
      List<S2Point> loopVertices = Lists.newArrayList();
      Set<S2Point> points = Sets.newHashSet();
      for (S2CellId id = S2CellId.begin(level); !id.equals(S2CellId.end(level)); id = id.next()) {
        S2Cell cell = new S2Cell(id);
        points.add(cell.getCenter());
        for (int k = 0; k < 4; ++k) {
          loopVertices.add(cell.getVertex(k));
          points.add(cell.getVertex(k));
        }
        loops.add(new S2Loop(loopVertices));
        loopVertices.clear();
      }
      for (S2Point point : points) {
        int count = 0;
        for (int j = 0; j < loops.size(); ++j) {
          if (loops.get(j).contains(point)) {
            ++count;
          }
        }
        assertEquals(count, 1);
      }
    }
  }

  public void testContainsMatchesRobustCrossing() {
    // This test demonstrates a former incompatibility between robustCrossing()
    // and contains(S2Point).  It constructs an S2Cell-based loop L and
    // an edge E from Origin to a0 that crosses exactly one edge of L.  Yet
    // previously, contains() returned false for both endpoints of E.
    //
    // The reason for the bug was that the loop bound was sometimes too tight.
    // The contains() code for a0 bailed out early because a0 was found not to
    // be inside the bound of L.

    // Start with a cell that ends up producing the problem.
    S2CellId cellId = S2CellId.fromPoint(new S2Point(1, 1, 1)).parent(21);

    S2Cell[] children = {new S2Cell(), new S2Cell(), new S2Cell(), new S2Cell()};
    new S2Cell(cellId).subdivide(children);

    S2Point[] points = new S2Point[4];
    for (int i = 0; i < 4; ++i) {
      // Note extra normalization. GetCenter() is already normalized.
      // The test results will no longer be inconsistent if the extra
      // Normalize() is removed.
      points[i] = S2Point.normalize(children[i].getCenter());
    }

    S2Loop loop = new S2Loop(Arrays.asList(points));

    // Get a vertex from a grandchild cell.  Mathematically, a0 should be the
    // same as points[0], but rounding errors make it slightly different.
    // +---------------+---------------+
    // |               |               |
    // |    points[3]  |   points[2]   |
    // |       v       |       v       |
    // |       +-------+------ +       |
    // |       |       |       |       |
    // |       |       |       |       |
    // |       |       |       |       |
    // +-------+-------+-------+-------+
    // |       |       |       |       |
    // |       |    <----------------------- grandchild_cell
    // |       |       |       |       |
    // |       +-------+------ +       |
    // |       ^       |       ^       | <-- cell
    // | points[0]/a0  |     points[1] |
    // |               |               |
    // +---------------+---------------+

    S2Cell grandchildCell = new S2Cell(cellId.child(0).child(2));
    S2Point a0 = grandchildCell.getVertex(0);

    // The edge from a0 to the origin crosses one boundary.
    assertEquals(-1, S2EdgeUtil.robustCrossing(a0, S2.origin(), loop.vertex(0), loop.vertex(1)));
    assertEquals(1, S2EdgeUtil.robustCrossing(a0, S2.origin(), loop.vertex(1), loop.vertex(2)));
    assertEquals(-1, S2EdgeUtil.robustCrossing(a0, S2.origin(), loop.vertex(2), loop.vertex(3)));
    assertEquals(-1, S2EdgeUtil.robustCrossing(a0, S2.origin(), loop.vertex(3), loop.vertex(4)));

    // Contains should return false for the origin, and true for a0.
    assertFalse(loop.contains(S2.origin()));
    assertTrue(loop.contains(a0));

    // Since a0 is inside the loop, it should be inside the bound.
    S2LatLngRect bound = loop.getRectBound();
    assertTrue(bound.contains(a0));
  }

  private S2CellId advance(S2CellId id, int n) {
    while (id.isValid() && --n >= 0) {
      id = id.next();
    }
    return id;
  }

  private S2Loop makeCellLoop(S2CellId begin, S2CellId end) {
    // Construct a CCW polygon whose boundary is the union of the cell ids
    // in the range [begin, end). We add the edges one by one, removing
    // any edges that are already present in the opposite direction.
    Map<S2Point, Set<S2Point>> edges = Maps.newHashMap();
    for (S2CellId id = begin; !id.equals(end); id = id.next()) {
      S2Cell cell = new S2Cell(id);
      for (int k = 0; k < 4; ++k) {
        S2Point a = cell.getVertex(k);
        S2Point b = cell.getVertex((k + 1) & 3);
        if (edges.get(b) == null) {
          edges.put(b, Sets.<S2Point>newHashSet());
        }

        // if a is in b's set, remove it and remove b's set if it's empty
        // otherwise, add b to a's set
        if (!edges.get(b).remove(a)) {
          if (edges.get(a) == null) {
            edges.put(a, Sets.<S2Point>newHashSet());
          }
          edges.get(a).add(b);
        } else if (edges.get(b).isEmpty()) {
          edges.remove(b);
        }
      }
    }

    // The remaining edges form a single loop. We simply follow it starting
    // at an arbitrary vertex and build up a list of vertices.
    List<S2Point> vertices = Lists.newArrayList();
    S2Point p = edges.keySet().iterator().next();
    while (!edges.isEmpty()) {
      assertEquals(1, edges.get(p).size());
      S2Point next = edges.get(p).iterator().next();
      vertices.add(p);
      edges.remove(p);
      p = next;
    }

    return new S2Loop(vertices);
  }

  /** Given a pair of loops where A contains B, check various identities. */
  private static void checkOneNestedPair(S2Loop a, S2Loop b) {
    assertTrue(a.contains(b));
    assertEquals(a.boundaryEquals(b), b.contains(a));
    assertEquals(!b.isEmpty(), a.intersects(b));
    assertEquals(!b.isEmpty(), b.intersects(a));
  }

  /** Given a pair of disjoint loops A and B, check various identities. */
  private static void checkOneDisjointPair(S2Loop a, S2Loop b) {
    assertFalse(a.intersects(b));
    assertFalse(b.intersects(a));
    assertEquals(b.isEmpty(), a.contains(b));
    assertEquals(a.isEmpty(), b.contains(a));
  }

  /** Given loops A and B whose union covers the sphere, check various identities. */
  private static void checkOneCoveringPair(S2Loop a, S2Loop b) {
    assertEquals(a.isFull(), a.contains(b));
    assertEquals(b.isFull(), b.contains(a));
    S2Loop a1 = invert(a);
    boolean complementary = a1.boundaryEquals(b);
    assertEquals(!complementary, a.intersects(b));
    assertEquals(!complementary, b.intersects(a));
  }

  /**
   * Given loops A and B such that both A and its complement intersect both B and its complement,
   * check various identities.
   */
  private static void checkOneOverlappingPair(S2Loop a, S2Loop b) {
    assertFalse(a.contains(b));
    assertFalse(b.contains(a));
    assertTrue(a.intersects(b));
    assertTrue(b.intersects(a));
  }

  /**
   * Given a pair of loops where A contains B, test various identities involving A, B, and their
   * complements.
   */
  private static void checkNestedPair(S2Loop a, S2Loop b) {
    S2Loop a1 = invert(a);
    S2Loop b1 = invert(b);
    checkOneNestedPair(a, b);
    checkOneNestedPair(b1, a1);
    checkOneDisjointPair(a1, b);
    checkOneCoveringPair(a, b1);
  }

  /**
   * Given a pair of disjoint loops A and B, test various identities involving A, B, and their
   * complements.
   */
  private static void checkDisjointPair(S2Loop a, S2Loop b) {
    S2Loop a1 = new S2Loop(a);
    a1.invert();
    checkNestedPair(a1, b);
  }

  /**
   * Given loops A and B whose union covers the sphere, test various identities involving A, B, and
   * their complements.
   */
  private static void checkCoveringPair(S2Loop a, S2Loop b) {
    S2Loop b1 = new S2Loop(b);
    b1.invert();
    checkNestedPair(a, b1);
  }

  /**
   * Given loops A and B such that both A and its complement intersect both B and its complement,
   * test various identities involving these four loops.
   */
  private static void checkOverlappingPair(S2Loop a, S2Loop b) {
    S2Loop a1 = invert(a);
    S2Loop b1 = invert(b);
    checkOneOverlappingPair(a, b);
    checkOneOverlappingPair(a1, b1);
    checkOneOverlappingPair(a1, b);
    checkOneOverlappingPair(a, b1);
  }

  /** A contains B. */
  private static final int CONTAINS = 0x1;
  /** B contains A. */
  private static final int CONTAINED = 0x2;
  /** A and B are disjoint (intersection is empty.) */
  private static final int DISJOINT = 0x4;
  /** (A union B) covers the entire sphere. */
  private static final int COVERS = 0x8;

  /** Returns true if any of the bits in needles are set in haystack. */
  private static final boolean test(int haystack, int needles) {
    return (haystack & needles) != 0;
  }

  /**
   * Verifies the relationship between two loops A and B.  "flags" is the set of RelationFlags that
   * apply.  "shared_edge" means that the loops share at least one edge (possibly reversed.)
   */
  private static void checkRelationImpl(S2Loop a, S2Loop b, int flags, boolean sharedEdge) {
    if (test(flags, CONTAINS)) {
      checkNestedPair(a, b);
    }
    if (test(flags, CONTAINED)) {
      checkNestedPair(b, a);
    }
    if (test(flags, COVERS)) {
      checkCoveringPair(a, b);
    }
    if (test(flags, DISJOINT)) {
      checkDisjointPair(a, b);
    } else if (!(test(flags, (CONTAINS | CONTAINED | COVERS)))) {
      checkOverlappingPair(a, b);
    }
    if (!sharedEdge && (test(flags, (CONTAINS | CONTAINED | DISJOINT)))) {
      assertEquals(a.contains(b), a.containsNested(b));
    }
    // A contains the boundary of B if either A contains B, or the two loops
    // contain each other's boundaries and there are no shared edges (since at
    // least one such edge must be reversed, and therefore is not considered to
    // be contained according to the rules of CompareBoundary).
    int comparison = 0;
    if (test(flags, CONTAINS) || ((test(flags, COVERS) && !sharedEdge))) {
      comparison = 1;
    }
    // Similarly, A excludes the boundary of B if either A and B are disjoint,
    // or B contains A and there are no shared edges (since A is considered to
    // contain such edges according to the rules of CompareBoundary).
    if (test(flags, DISJOINT) || (test(flags, CONTAINED) && !sharedEdge)) {
      comparison = -1;
    }
    // CompareBoundary requires that neither loop is empty.
    if (!a.isEmpty() && !b.isEmpty()) {
      assertEquals(comparison, a.compareBoundary(b));
    }
  }

  private void checkRelation(S2Loop a, S2Loop b, int flags, boolean sharedEdge) {
    try {
      checkRelationImpl(a, b, flags, sharedEdge);
    } catch (AssertionError e) {
      System.err.println("args "  + a + ", " + b);
      throw e;
    }
  }

  public void testLoopRelations() {
    // Check full and empty relationships with normal loops and each other.
    checkRelation(S2Loop.full(), S2Loop.full(), CONTAINS | CONTAINED | COVERS, true);
    checkRelation(S2Loop.full(), northHemi, CONTAINS | COVERS, false);
    checkRelation(S2Loop.full(), S2Loop.empty(), CONTAINS | DISJOINT | COVERS, false);
    checkRelation(northHemi, S2Loop.full(), CONTAINED | COVERS, false);
    checkRelation(northHemi, S2Loop.empty(), CONTAINS | DISJOINT, false);
    checkRelation(S2Loop.empty(), S2Loop.full(), CONTAINED | DISJOINT | COVERS, false);
    checkRelation(S2Loop.empty(), northHemi, CONTAINED | DISJOINT, false);
    checkRelation(S2Loop.empty(), S2Loop.empty(), CONTAINS | CONTAINED | DISJOINT, false);

    checkRelation(northHemi, northHemi, CONTAINS | CONTAINED, true);
    checkRelation(northHemi, southHemi, DISJOINT | COVERS, true);
    checkRelation(northHemi, eastHemi, 0, false);
    checkRelation(northHemi, arctic80, CONTAINS, false);
    checkRelation(northHemi, antarctic80, DISJOINT, false);
    checkRelation(northHemi, candyCane, 0, false);

    // We can't compare north_hemi3 vs. north_hemi or south_hemi because the
    // result depends on the "simulation of simplicity" implementation details.
    checkRelation(northHemi3, northHemi3, CONTAINS | CONTAINED, true);
    checkRelation(northHemi3, eastHemi, 0, false);
    checkRelation(northHemi3, arctic80, CONTAINS, false);
    checkRelation(northHemi3, antarctic80, DISJOINT, false);
    checkRelation(northHemi3, candyCane, 0, false);

    checkRelation(southHemi, northHemi, DISJOINT | COVERS, true);
    checkRelation(southHemi, southHemi, CONTAINS | CONTAINED, true);
    checkRelation(southHemi, farHemi, 0, false);
    checkRelation(southHemi, arctic80, DISJOINT, false);
    checkRelation(southHemi, antarctic80, CONTAINS, false);
    checkRelation(southHemi, candyCane, 0, false);

    checkRelation(candyCane, northHemi, 0, false);
    checkRelation(candyCane, southHemi, 0, false);
    checkRelation(candyCane, arctic80, DISJOINT, false);
    checkRelation(candyCane, antarctic80, DISJOINT, false);
    checkRelation(candyCane, candyCane, CONTAINS | CONTAINED, true);

    checkRelation(nearHemi, westHemi, 0, false);

    checkRelation(smallNeCw, southHemi, CONTAINS, false);
    checkRelation(smallNeCw, westHemi, CONTAINS, false);

    checkRelation(smallNeCw, northHemi, COVERS, false);
    checkRelation(smallNeCw, eastHemi, COVERS, false);

    checkRelation(loopA, loopA, CONTAINS | CONTAINED, true);
    checkRelation(loopA, loopB, 0, false);
    checkRelation(loopA, aIntersectB, CONTAINS, true);
    checkRelation(loopA, aUnionB, CONTAINED, true);
    checkRelation(loopA, aMinusB, CONTAINS, true);
    checkRelation(loopA, bMinusA, DISJOINT, true);

    checkRelation(loopB, loopA, 0, false);
    checkRelation(loopB, loopB, CONTAINS | CONTAINED, true);
    checkRelation(loopB, aIntersectB, CONTAINS, true);
    checkRelation(loopB, aUnionB, CONTAINED, true);
    checkRelation(loopB, aMinusB, DISJOINT, true);
    checkRelation(loopB, bMinusA, CONTAINS, true);

    checkRelation(aIntersectB, loopA, CONTAINED, true);
    checkRelation(aIntersectB, loopB, CONTAINED, true);
    checkRelation(aIntersectB, aIntersectB, CONTAINS | CONTAINED, true);
    checkRelation(aIntersectB, aUnionB, CONTAINED, false);
    checkRelation(aIntersectB, aMinusB, DISJOINT, true);
    checkRelation(aIntersectB, bMinusA, DISJOINT, true);

    checkRelation(aUnionB, loopA, CONTAINS, true);
    checkRelation(aUnionB, loopB, CONTAINS, true);
    checkRelation(aUnionB, aIntersectB, CONTAINS, false);
    checkRelation(aUnionB, aUnionB, CONTAINS | CONTAINED, true);
    checkRelation(aUnionB, aMinusB, CONTAINS, true);
    checkRelation(aUnionB, bMinusA, CONTAINS, true);

    checkRelation(aMinusB, loopA, CONTAINED, true);
    checkRelation(aMinusB, loopB, DISJOINT, true);
    checkRelation(aMinusB, aIntersectB, DISJOINT, true);
    checkRelation(aMinusB, aUnionB, CONTAINED, true);
    checkRelation(aMinusB, aMinusB, CONTAINS | CONTAINED, true);
    checkRelation(aMinusB, bMinusA, DISJOINT, false);

    checkRelation(bMinusA, loopA, DISJOINT, true);
    checkRelation(bMinusA, loopB, CONTAINED, true);
    checkRelation(bMinusA, aIntersectB, DISJOINT, true);
    checkRelation(bMinusA, aUnionB, CONTAINED, true);
    checkRelation(bMinusA, aMinusB, DISJOINT, false);
    checkRelation(bMinusA, bMinusA, CONTAINS | CONTAINED, true);
  }

  /**
   * Make sure the relations are correct if the loop crossing happens on two ends of a shared
   * boundary segment.
   */
  public void testLoopRelationsWhenSameExceptPiecesStickingOutAndIn() {
    checkRelation(loopA, loopC, 0, true);
    checkRelation(loopC, loopA, 0, true);
    checkRelation(loopA, loopD, CONTAINED, true);
    checkRelation(loopD, loopA, CONTAINS, true);
    checkRelation(loopE, loopF, DISJOINT, true);
    checkRelation(loopE, loopG, CONTAINS, true);
    checkRelation(loopE, loopH, 0, true);
    checkRelation(loopE, loopI, 0, false);
    checkRelation(loopF, loopG, DISJOINT, true);
    checkRelation(loopF, loopH, 0, true);
    checkRelation(loopF, loopI, 0, false);
    checkRelation(loopG, loopH, CONTAINED, true);
    checkRelation(loopH, loopG, CONTAINS, true);
    checkRelation(loopG, loopI, DISJOINT, true);
    checkRelation(loopH, loopI, CONTAINS, true);
  }

  /**
   * TODO(user, ericv) Fix this test. It fails sporadically.
   * <p>
   * The problem is not in this test, it is that
   * {@link S2#robustCCW(S2Point, S2Point, S2Point)} currently requires
   * arbitrary-precision arithmetic to be truly robust. That means it can give
   * the wrong answers in cases where we are trying to determine edge
   * intersections.
   * <p>
   * It seems the strictfp modifier here in java (required for correctness in
   * other areas of the library) restricts the size of temporary registers,
   * causing us to lose some of the precision that the C++ version gets.
   * <p>
   * This test fails when it randomly chooses a cell loop with nearly colinear
   * edges. That's where S2.robustCCW provides the wrong answer. Note that there
   * is an attempted workaround in {@link S2Loop#isValid()}, but it
   * does not cover all cases.
   */
  public void suppressedTestLoopRelations2() {
    // Construct polygons consisting of a sequence of adjacent cell ids
    // at some fixed level. Comparing two polygons at the same level
    // ensures that there are no T-vertices.
    for (int iter = 0; iter < TestPlatform.S2_LOOP_TEST_LOOP_RELATIONS2_ITERATIONS; ++iter) {
      long num = rand.nextLong();
      S2CellId begin = new S2CellId(num | 1);
      if (!begin.isValid()) {
        continue;
      }
      begin = begin.parent((int) Math.round(rand.nextDouble() * S2CellId.MAX_LEVEL));
      S2CellId aBegin = advance(begin, skewed(6));
      S2CellId aEnd = advance(aBegin, skewed(6) + 1);
      S2CellId bBegin = advance(begin, skewed(6));
      S2CellId bEnd = advance(bBegin, skewed(6) + 1);
      if (!aEnd.isValid() || !bEnd.isValid()) {
        continue;
      }

      S2Loop a = makeCellLoop(aBegin, aEnd);
      S2Loop b = makeCellLoop(bBegin, bEnd);
      boolean contained = (aBegin.lessOrEquals(bBegin) && bEnd.lessOrEquals(aEnd));
      boolean intersects = (aBegin.lessThan(bEnd) && bBegin.lessThan(aEnd));
      log.info(
          "Checking " + a.numVertices() + " vs. " + b.numVertices() + ", contained = " + contained
              + ", intersects = " + intersects);

      assertEquals(contained, a.contains(b));
      assertEquals(intersects, a.intersects(b));
    }
  }

  public void testBoundsForLoopContainment() {
    // To reliably test whether one loop contains another, the bounds of the
    // outer loop are expanded slightly.  This test constructs examples where
    // this expansion is necessary and verifies that it is sufficient.
    for (int iter = 0; iter < 1000; ++iter) {
      // We construct a triangle ABC such that A,B,C are nearly colinear, B is
      // the point of maximum latitude, and the edge AC passes very slightly
      // below B (i.e., ABC is CCW).
      S2Point b = S2Point.normalize(S2Point.add(randomPoint(), S2Point.Z_POS));
      S2Point v = S2Point.normalize(S2Point.crossProd(b, S2Point.Z_POS));
      S2Point a = S2EdgeUtil.interpolate(rand.nextDouble(), S2Point.neg(v), b);
      S2Point c = S2EdgeUtil.interpolate(rand.nextDouble(), b, v);
      if (S2.robustCCW(a, b, c) < 0) {
        --iter;
        continue;
      }
      // Now construct another point D directly below B, and create two loops
      // ABCD and ACD.
      S2Point d = S2Point.normalize(new S2Point(b.getX(), b.getY(), 0));
      List<S2Point> vertices = ImmutableList.of(c, d, a, b);  // Reordered for convenience
      S2Loop outer = new S2Loop(vertices);
      S2Loop inner = new S2Loop(vertices.subList(0, 3));
      // Now because the bounds calculation is less accurate when the maximum is
      // attained along an edge (rather than at a vertex), sometimes the inner
      // loop will have a *larger* bounding box than the outer loop.  We look
      // only for those cases.
      if (outer.getRectBound().contains(inner.getRectBound())) {
        --iter;
        continue;
      }
      assertTrue(outer.contains(inner));
    }
  }

  private static void testNear(String a_str, String b_str, double max_error, boolean expected) {
    S2Loop a = makeLoop(a_str);
    S2Loop b = makeLoop(b_str);
    assertEquals(expected, a.boundaryNear(b, max_error));
    assertEquals(expected, b.boundaryNear(a, max_error));
  }

  public void testBoundaryNear() {
    double degree = S1Angle.degrees(1).radians();

    testNear("0:0, 0:10, 5:5",
        "0:0.1, -0.1:9.9, 5:5.2",
        0.5 * degree, true);
    testNear("0:0, 0:3, 0:7, 0:10, 3:7, 5:5",
        "0:0, 0:10, 2:8, 5:5, 4:4, 3:3, 1:1",
        1e-3, true);

    // All vertices close to some edge, but not equivalent.
    testNear("0:0, 0:2, 2:2, 2:0",
        "0:0, 1.9999:1, 0:2, 2:2, 2:0",
        0.5 * degree, false);

    // Two triangles that backtrack a bit on different edges.  A simple
    // greedy matching algorithm would fail on this example.
    String t1 = "0.1:0, 0.1:1, 0.1:2, 0.1:3, 0.1:4, 1:4, 2:4, 3:4, " +
        "2:4.1, 1:4.1, 2:4.2, 3:4.2, 4:4.2, 5:4.2";
    String t2 = "0:0, 0:1, 0:2, 0:3, 0.1:2, 0.1:1, 0.2:2, 0.2:3, " +
        "0.2:4, 1:4.1, 2:4, 3:4, 4:4, 5:4";
    testNear(t1, t2, 1.5 * degree, true);
    testNear(t1, t2, 0.5 * degree, false);
  }

  static void checkEmptyFullSnapped(S2Loop loop, int level) {
    assertTrue(loop.isEmptyOrFull());
    S2CellId cellid = S2CellId.fromPoint(loop.vertex(0)).parent(level);
    S2Loop loop2 = new S2Loop(ImmutableList.of(cellid.toPoint()));
    assertTrue(loop.boundaryEquals(loop2));
    assertTrue(loop.boundaryApproxEquals(loop2));
    assertTrue(loop.boundaryNear(loop2));
  }

  // Test converting the empty/full loops to S2LatLng representations.  (We
  // don't bother testing E5/E6/E7 because that test is less demanding.)
  static void checkEmptyFullLatLng(S2Loop loop) {
    assertTrue(loop.isEmptyOrFull());
    S2Loop loop2 = new S2Loop(ImmutableList.of(new S2LatLng(loop.vertex(0)).toPoint()));
    assertTrue(loop.boundaryEquals(loop2));
    assertTrue(loop.boundaryApproxEquals(loop2));
    assertTrue(loop.boundaryNear(loop2));
  }

  static void checkEmptyFullConversions(S2Loop loop) {
    checkEmptyFullSnapped(loop, S2CellId.MAX_LEVEL);
    checkEmptyFullSnapped(loop, 1);  // Worst case for approximation
    checkEmptyFullSnapped(loop, 0);
    checkEmptyFullLatLng(loop);
  }

  public void testEmptyFullLossyConversions() {
    // Verify that the empty and full loops can be encoded lossily.
    checkEmptyFullConversions(S2Loop.empty());
    checkEmptyFullConversions(S2Loop.full());
  }

  /**
   * Tests that nearly colinear points pass S2Loop.isValid()
   */
  public void testRoundingError() {
    S2Point a = new S2Point(-0.9190364081111774, 0.17231932652084575, 0.35451111445694833);
    S2Point b = new S2Point(-0.92130667053206, 0.17274500072476123, 0.3483578383756171);
    S2Point c = new S2Point(-0.9257244057938284, 0.17357332608634282, 0.3360158106235289);
    S2Point d = new S2Point(-0.9278712595449962, 0.17397586116468677, 0.32982923679138537);

    assertTrue(S2Loop.isValid(Lists.newArrayList(a, b, c, d)));
  }

  /**
   * Tests {@link S2Loop#isValid()}.
   */
  public void testIsValid() {
    assertTrue(loopA.isValid());
    assertTrue(loopB.isValid());
    assertFalse(bowtie.isValid());
  }

  /**
   * Tests {@link S2Loop#compareTo(S2Loop)}.
   */
  public void testComparisons() {
    S2Loop abc = makeLoop("0:1, 0:2, 1:2");
    S2Loop abcd = makeLoop("0:1, 0:2, 1:2, 1:1");
    S2Loop abcde = makeLoop("0:1, 0:2, 1:2, 1:1, 1:0");
    assertTrue(abc.compareTo(abcd) < 0);
    assertTrue(abc.compareTo(abcde) < 0);
    assertTrue(abcd.compareTo(abcde) < 0);
    assertTrue(abcd.compareTo(abc) > 0);
    assertTrue(abcde.compareTo(abc) > 0);
    assertTrue(abcde.compareTo(abcd) > 0);

    S2Loop bcda = makeLoop("0:2, 1:2, 1:1, 0:1");
    assertEquals(0, abcd.compareTo(bcda));
    assertEquals(0, bcda.compareTo(abcd));

    S2Loop wxyz = makeLoop("10:11, 10:12, 11:12, 11:11");
    assertTrue(abcd.compareTo(wxyz) > 0);
    assertTrue(wxyz.compareTo(abcd) < 0);
  }

  public void testGetDistance() {
    // Error margin since we're doing numerical computations
    double epsilon = 1e-15;

    // A square with (lat,lng) vertices (0,1), (1,1), (1,2) and (0,2)
    // Tests the case where the shortest distance is along a normal to an edge,
    // onto a vertex
    S2Loop s1 = makeLoop("0:1, 1:1, 1:2, 0:2");

    // A square with (lat,lng) vertices (-1,1), (1,1), (1,2) and (-1,2)
    // Tests the case where the shortest distance is along a normal to an edge,
    // not onto a vertex
    S2Loop s2 = makeLoop("-1:1, 1:1, 1:2, -1:2");

    // A diamond with (lat,lng) vertices (1,0), (2,1), (3,0) and (2,-1)
    // Test the case where the shortest distance is NOT along a normal to an
    // edge
    S2Loop s3 = makeLoop("1:0, 2:1, 3:0, 2:-1");

    // All the vertices should be distance 0
    for (int i = 0; i < s1.numVertices(); i++) {
      assertEquals(0d, s1.getDistance(s1.vertex(i)).radians(), epsilon);
    }

    // A point on one of the edges should be distance 0
    assertEquals(0d, s1.getDistance(S2LatLng.fromDegrees(0.5, 1).toPoint()).radians(), epsilon);

    // In all three cases, the closest point to the origin is (0,1), which is at
    // a distance of 1 degree.
    // Note: all of these are intentionally distances measured along the
    // equator, since that makes the math significantly simpler. Otherwise, the
    // distance wouldn't actually be 1 degree.
    S2Point origin = S2LatLng.fromDegrees(0, 0).toPoint();
    assertEquals(1d, s1.getDistance(origin).degrees(), epsilon);
    assertEquals(1d, s2.getDistance(origin).degrees(), epsilon);
    assertEquals(1d, s3.getDistance(origin).degrees(), epsilon);
  }

  /**
   * Check that the turning angle is *identical* when the vertex order is rotated, and that the
   * sign is inverted when the vertices are reversed.
   */
  private static void checkTurningAngleInvariants(S2Loop loop) {
    double expected = loop.getTurningAngle();
    S2Loop loop_copy = new S2Loop(loop);
    for (int i = 0; i < loop.numVertices(); ++i) {
      loop_copy.invert();
      assertEquals(-expected, loop_copy.getTurningAngle());
      loop_copy.invert();
      loop_copy = rotate(loop_copy);
      assertEquals(expected, loop_copy.getTurningAngle());
    }
  }

  public void testTurningAngle() {
    assertEquals(2 * S2.M_PI, S2Loop.empty().getTurningAngle());
    assertEquals(-2 * S2.M_PI, S2Loop.full().getTurningAngle());

    assertDoubleNear(0, northHemi3.getTurningAngle(), 1e-15);
    checkTurningAngleInvariants(northHemi3);

    assertDoubleNear(0, westHemi.getTurningAngle(), 1e-15);
    checkTurningAngleInvariants(westHemi);

    // We don't have an easy way to estimate the turning angle of this loop, but
    // we can still check that the expected invariants hold.
    checkTurningAngleInvariants(candyCane);

    assertEquals(2 * S2.M_PI, lineTriangle.getTurningAngle());
    checkTurningAngleInvariants(lineTriangle);

    assertEquals(2 * S2.M_PI, skinnyChevron.getTurningAngle());
    checkTurningAngleInvariants(skinnyChevron);
  }

  /**
   * Checks that if a loop is normalized, it doesn't contain a point outside of it, and vice versa.
   */
  private static void checkNormalizeAndContains(S2Loop loop) {
    S2Point p = makePoint("40:40");
    S2Loop flip = new S2Loop(loop);

    flip.invert();
    assertTrue(loop.isNormalized() ^ loop.contains(p));
    assertTrue(flip.isNormalized() ^ flip.contains(p));
    assertTrue(loop.isNormalized() ^ flip.isNormalized());

    flip.normalize();
    assertFalse(flip.contains(p));
  }

  public void testNormalizedCompatibleWithContains() {
    checkNormalizeAndContains(lineTriangle);
    checkNormalizeAndContains(skinnyChevron);
  }

  public void testIsOriginInside() {
    assertEquals(eastHemi.contains(S2.origin()), eastHemi.isOriginInside());
    assertEquals(arctic80.contains(S2.origin()), arctic80.isOriginInside());
    assertFalse(S2Loop.empty().hasInterior());
    assertTrue(S2Loop.full().hasInterior());
    assertTrue(eastHemi.hasInterior());
  }

  /**
   * This function is useful for debugging.
   */
  @SuppressWarnings("unused")
  private void dumpCrossings(S2Loop loop) {

    System.out.println("Ortho(v1): " + S2.ortho(loop.vertex(1)));
    System.out.println("Contains(kOrigin): " + loop.contains(S2.origin()) + "\n");
    for (int i = 1; i <= loop.numVertices(); ++i) {
      S2Point a = S2.ortho(loop.vertex(i));
      S2Point b = loop.vertex(i - 1);
      S2Point c = loop.vertex(i + 1);
      S2Point o = loop.vertex(i);
      Platform.printf(System.out, "Vertex %d: [%.17g, %.17g, %.17g], "
          + "%d%dR=%d, %d%d%d=%d, R%d%d=%d, inside: %b\n",
          i,
          loop.vertex(i).x,
          loop.vertex(i).y,
          loop.vertex(i).z,
          i - 1,
          i,
          S2.robustCCW(b, o, a),
          i + 1,
          i,
          i - 1,
          S2.robustCCW(c, o, b),
          i,
          i + 1,
          S2.robustCCW(a, o, c),
          S2.orderedCCW(a, b, c, o));
    }
    for (int i = 0; i < loop.numVertices() + 2; ++i) {
      S2Point orig = S2.origin();
      S2Point dest;
      if (i < loop.numVertices()) {
        dest = loop.vertex(i);
        Platform.printf(System.out, "Origin->%d crosses:", i);
      } else {
        dest = new S2Point(0, 0, 1);
        if (i == loop.numVertices() + 1) {
          orig = loop.vertex(1);
        }
        Platform.printf(System.out, "Case %d:", i);
      }
      for (int j = 0; j < loop.numVertices(); ++j) {
        System.out.println(
            " " + S2EdgeUtil.edgeOrVertexCrossing(orig, dest, loop.vertex(j), loop.vertex(j + 1)));
      }
      System.out.println();
    }
    for (int i = 0; i <= 2; i += 2) {
      Platform.printf(System.out, "Origin->v1 crossing v%d->v1: ", i);
      S2Point a = S2.ortho(loop.vertex(1));
      S2Point b = loop.vertex(i);
      S2Point c = S2.origin();
      S2Point o = loop.vertex(1);
      Platform.printf(System.out, "%d1R=%d, M1%d=%d, R1M=%d, crosses: %b\n",
          i,
          S2.robustCCW(b, o, a),
          i,
          S2.robustCCW(c, o, b),
          S2.robustCCW(a, o, c),
          S2EdgeUtil.edgeOrVertexCrossing(c, o, b, a));
    }
  }
}
