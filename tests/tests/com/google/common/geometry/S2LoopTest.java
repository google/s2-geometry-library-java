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

import static com.google.common.geometry.S2.M_PI_2;
import static com.google.common.geometry.TestDataGenerator.makeLoop;
import static com.google.common.geometry.TestDataGenerator.makePoint;
import static com.google.common.geometry.TestDataGenerator.parseVertices;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.acos;
import static java.lang.Math.asin;
import static java.lang.Math.cos;
import static java.lang.Math.min;
import static java.lang.Math.sin;
import static java.lang.Math.tan;
import static java.lang.Math.toDegrees;
import static junit.framework.TestCase.assertNull;

import com.google.common.base.Predicate;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.io.BaseEncoding;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * Tests for {@link S2Loop}.
 *
 * @author shakusa@google.com (Steven Hakusa) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
public strictfp class S2LoopTest extends GeometryTestCase {

  private enum TestLoops {
    // A stripe that slightly over-wraps the equator.
    CANDY_CANE(makeLoop("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70")),

    // A small clockwise loop in the northern & eastern hemisperes.
    SMALL_NE_CW(makeLoop("35:20, 45:20, 40:25")),

    // Loop around the north pole at 80 degrees.
    ARCTIC_80(makeLoop("80:-150, 80:-30, 80:90")),

    // Loop around the south pole at 80 degrees.
    ANTARCTIC_80(makeLoop("-80:120, -80:0, -80:-120")),

    // A completely degenerate triangle along the equator that RobustCCW() considers to be CCW.
    LINE_TRIANGLE(makeLoop("0:1, 0:2, 0:3")),

    // A nearly-degenerate CCW chevron near the equator with very long sides (about 80 degrees).
    // Note that the precision is less than the C++ equivalent due to Java strictfp precision
    // issues, but it is precise enough to still have incredibly small area.
    SKINNY_CHEVRON(makeLoop("0:0, -1e-80:80, 0:1e-80, 1e-80:80")),

    // The northern hemisphere, defined using two pairs of antipodal points.
    NORTH_HEMI(makeLoop("0:-180, 0:-90, 0:0, 0:90")),

    // The northern hemisphere, defined using three points 120 degrees apart.
    NORTH_HEMI3(makeLoop("0:-180, 0:-60, 0:60")),

    // The western hemisphere, defined using two pairs of antipodal points.
    WEST_HEMI(makeLoop("0:-180, -90:0, 0:0, 90:0")),

    // The "near" hemisphere, defined using two pairs of antipodal points.
    NEAR_HEMI(makeLoop("0:-90, -90:0, 0:90, 90:0")),

    // A diamond-shaped loop around the point 0:180.
    LOOP_A(makeLoop("0:178, -1:180, 0:-179, 1:-180")),

    // Another diamond-shaped loop around the point 0:180.
    LOOP_B(makeLoop("0:179, -1:180, 0:-178, 1:-180")),

    // A shape gotten from A by adding a triangle to one edge, and subtracting a triangle from the
    // opposite edge.
    LOOP_C(makeLoop("0:178, 0:180, -1:180, 0:-179, 1:-179, 1:-180")),

    // A shape gotten from A by adding a triangle to one edge, and adding another triangle to the
    // opposite edge.
    LOOP_D(makeLoop("0:178, -1:178, -1:180, 0:-179, 1:-179, 1:-180")),

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
    // Important: It is not okay to skip over collinear vertices when defining these loops (e.g. to
    // define loop E as "0,1,2,3") because S2 uses symbolic perturbations to ensure that no three
    // vertices are *ever* considered collinear (e.g., vertices 0, 6, 9 are not collinear).  In
    // other words, it is unpredictable (modulo knowing the details of the symbolic perturbations)
    // whether 0123 contains 06123, for example.
    //
    // Loop E:  0,6,9,a,d,1,2,3
    // Loop F:  0,4,5,1,d,a,9,6
    // Loop G:  0,6,7,8,9,a,b,c,d,1,2,3
    // Loop H:  0,6,f,e,9,a,b,c,d,1,2,3
    // Loop I:  7,6,f,e,9,8
    LOOP_E(makeLoop("0:30, 0:34, 0:36, 0:39, 0:41, 0:44, 30:44, 30:30")),
    LOOP_F(makeLoop("0:30, -30:30, -30:44, 0:44, 0:41, 0:39, 0:36, 0:34")),
    LOOP_G(
        makeLoop("0:30, 0:34, 10:34, 10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30")),
    LOOP_H(
        makeLoop("0:30, 0:34, -10:34, -10:36, 0:36, 0:39, 10:39, 10:41, 0:41, 0:44, 30:44, 30:30")),
    LOOP_I(makeLoop("10:34, 0:34, -10:34, -10:36, 0:36, 10:36")),

    // Like loop_a, but the vertices are at leaf cell centers.
    SNAPPED_LOOP_A(makeLoop("0:178, -1:180, 0:-179, 1:-180", S2CellId.MAX_LEVEL)),

    // The intersection of A and B.
    A_INTERSECT_B(makeLoop("0:179, -1:180, 0:-179, 1:-180")),

    // The union of A and B.
    A_UNION_B(makeLoop("0:178, -1:180, 0:-178, 1:-180")),

    // A minus B (concave)
    A_MINUS_B(makeLoop("0:178, -1:180, 0:179, 1:-180")),

    // B minus A (concave)
    B_MINUS_A(makeLoop("0:-179, -1:180, 0:-178, 1:-180")),

    // A self-crossing loop with a duplicated vertex
    BOWTIE(makeLoop("0:0, 2:0, 1:1, 0:2, 2:2, 1:1")),

    SOUTH_HEMI(invert(NORTH_HEMI.loop)),
    EAST_HEMI(invert(WEST_HEMI.loop)),
    FAR_HEMI(invert(NEAR_HEMI.loop));

    final S2Loop loop;

    TestLoops(S2Loop loop) {
      this.loop = loop;
    }
  }

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

  public void testIsEmptyOrIsFull() {
    assertTrue(S2Loop.empty().isEmpty());
    assertFalse(S2Loop.empty().isFull());

    assertFalse(S2Loop.full().isEmpty());
    assertTrue(S2Loop.full().isFull());

    assertFalse(TestLoops.CANDY_CANE.loop.isEmpty());
    assertFalse(TestLoops.CANDY_CANE.loop.isFull());
  }

  public void testBounds() {
    assertTrue(S2Loop.empty().getRectBound().isEmpty());
    assertTrue(S2Loop.full().getRectBound().isFull());
    assertTrue(TestLoops.CANDY_CANE.loop.getRectBound().lng().isFull());
    assertTrue(TestLoops.CANDY_CANE.loop.getRectBound().latLo().degrees() < -20);
    assertTrue(TestLoops.CANDY_CANE.loop.getRectBound().latHi().degrees() > 10);
    assertTrue(TestLoops.SMALL_NE_CW.loop.getRectBound().isFull());
    // Permit small latitude error since we pad the loop bounds a bit, but longitude must be exact.
    double latError = 2E-15;
    assertBoundsEqual(
        TestLoops.ARCTIC_80.loop.getRectBound(),
        new S2LatLngRect(S2LatLng.fromDegrees(80, -180), S2LatLng.fromDegrees(90, 180)),
        S2LatLng.fromRadians(latError, 0));
    assertBoundsEqual(
        TestLoops.ANTARCTIC_80.loop.getRectBound(),
        new S2LatLngRect(S2LatLng.fromDegrees(-90, -180), S2LatLng.fromDegrees(-80, 180)),
        S2LatLng.fromRadians(latError, 0));

    TestLoops.ARCTIC_80.loop.invert();
    // The highest latitude of each edge is attained at its midpoint.
    S2Point mid =
        TestLoops.ARCTIC_80.loop.vertex(0).add(TestLoops.ARCTIC_80.loop.vertex(1)).mul(0.5);
    assertDoubleNear(
        TestLoops.ARCTIC_80.loop.getRectBound().latHi().radians(),
        new S2LatLng(mid).lat().radians());
    TestLoops.ARCTIC_80.loop.invert();

    assertTrue(TestLoops.SOUTH_HEMI.loop.getRectBound().lng().isFull());
    assertBoundsEqual(
        TestLoops.SOUTH_HEMI.loop.getRectBound().lat(), new R1Interval(-M_PI_2, 0), latError);
  }

  public void testFastConstructor() {
    List<S2Point> vertices = Lists.newArrayList();
    S2LatLngRect bound = parseVertices("-80:120, -80:0, -80:-120", vertices);
    S2Loop loop = S2Loop.newLoopWithTrustedDetails(vertices, false, bound);
    assertTrue(loop.isValid());
    assertEquals(0, loop.depth());
  }

  public void testAreaCentroid() {
    assertExactly(0.0, S2Loop.empty().getArea());
    assertExactly(4 * PI, S2Loop.full().getArea());
    assertNull(S2Loop.empty().getCentroid());
    assertNull(S2Loop.full().getCentroid());

    assertDoubleEquals(TestLoops.NORTH_HEMI.loop.getArea(), 2 * PI);
    assertDoubleEquals(TestLoops.EAST_HEMI.loop.getArea(), 2 * PI);

    // Construct spherical caps of random height, and approximate their boundary with closely
    // spaced vertices. Then check that the area and centroid are correct.

    for (int i = 0; i < 100; ++i) {
      // Choose a coordinate frame for the spherical cap.
      S2Point x = data.getRandomPoint();
      S2Point y = x.crossProd(data.getRandomPoint()).normalize();
      S2Point z = x.crossProd(y).normalize();

      // Given two points at latitude phi and whose longitudes differ by dtheta, the geodesic
      // between the two points has a maximum latitude of atan(tan(phi) / cos(dtheta/2)). This can
      // be derived by positioning the two points at (-dtheta/2, phi) and (dtheta/2, phi).
      //
      // We want to position the vertices close enough together so that their maximum distance from
      // the boundary of the spherical cap is kMaxDist. Thus we want:
      //   fabs(atan(tan(phi) / cos(dtheta/2)) - phi) <= kMaxDist.
      double kMaxDist = TestPlatform.S2_LOOP_TEST_MAX_ANGLE_DELTA;
      double height = 2 * data.nextDouble();
      double phi = asin(1 - height);
      double maxDtheta = 2 * acos(tan(abs(phi)) / tan(abs(phi) + kMaxDist));
      maxDtheta = min(PI, maxDtheta); // At least 3 vertices.

      List<S2Point> vertices = Lists.newArrayList();
      for (double theta = 0; theta < 2 * PI; theta += data.nextDouble() * maxDtheta) {

        S2Point xCosThetaCosPhi = x.mul(cos(theta) * cos(phi));
        S2Point ySinThetaCosPhi = y.mul(sin(theta) * cos(phi));
        S2Point zSinPhi = z.mul(sin(phi));

        S2Point sum = xCosThetaCosPhi.add(ySinThetaCosPhi).add(zSinPhi);

        vertices.add(sum);
      }

      S2Loop loop = new S2Loop(vertices);
      S2AreaCentroid areaCentroid = loop.getAreaAndCentroid();

      double area = loop.getArea();
      S2Point centroid = loop.getCentroid();
      double expectedArea = 2 * PI * height;
      assertTrue(areaCentroid.getArea() == area);
      assertTrue(centroid.equals(areaCentroid.getCentroid()));
      assertTrue(abs(area - expectedArea) <= 2 * PI * kMaxDist);

      // high probability
      assertTrue(abs(area - expectedArea) >= 0.01 * kMaxDist);

      S2Point expectedCentroid = z.mul(expectedArea * (1 - 0.5 * height));

      assertTrue(centroid.sub(expectedCentroid).norm() <= 2 * kMaxDist);
    }
  }

  public void testAreaConsistentWithTurningAngle() {
    // Check that the area computed using getArea() is consistent with the turning angle of the loop
    // computed using getTurnAngle(). According to the Gauss-Bonnet theorem, the area of the loop
    // should be equal to 2*Pi minus its turning angle.
    for (TestLoops test : TestLoops.values()) {
      S2Loop loop = test.loop;
      if (loop.isValid()) {
        double area = loop.getArea();
        double gaussArea = 2 * PI - loop.getTurningAngle();
        // TODO(user): The error bound below is much larger than it should be. Need to
        // improve the error minimization analysis in S2.area().
        assertTrue(
            "Failed loop: " + loop + "\nArea = " + area + ", Gauss Area = " + gaussArea,
            abs(area - gaussArea) <= 1e-9);
      }
    }
  }

  public void testGetAreaConsistentWithRobustCCW() {
    // Test that getArea() returns an area near 0 for degenerate loops that contain almost no
    // points, and an area near 4*Pi for degenerate loops that contain almost all points.
    final int maxVertices = 6;
    for (int i = 0; i < 50; i++) {
      int numVertices = 3 + data.random(maxVertices - 3 + 1);
      // Repeatedly choose N vertices that are exactly on the equator until we find some that form a
      // valid loop.
      S2Loop loop;
      do {
        List<S2Point> vertices = new ArrayList<>();
        for (int j = 0; j < numVertices; j++) {
          // We limit longitude to the range [0, 90] to ensure that the loop is degenerate (as
          // opposed to following the entire equator).
          vertices.add(S2LatLng.fromRadians(0, data.uniform(0, M_PI_2)).toPoint());
        }
        loop = new S2Loop(vertices);
      } while (!loop.isValid());
      boolean ccw = loop.isNormalized();
      // TODO(user): The error bound below is much larger than it should be. Need to improve
      // the error minimization analysis in S2.area().
      assertEquals("Failed loop " + i + ": " + loop, ccw ? 0 : 4 * PI, loop.getArea(), 1e-8);
      assertEquals(!ccw, loop.contains(S2Point.Z_POS));
    }
  }

  public void testGetAreaAccuracy() {
    // TODO(user): Test that GetArea() has an accuracy significantly better than 1e-15 on
    // loops whose area is small.
  }

  private static S2Loop rotate(S2Loop loop) {
    List<S2Point> vertices = Lists.newArrayList();
    for (int i = 1; i <= loop.numVertices(); ++i) {
      vertices.add(loop.vertex(i));
    }
    return new S2Loop(vertices);
  }

  public void testContains() {
    // Check the full and empty loops have the correct containment relationship with the special
    // "vertex" that defines them.
    assertFalse(S2Loop.empty().contains(S2Loop.EMPTY_VERTEX));
    assertTrue(S2Loop.full().contains(S2Loop.FULL_VERTEX));

    assertTrue(TestLoops.CANDY_CANE.loop.contains(S2LatLng.fromDegrees(5, 71).toPoint()));
    S2Loop northHemi = TestLoops.NORTH_HEMI.loop;
    S2Loop southHemi = TestLoops.SOUTH_HEMI.loop;
    S2Loop eastHemi = TestLoops.EAST_HEMI.loop;
    S2Loop westHemi = TestLoops.WEST_HEMI.loop;
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

    // This code checks each cell vertex is contained by exactly one of the adjacent cells.
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
        assertEquals(1, count);
      }
    }
  }

  public void testContainsMatchesRobustCrossing() {
    // This test demonstrates a former incompatibility between robustCrossing() and
    // contains(S2Point).  It constructs an S2Cell-based loop L and an edge E from Origin to a0 that
    // crosses exactly one edge of L. Yet previously, contains() returned false for both endpoints
    // of E.
    //
    // The reason for the bug was that the loop bound was sometimes too tight. The contains() code
    // for a0 bailed out early because a0 was found not to be inside the bound of L.

    // Start with a cell that ends up producing the problem.
    S2CellId cellId = S2CellId.fromPoint(new S2Point(1, 1, 1)).parent(21);

    S2Cell[] children = {new S2Cell(), new S2Cell(), new S2Cell(), new S2Cell()};
    new S2Cell(cellId).subdivide(children);

    S2Point[] points = new S2Point[4];
    for (int i = 0; i < 4; ++i) {
      // Note extra normalization. GetCenter() is already normalized. The test results will no
      // longer be inconsistent if the extra normalize() is removed.
      points[i] = children[i].getCenter().normalize();
    }

    S2Loop loop = new S2Loop(Arrays.asList(points));

    // Get a vertex from a grandchild cell.  Mathematically, a0 should be the same as points[0], but
    // rounding errors make it slightly different.
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
    // Construct a CCW polygon whose boundary is the union of the cell ids in the range
    // [begin, end). We add the edges one by one, removing any edges that are already present in the
    // opposite direction.
    Map<S2Point, Set<S2Point>> edges = Maps.newHashMap();
    for (S2CellId id = begin; !id.equals(end); id = id.next()) {
      S2Cell cell = new S2Cell(id);
      for (int k = 0; k < 4; ++k) {
        S2Point a = cell.getVertex(k);
        S2Point b = cell.getVertex((k + 1) & 3);
        if (edges.get(b) == null) {
          edges.put(b, Sets.<S2Point>newHashSet());
        }

        // If a is in b's set, remove it and remove b's set if it's empty. Otherwise, add b to a's
        // set.
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

    // The remaining edges form a single loop. We simply follow it starting at an arbitrary vertex
    // and build up a list of vertices.
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
   * Verifies the relationship between two loops A and B. "flags" is the set of RelationFlags that
   * apply. "shared_edge" means that the loops share at least one edge (possibly reversed.)
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
    } else if (!test(flags, (CONTAINS | CONTAINED | COVERS))) {
      checkOverlappingPair(a, b);
    }
    if (!sharedEdge && test(flags, (CONTAINS | CONTAINED | DISJOINT))) {
      assertEquals(a.contains(b), a.containsNested(b));
    }
    // A contains the boundary of B if either A contains B, or the two loops
    // contain each other's boundaries and there are no shared edges (since at
    // least one such edge must be reversed, and therefore is not considered to
    // be contained according to the rules of CompareBoundary).
    int comparison = 0;
    if (test(flags, CONTAINS) || (test(flags, COVERS) && !sharedEdge)) {
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
      System.err.println("args " + a + ", " + b);
      throw e;
    }
  }

  public void testLoopRelations() {
    // Check full and empty relationships with normal loops and each other.
    checkRelation(S2Loop.full(), S2Loop.full(), CONTAINS | CONTAINED | COVERS, true);
    checkRelation(S2Loop.full(), TestLoops.NORTH_HEMI.loop, CONTAINS | COVERS, false);
    checkRelation(S2Loop.full(), S2Loop.empty(), CONTAINS | DISJOINT | COVERS, false);
    checkRelation(TestLoops.NORTH_HEMI.loop, S2Loop.full(), CONTAINED | COVERS, false);
    checkRelation(TestLoops.NORTH_HEMI.loop, S2Loop.empty(), CONTAINS | DISJOINT, false);
    checkRelation(S2Loop.empty(), S2Loop.full(), CONTAINED | DISJOINT | COVERS, false);
    checkRelation(S2Loop.empty(), TestLoops.NORTH_HEMI.loop, CONTAINED | DISJOINT, false);
    checkRelation(S2Loop.empty(), S2Loop.empty(), CONTAINS | CONTAINED | DISJOINT, false);

    checkRelation(TestLoops.NORTH_HEMI.loop, TestLoops.NORTH_HEMI.loop, CONTAINS | CONTAINED, true);
    checkRelation(TestLoops.NORTH_HEMI.loop, TestLoops.SOUTH_HEMI.loop, DISJOINT | COVERS, true);
    checkRelation(TestLoops.NORTH_HEMI.loop, TestLoops.EAST_HEMI.loop, 0, false);
    checkRelation(TestLoops.NORTH_HEMI.loop, TestLoops.ARCTIC_80.loop, CONTAINS, false);
    checkRelation(TestLoops.NORTH_HEMI.loop, TestLoops.ANTARCTIC_80.loop, DISJOINT, false);
    checkRelation(TestLoops.NORTH_HEMI.loop, TestLoops.CANDY_CANE.loop, 0, false);

    // We can't compare north_hemi3 vs. north_hemi or south_hemi because the result depends on the
    // "simulation of simplicity" implementation details.

    checkRelation(
        TestLoops.NORTH_HEMI3.loop, TestLoops.NORTH_HEMI3.loop, CONTAINS | CONTAINED, true);
    checkRelation(TestLoops.NORTH_HEMI3.loop, TestLoops.EAST_HEMI.loop, 0, false);
    checkRelation(TestLoops.NORTH_HEMI3.loop, TestLoops.ARCTIC_80.loop, CONTAINS, false);
    checkRelation(TestLoops.NORTH_HEMI3.loop, TestLoops.ANTARCTIC_80.loop, DISJOINT, false);
    checkRelation(TestLoops.NORTH_HEMI3.loop, TestLoops.CANDY_CANE.loop, 0, false);

    checkRelation(TestLoops.SOUTH_HEMI.loop, TestLoops.NORTH_HEMI.loop, DISJOINT | COVERS, true);
    checkRelation(TestLoops.SOUTH_HEMI.loop, TestLoops.SOUTH_HEMI.loop, CONTAINS | CONTAINED, true);
    checkRelation(TestLoops.SOUTH_HEMI.loop, TestLoops.FAR_HEMI.loop, 0, false);
    checkRelation(TestLoops.SOUTH_HEMI.loop, TestLoops.ARCTIC_80.loop, DISJOINT, false);
    checkRelation(TestLoops.SOUTH_HEMI.loop, TestLoops.ANTARCTIC_80.loop, CONTAINS, false);
    checkRelation(TestLoops.SOUTH_HEMI.loop, TestLoops.CANDY_CANE.loop, 0, false);

    checkRelation(TestLoops.CANDY_CANE.loop, TestLoops.NORTH_HEMI.loop, 0, false);
    checkRelation(TestLoops.CANDY_CANE.loop, TestLoops.SOUTH_HEMI.loop, 0, false);
    checkRelation(TestLoops.CANDY_CANE.loop, TestLoops.ARCTIC_80.loop, DISJOINT, false);
    checkRelation(TestLoops.CANDY_CANE.loop, TestLoops.ANTARCTIC_80.loop, DISJOINT, false);
    checkRelation(TestLoops.CANDY_CANE.loop, TestLoops.CANDY_CANE.loop, CONTAINS | CONTAINED, true);

    checkRelation(TestLoops.NEAR_HEMI.loop, TestLoops.WEST_HEMI.loop, 0, false);

    checkRelation(TestLoops.SMALL_NE_CW.loop, TestLoops.SOUTH_HEMI.loop, CONTAINS, false);
    checkRelation(TestLoops.SMALL_NE_CW.loop, TestLoops.WEST_HEMI.loop, CONTAINS, false);

    checkRelation(TestLoops.SMALL_NE_CW.loop, TestLoops.NORTH_HEMI.loop, COVERS, false);
    checkRelation(TestLoops.SMALL_NE_CW.loop, TestLoops.EAST_HEMI.loop, COVERS, false);

    checkRelation(TestLoops.LOOP_A.loop, TestLoops.LOOP_A.loop, CONTAINS | CONTAINED, true);
    checkRelation(TestLoops.LOOP_A.loop, TestLoops.LOOP_B.loop, 0, false);
    checkRelation(TestLoops.LOOP_A.loop, TestLoops.A_INTERSECT_B.loop, CONTAINS, true);
    checkRelation(TestLoops.LOOP_A.loop, TestLoops.A_UNION_B.loop, CONTAINED, true);
    checkRelation(TestLoops.LOOP_A.loop, TestLoops.A_MINUS_B.loop, CONTAINS, true);
    checkRelation(TestLoops.LOOP_A.loop, TestLoops.B_MINUS_A.loop, DISJOINT, true);

    checkRelation(TestLoops.LOOP_B.loop, TestLoops.LOOP_A.loop, 0, false);
    checkRelation(TestLoops.LOOP_B.loop, TestLoops.LOOP_B.loop, CONTAINS | CONTAINED, true);
    checkRelation(TestLoops.LOOP_B.loop, TestLoops.A_INTERSECT_B.loop, CONTAINS, true);
    checkRelation(TestLoops.LOOP_B.loop, TestLoops.A_UNION_B.loop, CONTAINED, true);
    checkRelation(TestLoops.LOOP_B.loop, TestLoops.A_MINUS_B.loop, DISJOINT, true);
    checkRelation(TestLoops.LOOP_B.loop, TestLoops.B_MINUS_A.loop, CONTAINS, true);

    checkRelation(TestLoops.A_INTERSECT_B.loop, TestLoops.LOOP_A.loop, CONTAINED, true);
    checkRelation(TestLoops.A_INTERSECT_B.loop, TestLoops.LOOP_B.loop, CONTAINED, true);
    checkRelation(
        TestLoops.A_INTERSECT_B.loop, TestLoops.A_INTERSECT_B.loop, CONTAINS | CONTAINED, true);
    checkRelation(TestLoops.A_INTERSECT_B.loop, TestLoops.A_UNION_B.loop, CONTAINED, false);
    checkRelation(TestLoops.A_INTERSECT_B.loop, TestLoops.A_MINUS_B.loop, DISJOINT, true);
    checkRelation(TestLoops.A_INTERSECT_B.loop, TestLoops.B_MINUS_A.loop, DISJOINT, true);

    checkRelation(TestLoops.A_UNION_B.loop, TestLoops.LOOP_A.loop, CONTAINS, true);
    checkRelation(TestLoops.A_UNION_B.loop, TestLoops.LOOP_B.loop, CONTAINS, true);
    checkRelation(TestLoops.A_UNION_B.loop, TestLoops.A_INTERSECT_B.loop, CONTAINS, false);
    checkRelation(TestLoops.A_UNION_B.loop, TestLoops.A_UNION_B.loop, CONTAINS | CONTAINED, true);
    checkRelation(TestLoops.A_UNION_B.loop, TestLoops.A_MINUS_B.loop, CONTAINS, true);
    checkRelation(TestLoops.A_UNION_B.loop, TestLoops.B_MINUS_A.loop, CONTAINS, true);

    checkRelation(TestLoops.A_MINUS_B.loop, TestLoops.LOOP_A.loop, CONTAINED, true);
    checkRelation(TestLoops.A_MINUS_B.loop, TestLoops.LOOP_B.loop, DISJOINT, true);
    checkRelation(TestLoops.A_MINUS_B.loop, TestLoops.A_INTERSECT_B.loop, DISJOINT, true);
    checkRelation(TestLoops.A_MINUS_B.loop, TestLoops.A_UNION_B.loop, CONTAINED, true);
    checkRelation(TestLoops.A_MINUS_B.loop, TestLoops.A_MINUS_B.loop, CONTAINS | CONTAINED, true);
    checkRelation(TestLoops.A_MINUS_B.loop, TestLoops.B_MINUS_A.loop, DISJOINT, false);

    checkRelation(TestLoops.B_MINUS_A.loop, TestLoops.LOOP_A.loop, DISJOINT, true);
    checkRelation(TestLoops.B_MINUS_A.loop, TestLoops.LOOP_B.loop, CONTAINED, true);
    checkRelation(TestLoops.B_MINUS_A.loop, TestLoops.A_INTERSECT_B.loop, DISJOINT, true);
    checkRelation(TestLoops.B_MINUS_A.loop, TestLoops.A_UNION_B.loop, CONTAINED, true);
    checkRelation(TestLoops.B_MINUS_A.loop, TestLoops.A_MINUS_B.loop, DISJOINT, false);
    checkRelation(TestLoops.B_MINUS_A.loop, TestLoops.B_MINUS_A.loop, CONTAINS | CONTAINED, true);
  }

  /**
   * Make sure the relations are correct if the loop crossing happens on two ends of a shared
   * boundary segment.
   */
  public void testLoopRelationsWhenSameExceptPiecesStickingOutAndIn() {
    checkRelation(TestLoops.LOOP_A.loop, TestLoops.LOOP_C.loop, 0, true);
    checkRelation(TestLoops.LOOP_C.loop, TestLoops.LOOP_A.loop, 0, true);
    checkRelation(TestLoops.LOOP_A.loop, TestLoops.LOOP_D.loop, CONTAINED, true);
    checkRelation(TestLoops.LOOP_D.loop, TestLoops.LOOP_A.loop, CONTAINS, true);
    checkRelation(TestLoops.LOOP_E.loop, TestLoops.LOOP_F.loop, DISJOINT, true);
    checkRelation(TestLoops.LOOP_E.loop, TestLoops.LOOP_G.loop, CONTAINS, true);
    checkRelation(TestLoops.LOOP_E.loop, TestLoops.LOOP_H.loop, 0, true);
    checkRelation(TestLoops.LOOP_E.loop, TestLoops.LOOP_I.loop, 0, false);
    checkRelation(TestLoops.LOOP_F.loop, TestLoops.LOOP_G.loop, DISJOINT, true);
    checkRelation(TestLoops.LOOP_F.loop, TestLoops.LOOP_H.loop, 0, true);
    checkRelation(TestLoops.LOOP_F.loop, TestLoops.LOOP_I.loop, 0, false);
    checkRelation(TestLoops.LOOP_G.loop, TestLoops.LOOP_H.loop, CONTAINED, true);
    checkRelation(TestLoops.LOOP_H.loop, TestLoops.LOOP_G.loop, CONTAINS, true);
    checkRelation(TestLoops.LOOP_G.loop, TestLoops.LOOP_I.loop, DISJOINT, true);
    checkRelation(TestLoops.LOOP_H.loop, TestLoops.LOOP_I.loop, CONTAINS, true);
  }

  public void testLoopRelations2() {
    // Construct polygons consisting of a sequence of adjacent cell ids at some fixed level.
    // Comparing two polygons at the same level ensures that there are no T-vertices.
    for (int iter = 0; iter < TestPlatform.S2_LOOP_TEST_LOOP_RELATIONS2_ITERATIONS; ++iter) {
      long num = data.nextLong();
      S2CellId begin = new S2CellId(num | 1);
      if (!begin.isValid()) {
        continue;
      }
      begin = begin.parent(data.random(S2CellId.MAX_LEVEL));
      S2CellId aBegin = advance(begin, data.skewed(6));
      S2CellId aEnd = advance(aBegin, data.skewed(6) + 1);
      S2CellId bBegin = advance(begin, data.skewed(6));
      S2CellId bEnd = advance(bBegin, data.skewed(6) + 1);
      if (!aEnd.isValid() || !bEnd.isValid()) {
        continue;
      }

      S2Loop a = makeCellLoop(aBegin, aEnd);
      S2Loop b = makeCellLoop(bBegin, bEnd);
      boolean contained = (aBegin.lessOrEquals(bBegin) && bEnd.lessOrEquals(aEnd));
      boolean intersects = (aBegin.lessThan(bEnd) && bBegin.lessThan(aEnd));
      String message = "Failed with " + a.numVertices() + " vs. " + b.numVertices();
      assertEquals(message, contained, a.contains(b));
      assertEquals(message, intersects, a.intersects(b));
    }
  }

  /**
   * This unit test helped us find an error in {@link S2EdgeQuery#splitBound} where {@code
   * edgeBound} was being used, instead of copies of {@code edgeBound}. This test fails at iteration
   * 21, where {@code numVertices} = 512, if the error is reintroduced.
   */
  public void testCompareBoundaryCrossingTestLoops() {
    for (int numVertices = 8; numVertices <= 512; numVertices *= 8) {
      for (int i = 0; i < 100; ++i) {
        List<S2Loop> loops = data.makeCrossingLoopPairDefault(numVertices,
            data.getRandomPoint(), data.getRandomPoint());
        S2Loop a = loops.get(0);
        S2Loop b = loops.get(1);
        assertEquals(0, a.compareBoundary(b));
      }
    }
  }

  public void testBoundsForLoopContainment() {
    // To reliably test whether one loop contains another, the bounds of the outer loop are expanded
    // slightly.  This test constructs examples where this expansion is necessary and verifies that
    // it is sufficient.
    for (int iter = 0; iter < 1000; ++iter) {
      // We construct a triangle ABC such that A,B,C are nearly colinear, B is the point of maximum
      // latitude, and the edge AC passes very slightly below B (i.e., ABC is CCW).
      S2Point b = data.getRandomPoint().add(S2Point.Z_POS).normalize();
      S2Point v = b.crossProd(S2Point.Z_POS).normalize();
      S2Point a = S2EdgeUtil.interpolate(data.nextDouble(), v.neg(), b);
      S2Point c = S2EdgeUtil.interpolate(data.nextDouble(), b, v);
      if (S2Predicates.sign(a, b, c) < 0) {
        --iter;
        continue;
      }
      // Now construct another point D directly below B, and create two loops ABCD and ACD.
      S2Point d = new S2Point(b.getX(), b.getY(), 0).normalize();
      List<S2Point> vertices = ImmutableList.of(c, d, a, b); // Reordered for convenience
      S2Loop outer = new S2Loop(vertices);
      S2Loop inner = new S2Loop(vertices.subList(0, 3));
      // Now because the bounds calculation is less accurate when the maximum is attained along an
      // edge (rather than at a vertex), sometimes the inner loop will have a *larger* bounding box
      // than the outer loop.  We look only for those cases.
      if (outer.getRectBound().contains(inner.getRectBound())) {
        --iter;
        continue;
      }
      assertTrue(outer.contains(inner));
    }
  }

  private static void checkNear(String aStr, String bStr, double maxError, boolean expected) {
    S2Loop a = makeLoop(aStr);
    S2Loop b = makeLoop(bStr);
    assertEquals(expected, a.boundaryNear(b, maxError));
    assertEquals(expected, b.boundaryNear(a, maxError));
  }

  public void testBoundaryNear() {
    double degree = S1Angle.degrees(1).radians();

    checkNear("0:0, 0:10, 5:5", "0:0.1, -0.1:9.9, 5:5.2", 0.5 * degree, true);
    checkNear("0:0, 0:3, 0:7, 0:10, 3:7, 5:5", "0:0, 0:10, 2:8, 5:5, 4:4, 3:3, 1:1", 1e-3, true);

    // All vertices close to some edge, but not equivalent.
    checkNear("0:0, 0:2, 2:2, 2:0", "0:0, 1.9999:1, 0:2, 2:2, 2:0", 0.5 * degree, false);

    // Two triangles that backtrack a bit on different edges. A simple greedy matching algorithm
    // would fail on this example.
    String t1 =
        "0.1:0, 0.1:1, 0.1:2, 0.1:3, 0.1:4, 1:4, 2:4, 3:4, "
            + "2:4.1, 1:4.1, 2:4.2, 3:4.2, 4:4.2, 5:4.2";
    String t2 =
        "0:0, 0:1, 0:2, 0:3, 0.1:2, 0.1:1, 0.2:2, 0.2:3, " + "0.2:4, 1:4.1, 2:4, 3:4, 4:4, 5:4";
    checkNear(t1, t2, 1.5 * degree, true);
    checkNear(t1, t2, 0.5 * degree, false);
  }

  private static void checkIdentical(S2Loop a, S2Loop b) {
    assertEquals(a.numVertices(), b.numVertices());
    for (int i = 0; i < a.numVertices(); i++) {
      assertEquals(a.vertex(i), b.vertex(i));
    }
    assertEquals(a.depth(), b.depth());
    assertEquals(a.isEmpty(), b.isEmpty());
    assertEquals(a.isFull(), b.isFull());
    assertEquals(a.isNormalized(), b.isNormalized());
    assertEquals(a.contains(S2Point.ORIGIN), b.contains(S2Point.ORIGIN));
    assertEquals(a.getRectBound(), b.getRectBound());
  }

  private static byte[] encodeCompressed(S2Loop loop, int level) throws IOException {
    ByteArrayOutputStream outputStream = new ByteArrayOutputStream();
    LittleEndianOutput encoder = new LittleEndianOutput(outputStream);
    loop.encodeCompressed(level, encoder);
    encoder.close();
    return outputStream.toByteArray();
  }

  private static S2Loop decodeCompressed(int level, byte[] encoded) throws IOException {
    LittleEndianInput decoder = new LittleEndianInput(new ByteArrayInputStream(encoded));
    S2Loop loop = S2Loop.decodeCompressed(level, decoder);
    assertNotNull(loop);
    return loop;
  }

  private static void encodeDecodeCompressed(S2Loop loop, int level) throws IOException {
    byte[] encoded = encodeCompressed(loop, level);
    S2Loop decodedLoop = decodeCompressed(level, encoded);
    checkIdentical(loop, decodedLoop);
  }

  public void testEncodeDecodeCompressed() throws IOException {
    // Empty loop.
    encodeDecodeCompressed(S2Loop.empty(), S2CellId.MAX_LEVEL);

    // Full loop.
    encodeDecodeCompressed(S2Loop.full(), S2CellId.MAX_LEVEL);

    // Unsnapped loop.
    S2Loop l = makeLoop("30:20, 40:20, 39:43, 33:35");
    l.setDepth(3);
    encodeDecodeCompressed(l, S2CellId.MAX_LEVEL);

    // Snapped loop.
    S2Loop snappedLoopAClone = new S2Loop(TestLoops.SNAPPED_LOOP_A.loop);
    snappedLoopAClone.setDepth(3);
    encodeDecodeCompressed(snappedLoopAClone, S2CellId.MAX_LEVEL);
  }

  public void testFourVertexCompressedLoopSize() throws IOException {
    byte[] encoded = encodeCompressed(TestLoops.SNAPPED_LOOP_A.loop, S2CellId.MAX_LEVEL);

    // 1 byte for numVertices
    // 1 byte for originInside and boolean indicating we did not encode the bound
    // 1 byte for depth
    // Vertices:
    // 1 byte for faces
    // 8 bytes for each vertex.
    // 1 byte indicating that there is no unsnapped vertex.
    assertEquals(37, encoded.length);
  }

  static void checkEmptyFullSnapped(S2Loop loop, int level) {
    assertTrue(loop.isEmptyOrFull());
    S2CellId cellid = S2CellId.fromPoint(loop.vertex(0)).parent(level);
    S2Loop loop2 = new S2Loop(ImmutableList.of(cellid.toPoint()));
    assertTrue(loop.boundaryEquals(loop2));
    assertTrue(loop.boundaryApproxEquals(loop2));
    assertTrue(loop.boundaryNear(loop2));
  }

  // Test converting the empty/full loops to S2LatLng representations. (We don't bother testing
  // E5/E6/E7 because that test is less demanding.)
  static void checkEmptyFullLatLng(S2Loop loop) {
    assertTrue(loop.isEmptyOrFull());
    S2Loop loop2 = new S2Loop(ImmutableList.of(new S2LatLng(loop.vertex(0)).toPoint()));
    assertTrue(loop.boundaryEquals(loop2));
    assertTrue(loop.boundaryApproxEquals(loop2));
    assertTrue(loop.boundaryNear(loop2));
  }

  static void checkEmptyFullConversions(S2Loop loop) {
    checkEmptyFullSnapped(loop, S2CellId.MAX_LEVEL);
    checkEmptyFullSnapped(loop, 1); // Worst case for approximation
    checkEmptyFullSnapped(loop, 0);
    checkEmptyFullLatLng(loop);
  }

  public void testEmptyFullLossyConversions() {
    // Verify that the empty and full loops can be encoded lossily.
    checkEmptyFullConversions(S2Loop.empty());
    checkEmptyFullConversions(S2Loop.full());
  }

  /** Tests that nearly colinear points pass S2Loop.isValid() */
  public void testRoundingError() {
    S2Point a = new S2Point(-0.9190364081111774, 0.17231932652084575, 0.35451111445694833);
    S2Point b = new S2Point(-0.92130667053206, 0.17274500072476123, 0.3483578383756171);
    S2Point c = new S2Point(-0.9257244057938284, 0.17357332608634282, 0.3360158106235289);
    S2Point d = new S2Point(-0.9278712595449962, 0.17397586116468677, 0.32982923679138537);

    assertTrue(S2Loop.isValid(Lists.newArrayList(a, b, c, d)));
  }

  /** Tests {@link S2Loop#isValid()}. */
  public void testIsValid() {
    assertTrue(TestLoops.LOOP_A.loop.isValid());
    assertTrue(TestLoops.LOOP_B.loop.isValid());
    assertFalse(TestLoops.BOWTIE.loop.isValid());
  }

  /** Tests {@link S2Loop#compareTo(S2Loop)}. */
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

    // A square with (lat,lng) vertices (0,1), (1,1), (1,2) and (0,2).
    // Tests the case where the shortest distance is along a normal to an edge, onto a vertex
    S2Loop s1 = makeLoop("0:1, 1:1, 1:2, 0:2");

    // A square with (lat,lng) vertices (-1,1), (1,1), (1,2) and (-1,2).
    // Tests the case where the shortest distance is along a normal to an edge, not onto a vertex
    S2Loop s2 = makeLoop("-1:1, 1:1, 1:2, -1:2");

    // A diamond with (lat,lng) vertices (1,0), (2,1), (3,0) and (2,-1).
    // Test the case where the shortest distance is NOT along a normal to an edge
    S2Loop s3 = makeLoop("1:0, 2:1, 3:0, 2:-1");

    // All the vertices should be distance 0
    for (int i = 0; i < s1.numVertices(); i++) {
      assertEquals(0d, s1.getDistance(s1.vertex(i)).radians(), epsilon);
    }

    // A point on one of the edges should be distance 0
    assertEquals(0d, s1.getDistance(S2LatLng.fromDegrees(0.5, 1).toPoint()).radians(), epsilon);

    // In all three cases, the closest point to the origin is (0,1), which is at a distance of 1
    // degree.
    // Note: all of these are intentionally distances measured along the equator, since that makes
    // the math significantly simpler. Otherwise, the distance wouldn't actually be 1 degree.
    S2Point origin = S2LatLng.fromDegrees(0, 0).toPoint();
    assertEquals(1d, s1.getDistance(origin).degrees(), epsilon);
    assertEquals(1d, s2.getDistance(origin).degrees(), epsilon);
    assertEquals(1d, s3.getDistance(origin).degrees(), epsilon);
  }

  /**
   * Check that the turning angle is *identical* when the vertex order is rotated, and that the sign
   * is inverted when the vertices are reversed.
   */
  private static void checkTurningAngleInvariants(S2Loop loop) {
    double expected = loop.getTurningAngle();
    S2Loop loopCopy = new S2Loop(loop);
    for (int i = 0; i < loop.numVertices(); ++i) {
      loopCopy.invert();
      assertEquals(-expected, loopCopy.getTurningAngle(), 0.0);
      loopCopy.invert();
      loopCopy = rotate(loopCopy);
      assertEquals(expected, loopCopy.getTurningAngle(), 0.0);
    }
  }

  public void testTurningAngle() {
    assertEquals(2 * PI, S2Loop.empty().getTurningAngle(), 0.0);
    assertEquals(-2 * PI, S2Loop.full().getTurningAngle(), 0.0);

    assertDoubleNear(0, TestLoops.NORTH_HEMI3.loop.getTurningAngle(), 1e-15);
    checkTurningAngleInvariants(TestLoops.NORTH_HEMI3.loop);

    assertDoubleNear(0, TestLoops.WEST_HEMI.loop.getTurningAngle(), 1e-15);
    checkTurningAngleInvariants(TestLoops.WEST_HEMI.loop);

    // We don't have an easy way to estimate the turning angle of this loop, but we can still check
    // that the expected invariants hold.
    checkTurningAngleInvariants(TestLoops.CANDY_CANE.loop);

    assertEquals(2 * PI, TestLoops.LINE_TRIANGLE.loop.getTurningAngle(), 1e-14);
    checkTurningAngleInvariants(TestLoops.LINE_TRIANGLE.loop);

    assertEquals(2 * PI, TestLoops.SKINNY_CHEVRON.loop.getTurningAngle(), 1e-14);
    checkTurningAngleInvariants(TestLoops.SKINNY_CHEVRON.loop);

    // Builds a narrow spiral loop starting at the north pole.  This is designed to test that the
    // error in getTurningAngle is linear in the number of vertices even when the partial sum of the
    // turning angles gets very large. The spiral consists of two "arms" defining opposite sides of
    // the loop.
    int armVertexCount = 10000;
    double spiralRadius = 0.01;
    S2Point[] vertices = new S2Point[2 * armVertexCount];
    vertices[armVertexCount] = S2Point.Z_POS;
    for (int i = 0; i < armVertexCount; i++) {
      double angle = (2 * PI / 3) * i;
      double x = cos(angle);
      double y = sin(angle);
      double r1 = i * spiralRadius / armVertexCount;
      double r2 = (i + 1.5) * spiralRadius / armVertexCount;
      vertices[armVertexCount - i - 1] = new S2Point(r1 * x, r1 * y, 1).normalize();
      vertices[armVertexCount + i] = new S2Point(r2 * x, r2 * y, 1).normalize();
    }

    // Check that getTurningAngle() is consistent with getArea() to within the error bound of the
    // former.  We actually use a tiny fraction of the worst-case error bound, since the worst case
    // only happens when all the roundoff errors happen in the same direction and this test is not
    // designed to achieve that. The error in getArea() can be ignored for the purposes of this test
    // since it is generally much smaller.
    S2Loop spiral = new S2Loop(Arrays.asList(vertices));
    assertEquals(
        2 * PI - spiral.getArea(),
        spiral.getTurningAngle(),
        0.01 * S2.getTurningAngleMaxError(spiral.numVertices()));
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
    checkNormalizeAndContains(TestLoops.LINE_TRIANGLE.loop);
    checkNormalizeAndContains(TestLoops.SKINNY_CHEVRON.loop);
  }

  public void testIsOriginInside() {
    assertEquals(
        TestLoops.EAST_HEMI.loop.contains(S2.origin()), TestLoops.EAST_HEMI.loop.isOriginInside());
    assertEquals(
        TestLoops.ARCTIC_80.loop.contains(S2.origin()), TestLoops.ARCTIC_80.loop.isOriginInside());
  }

  // Some literals here cannot be precisely represented as doubles, but they're the same as the
  // literals used in the C++ tests.
  @SuppressWarnings("FloatingPointLiteralPrecision")
  public void testMakeRegularLoop() {
    S2Point center = S2LatLng.fromDegrees(80, 135).toPoint();
    S1Angle radius = S1Angle.degrees(20);
    S2Loop loop = S2Loop.makeRegularLoop(center, radius, 4);

    assertEquals(4, loop.numVertices());
    S2Point p0 = loop.vertex(0);
    S2Point p1 = loop.vertex(1);
    S2Point p2 = loop.vertex(2);
    S2Point p3 = loop.vertex(3);

    // Make sure that the radius is correct.
    assertDoubleEquals(20.0, toDegrees(center.angle(p0)));
    assertDoubleEquals(20.0, toDegrees(center.angle(p1)));
    assertDoubleEquals(20.0, toDegrees(center.angle(p2)));
    assertDoubleEquals(20.0, toDegrees(center.angle(p3)));

    // Make sure that all the angles of the polygon are the same.
    assertDoubleEquals(M_PI_2, p1.sub(p0).angle(p3.sub(p0)));
    assertDoubleEquals(M_PI_2, p2.sub(p1).angle(p0.sub(p1)));
    assertDoubleEquals(M_PI_2, p3.sub(p2).angle(p1.sub(p2)));
    assertDoubleEquals(M_PI_2, p0.sub(p3).angle(p2.sub(p3)));

    // Make sure that all the edges of the polygon have the same length.
    assertDoubleEquals(27.990890717782829, toDegrees(p0.angle(p1)));
    assertDoubleEquals(27.990890717782829, toDegrees(p1.angle(p2)));
    assertDoubleEquals(27.990890717782829, toDegrees(p2.angle(p3)));
    assertDoubleEquals(27.990890717782829, toDegrees(p3.angle(p0)));

    // Check actual coordinates. This may change if we switch the algorithm intentionally. These
    // values are exactly the same as the C++ versions, but with slight tolerances in a couple of
    // places where we lose precision at different points due to using strictfp.
    assertExactly(62.162880741097204, new S2LatLng(p0).lat().degrees());
    assertExactly(103.11051028343408, new S2LatLng(p0).lng().degrees());
    assertExactly(61.955157772928345, new S2LatLng(p1).lat().degrees());
    assertExactly(165.25681963683536, new S2LatLng(p1).lng().degrees());
    assertEquals(75.139812547718478, new S2LatLng(p2).lat().degrees(), 2e-14);
    assertExactly(-119.13042521187423, new S2LatLng(p2).lng().degrees());
    assertEquals(75.524190079054392, new S2LatLng(p3).lat().degrees(), 2e-14);
    assertExactly(26.392175948257943, new S2LatLng(p3).lng().degrees());
  }

  public void testSimplify() {
    // Ensure a redundant point in the middle is not discarded when outside the tolerance.
    S2Loop loop = makeLoop("0:0, 1e-10:5, 0:10, 5:5");
    assertEquals(loop, loop.simplify(S1Angle.radians(1e-12), null));
    // Ensure the point is simplified away when inside the tolerance.
    assertEquals(makeLoop("0:0, 0:10, 5:5"), loop.simplify(S1Angle.radians(1e-8), null));
    // Unless it's in the vertex merge filter.
    assertEquals(
        loop,
        loop.simplify(
            S1Angle.radians(1e-8),
            new Predicate<S2Point>() {
              @Override
              public boolean apply(S2Point vertex) {
                return makePoint("1e-10:5").equalsPoint(vertex);
              }
            }));
  }

  public void testEncodeDecode() throws Exception {
    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    S2Loop loop = makeLoop("1:2, 3:4, 5:6");
    loop.encode(new LittleEndianOutput(baos));

    LittleEndianInput decoder = new LittleEndianInput(new ByteArrayInputStream(baos.toByteArray()));
    S2Loop decodedLoop = S2Loop.decode(decoder);
    assertEquals(loop, decodedLoop);
    assertEquals(decodedLoop.depth(), loop.depth());
    assertFalse(loop.isOriginInside());
    assertEquals(loop.getRectBound(), decodedLoop.getRectBound());
  }

  public void testDecodeEmptyLoopLossless() throws IOException {
    String encodedBytesHexString =
        "010100000000000000000000000000000000000000000000000000F03F0000000000010000"
            + "00000000F03F0000000000000000182D4454FB210940182D4454FB2109C0";
    ByteArrayInputStream bais =
        new ByteArrayInputStream(BaseEncoding.base16().decode(encodedBytesHexString));
    LittleEndianInput decoder = new LittleEndianInput(bais);
    assertTrue(S2Loop.decode(decoder).isEmpty());
  }

  public void testDecodeFullLoopLossless() throws IOException {
    String encodedBytesHexString =
        "010100000000000000000000000000000000000000000000000000F0BF010000000001182D"
            + "4454FB21F9BF182D4454FB21F93F182D4454FB2109C0182D4454FB210940";
    ByteArrayInputStream bais =
        new ByteArrayInputStream(BaseEncoding.base16().decode(encodedBytesHexString));
    LittleEndianInput decoder = new LittleEndianInput(bais);
    assertTrue(S2Loop.decode(decoder).isFull());
  }

  public void testDecodeWellknownLoopLossless() throws IOException {
    String encodedBytesHexString =
        "0108000000D44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA1BFB4825F3C81FDEF3F"
            + "27DCF7C958DE913F1EDD892B0BDF91BFB4825F3C81FDEF3F27DCF7C958DE913F1EDD892B0B"
            + "DF913FD44A8442C3F9EF3F7EDA2AB341DC913F27DCF7C958DEA13FD44A8442C3F9EF3F7EDA"
            + "2AB341DC91BF27DCF7C958DEA13FB4825F3C81FDEF3F27DCF7C958DE91BF1EDD892B0BDF91"
            + "3FB4825F3C81FDEF3F27DCF7C958DE91BF1EDD892B0BDF91BFD44A8442C3F9EF3F7EDA2AB3"
            + "41DC91BF27DCF7C958DEA1BF0000000000013EFC10E8F8DFA1BF3EFC10E8F8DFA13F389D52"
            + "A246DF91BF389D52A246DF913F";
    ByteArrayInputStream bais =
        new ByteArrayInputStream(BaseEncoding.base16().decode(encodedBytesHexString));
    S2Loop result = S2Loop.decode(new LittleEndianInput(bais));
    assertEquals(makeLoop("-2:1, -1:1, 1:1, 2:1, 2:-1, 1:-1, -1:-1, -2:-1"), result);
  }

  public void testEmptyShape() {
    S2Shape shape = S2Loop.empty();
    assertTrue(shape.hasInterior());
    assertFalse(shape.containsOrigin());
    assertEquals(0, shape.numEdges());
    assertEquals(0, shape.numChains());
    assertEquals(2, shape.dimension());
  }

  public void testFullShape() {
    S2Shape shape = S2Loop.full();
    assertTrue(shape.hasInterior());
    assertTrue(shape.containsOrigin());
    assertEquals(0, shape.numEdges());
    assertEquals(1, shape.numChains());
    checkFirstNChainStarts(shape, 0);
    checkFirstNChainLengths(shape, 0);
    assertEquals(2, shape.dimension());
  }

  public void testThreeEdgeShape() {
    assertTrue(TestLoops.SMALL_NE_CW.loop.hasInterior());
    assertTrue(TestLoops.SMALL_NE_CW.loop.containsOrigin());
    assertEquals(3, TestLoops.SMALL_NE_CW.loop.numEdges());
    checkFirstNEdges(TestLoops.SMALL_NE_CW.loop, "35:20|45:20, 45:20|40:25, 40:25|35:20");
    assertEquals(1, TestLoops.SMALL_NE_CW.loop.numChains());
    checkFirstNChainStarts(TestLoops.SMALL_NE_CW.loop, 0);
    checkFirstNChainLengths(TestLoops.SMALL_NE_CW.loop, 3);
    checkFirstNChainEdges(TestLoops.SMALL_NE_CW.loop, 0, "35:20|45:20, 45:20|40:25, 40:25|35:20");
    assertEquals(2, TestLoops.SMALL_NE_CW.loop.dimension());
  }

  public void testFourEdgeShape() {
    assertTrue(TestLoops.EAST_HEMI.loop.hasInterior());
    assertTrue(TestLoops.EAST_HEMI.loop.containsOrigin());
    assertEquals(4, TestLoops.EAST_HEMI.loop.numEdges());
    checkFirstNEdges(TestLoops.EAST_HEMI.loop, "90:0|0:0, 0:0|-90:0, -90:0|0:-180, 0:-180|90:0");
    assertEquals(1, TestLoops.EAST_HEMI.loop.numChains());
    checkFirstNChainStarts(TestLoops.EAST_HEMI.loop, 0);
    checkFirstNChainLengths(TestLoops.EAST_HEMI.loop, 4);
    checkFirstNChainEdges(
        TestLoops.EAST_HEMI.loop, 0, "90:0|0:0, 0:0|-90:0, -90:0|0:-180, 0:-180|90:0");
    assertEquals(2, TestLoops.EAST_HEMI.loop.dimension());
  }

  /** This function is useful for debugging. */
  @SuppressWarnings("unused")
  private void dumpCrossings(S2Loop loop) {

    System.out.println("Ortho(v1): " + S2.ortho(loop.vertex(1)));
    System.out.println("Contains(kOrigin): " + loop.contains(S2.origin()) + "\n");
    for (int i = 1; i <= loop.numVertices(); ++i) {
      S2Point a = S2.ortho(loop.vertex(i));
      S2Point b = loop.vertex(i - 1);
      S2Point c = loop.vertex(i + 1);
      S2Point o = loop.vertex(i);
      Platform.printf(
          System.out,
          "Vertex %d: [%.17g, %.17g, %.17g], " + "%d%dR=%d, %d%d%d=%d, R%d%d=%d, inside: %b\n",
          i,
          loop.vertex(i).x,
          loop.vertex(i).y,
          loop.vertex(i).z,
          i - 1,
          i,
          S2Predicates.sign(b, o, a),
          i + 1,
          i,
          i - 1,
          S2Predicates.sign(c, o, b),
          i,
          i + 1,
          S2Predicates.sign(a, o, c),
          S2Predicates.orderedCCW(a, b, c, o));
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
      Platform.printf(
          System.out,
          "%d1R=%d, M1%d=%d, R1M=%d, crosses: %b\n",
          i,
          S2Predicates.sign(b, o, a),
          i,
          S2Predicates.sign(c, o, b),
          S2Predicates.sign(a, o, c),
          S2EdgeUtil.edgeOrVertexCrossing(c, o, b, a));
    }
  }
}
