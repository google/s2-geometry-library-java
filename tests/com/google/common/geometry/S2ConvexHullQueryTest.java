/*
 * Copyright 2016 Google Inc.
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
import com.google.common.base.Preconditions;
import java.util.ArrayList;

@GwtCompatible
public final strictfp class S2ConvexHullQueryTest extends GeometryTestCase {
  public void testNoPoints() {
    S2ConvexHullQuery query = new S2ConvexHullQuery();
    assertTrue(query.getConvexHull().isEmpty());
  }

  private static boolean loopHasVertex(S2Loop loop, S2Point p) {
    for (int i = 0; i < loop.numVertices(); ++i) {
      if (loop.vertex(i).equals(p)) {
        return true;
      }
    }
    return false;
  }

  public void testOnePoint() {
    S2ConvexHullQuery query = new S2ConvexHullQuery();
    S2Point p = new S2Point(0, 0, 1);
    query.addPoint(p);
    S2Loop result = query.getConvexHull();
    assertEquals(3, result.numVertices());
    assertTrue(result.isNormalized());
    assertTrue(loopHasVertex(result, p));
    // Add some duplicate points and check that the result is the same.
    query.addPoint(p);
    query.addPoint(p);
    S2Loop result2 = query.getConvexHull();
    assertTrue(result2.equals(result));
  }

  public void testTwoPoints() {
    S2ConvexHullQuery query = new S2ConvexHullQuery();
    S2Point p = new S2Point(0, 0, 1);
    S2Point q = new S2Point(0, 1, 0);
    query.addPoint(p);
    query.addPoint(q);
    S2Loop result = query.getConvexHull();
    assertEquals(3, result.numVertices());
    assertTrue(result.isNormalized());
    assertTrue(loopHasVertex(result, p));
    assertTrue(loopHasVertex(result, q));
    // Add some duplicate points and check that the result is the same.
    query.addPoint(q);
    query.addPoint(p);
    query.addPoint(p);
    S2Loop result2 = query.getConvexHull();
    assertTrue(result2.equals(result));
  }

  public void testTwoAntipodalPoints() {
    S2ConvexHullQuery query = new S2ConvexHullQuery();
    query.addPoint(S2Point.Z_POS);
    query.addPoint(S2Point.Z_NEG);
    S2Loop result = query.getConvexHull();
    assertTrue(result.isFull());
  }

  public void testEmptyLoop() {
    S2ConvexHullQuery query = new S2ConvexHullQuery();
    query.addLoop(S2Loop.empty());
    assertTrue(query.getConvexHull().isEmpty());
  }

  public void testFullLoop() {
    S2ConvexHullQuery query = new S2ConvexHullQuery();
    S2Loop full = S2Loop.full();
    query.addLoop(full);
    assertTrue(query.getConvexHull().isFull());
  }

  public void testEmptyPolygon() {
    S2ConvexHullQuery query = new S2ConvexHullQuery();
    S2Polygon empty = new S2Polygon(new ArrayList<S2Loop>());
    query.addPolygon(empty);
    assertTrue(query.getConvexHull().isEmpty());
  }

  public void testNonConvexPoints() {
    // Generate a point set such that the only convex region containing them is
    // the entire sphere.  In other words, you can generate any point on the
    // sphere by repeatedly linearly interpolating between the points.  (The
    // four points of a tetrahedron would also work, but this is easier.)
    S2ConvexHullQuery query = new S2ConvexHullQuery();
    for (int face = 0; face < 6; ++face) {
      query.addPoint(S2CellId.fromFace(face).toPoint());
    }
    assertTrue(query.getConvexHull().isFull());
  }

  public void testSimplePolyline() {
    // A polyline is handled identically to a point set, so there is no need
    // for special testing other than code coverage.
    S2Polyline polyline = makePolyline("0:1, 0:9, 1:6, 2:6, 3:10, 4:10, 5:5, 4:0, 3:0, 2:5, 1:5");
    S2ConvexHullQuery query = new S2ConvexHullQuery();
    query.addPolyline(polyline);
    S2Loop result = query.getConvexHull();
    S2Loop expectedResult = makeLoop("0:1, 0:9, 3:10, 4:10, 5:5, 4:0, 3:0");
    assertTrue(result.boundaryEquals(expectedResult));
  }

  private static void checkNorthPoleLoop(S1Angle radius, int numVertices) {
    // If the radius is very close to 90, then it's hard to predict whether the
    // result will be the full loop or not.
    Preconditions.checkArgument(Math.abs(radius.radians() - Math.PI / 2.0) >= 1e-15);

    S2ConvexHullQuery query = new S2ConvexHullQuery();
    S2Loop loop = S2Loop.makeRegularLoop(new S2Point(0, 0, 1), radius, numVertices);
    query.addLoop(loop);
    S2Loop result = query.getConvexHull();
    if (radius.radians() > Math.PI / 2.0) {
      assertTrue(result.isFull());
    } else {
      assertTrue(result.boundaryEquals(loop));
    }
  }

  public void testLoopsAroundNorthPole() {
    // Test loops of various sizes around the north pole.
    checkNorthPoleLoop(S1Angle.degrees(1), 3);
    checkNorthPoleLoop(S1Angle.degrees(89), 3);

    // The following two loops should yield the full loop.
    checkNorthPoleLoop(S1Angle.degrees(91), 3);
    checkNorthPoleLoop(S1Angle.degrees(179), 3);

    checkNorthPoleLoop(S1Angle.degrees(10), 100);
    checkNorthPoleLoop(S1Angle.degrees(89), 1000);
  }

  public void testPointsInsideHull() {
    // Repeatedly build the convex hull of a set of points, then add more points
    // inside that loop and build the convex hull again.  The result should
    // always be the same.
    int kIters = 1000;
    for (int iter = 0; iter < kIters; ++iter) {
      rand.setSeed(iter + 1); // Easier to reproduce a specific case.

      // Choose points from within a cap of random size, up to but not including
      // an entire hemisphere.
      S2Cap cap = getRandomCap(1e-15, 1.999 * Math.PI);
      S2ConvexHullQuery query = new S2ConvexHullQuery();
      int numPoints1 = rand.nextInt(100) + 3;
      for (int i = 0; i < numPoints1; ++i) {
        query.addPoint(samplePoint(cap));
      }
      S2Loop hull = query.getConvexHull();

      // When the convex hull is nearly a hemisphere, the algorithm sometimes
      // returns a full cap instead.  This is because it first computes a
      // bounding rectangle for all the input points/edges and then converts it
      // to a bounding cap, which sometimes yields a non-convex cap (radius
      // larger than 90 degrees).  This should not be a problem in practice
      // (since most convex hulls are not hemispheres), but in order to make this
      // test pass reliably it means that we need to reject convex hulls whose
      // bounding cap (when computed from a bounding rectangle) is not convex.
      //
      // TODO(user): This test can still fail (about 1 iteration in 500,000)
      // because the S2LatLngRect::GetCapBound implementation does not guarantee
      // that A.Contains(B) implies A.GetCapBound().Contains(B.GetCapBound()).
      if (hull.getCapBound().height() >= 1) {
        continue;
      }

      // Otherwise, add more points inside the convex hull.
      int numPoints2 = 1000;
      for (int i = 0; i < numPoints2; ++i) {
        S2Point p = samplePoint(cap);
        if (hull.contains(p)) {
          query.addPoint(p);
        }
      }
      // Finally, build a new convex hull and check that it hasn't changed.
      S2Loop hull2 = query.getConvexHull();
      assertTrue(hull2.boundaryEquals(hull));
    }
  }
}
