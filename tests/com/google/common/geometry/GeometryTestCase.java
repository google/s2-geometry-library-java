/*
 * Copyright 2005 Google Inc.
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
import com.google.common.annotations.GwtIncompatible;
import com.google.common.base.Splitter;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

import junit.framework.TestCase;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * Common code for geometry tests.
 */
@GwtCompatible(emulated = true)
public strictfp class GeometryTestCase extends TestCase {

  public Random rand;

  @Override
  protected void setUp() {
    rand = new Random(123455);
  }

  public void assertDoubleNear(double a, double b) {
    assertDoubleNear(a, b, 1e-9);
  }

  public void assertDoubleNear(double a, double b, double error) {
    assertTrue(a + error > b);
    assertTrue(a < b + error);
  }

  /**
   * Checks that the 3D distance between {@code expected} and {@code actual} is
   * at most {@code eps} units.
   */
  public void assertEquals(S2Point expected, S2Point actual, double eps) {
    assertTrue("expected: " + expected + " but was: " + actual,
        expected.getDistance2(actual) < eps * eps);
  }

  // maybe these should be put in a special testing util class
  /** Return a random unit-length vector. */
  public S2Point randomPoint() {
    return S2Point.normalize(new S2Point(
        2 * rand.nextDouble() - 1,
        2 * rand.nextDouble() - 1,
        2 * rand.nextDouble() - 1));
  }
  
  /**
   * Return a random cell id at the given level or at a randomly chosen level.
   * The distribution is uniform over the space of cell ids, but only
   * approximately uniform over the surface of the sphere.
   */
  public S2CellId getRandomCellId(int level) {
    int face = random(S2CellId.NUM_FACES);
    long pos = rand.nextLong() & ((1L << (2 * S2CellId.MAX_LEVEL)) - 1);
    return S2CellId.fromFacePosLevel(face, pos, level);
  }

  public S2CellId getRandomCellId() {
    return getRandomCellId(random(S2CellId.MAX_LEVEL + 1));
  }

  /** Returns a randomly selected value between a and b. */
  public double uniform(double a, double b) {
    return a + (b - a) * rand.nextDouble();
  }

  /** Returns true randomly, approximately one time in {@code n} calls. */
  public boolean oneIn(int n) {
    return rand.nextInt(n) == 0;
  }

  int random(int n) {
    if (n == 0) {
      return 0;
    }
    return rand.nextInt(n);
  }

  /**
   * Picks a "base" uniformly from range [0,maxLog] and then return "base" random bits. The effect
   * is to pick a number in the range [0,2^maxLog-1] with bias towards smaller numbers.
   */
  int skewed(int maxLog) {
    final int base = Math.abs(rand.nextInt()) % (maxLog + 1);
    // if (!base) return 0; // if 0==base, we & with 0 below.
    //
    // this distribution differs slightly from ACMRandom's Skewed,
    // since 0 occurs approximately 3 times more than 1 here, and
    // ACMRandom's Skewed never outputs 0.
    return rand.nextInt() & ((1 << base) - 1);
  }

  /**
   * As {@link #checkCovering(S2Region, S2CellUnion, boolean, S2CellId)}, but creates a default and
   * invalid S2CellId for the last argument.
   */
  void checkCovering(S2Region region, S2CellUnion covering, boolean checkTight) {
    checkCovering(region, covering, checkTight, new S2CellId());
  }

  /**
   * Checks that "covering" completely covers the given region. If "checkTight" is true, also checks
   * that it does not contain any cells that do not intersect the given region. ("id" is only used
   * internally.)
   */
  void checkCovering(S2Region region, S2CellUnion covering, boolean checkTight, S2CellId id) {
    if (!id.isValid()) {
      for (int face = 0; face < 6; ++face) {
        checkCovering(region, covering, checkTight, S2CellId.fromFacePosLevel(face, 0, 0));
      }
      return;
    }

    if (!region.mayIntersect(new S2Cell(id))) {
      // If region does not intersect id, then neither should the covering.
      if (checkTight) {
        assertTrue(!covering.intersects(id));
      }

    } else if (!covering.contains(id)) {
      // The region may intersect id, but we can't assert that the covering
      // intersects id because we may discover that the region does not actually
      // intersect upon further subdivision. (MayIntersect is not exact.)
      assertTrue(!region.contains(new S2Cell(id)));
      assertTrue(!id.isLeaf());
      S2CellId end = id.childEnd();
      for (S2CellId child = id.childBegin(); !child.equals(end); child = child.next()) {
        checkCovering(region, covering, checkTight, child);
      }
    }
  }

  S2Cap getRandomCap(double minArea, double maxArea) {
    double capArea = maxArea
      * Math.pow(minArea / maxArea, rand.nextDouble());
    assertTrue(capArea >= minArea && capArea <= maxArea);

    // The surface area of a cap is 2*Pi times its height.
    return S2Cap.fromAxisArea(randomPoint(), capArea);
  }

  /** Returns a polygon with given center, number of concentric loops, and vertices per loop. */
  static S2Polygon concentricLoopsPolygon(S2Point center, int numLoops, int numVerticesPerLoop) {
    Matrix3x3 m = Matrix3x3.fromCols(S2.getFrame(center));
    List<S2Loop> loops = new ArrayList<S2Loop>(numLoops);
    for (int li = 0; li < numLoops; ++li) {
      List<S2Point> vertices = new ArrayList<S2Point>(numVerticesPerLoop);
      double radius = 0.005 * (li + 1) / numLoops;
      double radianStep = 2 * S2.M_PI / numVerticesPerLoop;
      for (int vi = 0; vi < numVerticesPerLoop; ++vi) {
        double angle = vi * radianStep;
        S2Point p = new S2Point(radius * Math.cos(angle), radius * Math.sin(angle), 1);
        vertices.add(S2.rotate(p, m));
      }
      loops.add(new S2Loop(vertices));
    }
    return new S2Polygon(loops);
  }

  S2Point samplePoint(S2Cap cap) {
    // We consider the cap axis to be the "z" axis. We choose two other axes to
    // complete the coordinate frame.

    S2Point z = cap.axis();
    S2Point x = z.ortho();
    S2Point y = S2Point.crossProd(z, x);

    // The surface area of a spherical cap is directly proportional to its
    // height. First we choose a random height, and then we choose a random
    // point along the circle at that height.

    double h = rand.nextDouble() * cap.height();
    double theta = 2 * S2.M_PI * rand.nextDouble();
    double r = Math.sqrt(h * (2 - h)); // Radius of circle.

    // (cos(theta)*r*x + sin(theta)*r*y + (1-h)*z).Normalize()
    return S2Point.normalize(S2Point.add(
        S2Point.add(S2Point.mul(x, Math.cos(theta) * r), S2Point.mul(y, Math.sin(theta) * r)),
        S2Point.mul(z, (1 - h))));
  }
  
  
  /** Return a random point within the given S2LatLngRect. */
  S2Point samplePoint(S2LatLngRect rect) {
    // First choose a latitude uniformly with respect to area on the sphere.
    double sinLo = Math.sin(rect.lat().lo());
    double sinHi = Math.sin(rect.lat().hi());
    double lat = Math.asin(sinLo + rand.nextDouble() * (sinHi - sinLo));
    
    // Now choose longitude uniformly within the given range.
    double lng = rect.lng().lo() + rand.nextDouble() * rect.lng().getLength();
    
    return S2LatLng.fromRadians(lat, lng).normalized().toPoint();
  }

  static S2LatLngRect parseVertices(String str, List<S2Point> vertices) {
    if (str == null) {
      return null;
    }

    S2LatLngRect.Builder builder = S2LatLngRect.Builder.empty();
    for (String token : Splitter.on(',').trimResults().omitEmptyStrings().split(str)) {
      int colon = token.indexOf(':');
      if (colon == -1) {
        throw new IllegalArgumentException(
            "Illegal string:" + token + ". Should look like '35:20'");
      }
      double lat = Double.parseDouble(token.substring(0, colon));
      double lng = Double.parseDouble(token.substring(colon + 1));
      vertices.add(S2LatLng.fromDegrees(lat, lng).toPoint());
      builder.addPoint(S2LatLng.fromDegrees(lat, lng));
    }

    return builder.build();
  }

  static S2Point makePoint(String str) {
    List<S2Point> vertices = Lists.newArrayList();
    parseVertices(str, vertices);
    return Iterables.getOnlyElement(vertices);
  }

  /**
   * Given a string of latitude-longitude coordinates in degrees, returns a newly allocated loop.
   * Example of the input format: {@code -20:150, 10:-120, 0.123:-170.652}.
   *
   * <p>The exact strings "empty" or "full" create the empty or full loop, respectively.
   */
  static S2Loop makeLoop(String str) {
    if ("empty".equalsIgnoreCase(str)) {
      return S2Loop.empty();
    } else if ("full".equalsIgnoreCase(str)) {
      return S2Loop.full();
    }
    List<S2Point> vertices = Lists.newArrayList();
    parseVertices(str, vertices);
    return new S2Loop(vertices);
  }

  /** Convert a distance on the Earth's surface to an angle. */
  static S1Angle kmToAngle(double km) {
    return metersToAngle(1000 * km);
  }

  /** Convert a distance on the Earth's surface to an angle. */
  static S1Angle metersToAngle(double meters) {
    return S1Angle.radians(meters / S2LatLng.EARTH_RADIUS_METERS);
  }

  /**
   * Given a sequence of loops separated by semicolons, returns a newly allocated polygon.  Loops
   * are automatically normalized by inverting them if necessary so that they enclose at most half
   * of the unit sphere.
   *
   * <p>Historically this was once a requirement of polygon loops.  It also hides the problem that
   * if the user thinks of the coordinates as X:Y rather than LAT:LNG, it yields a loop with the
   * opposite orientation.
   *
   * <p>Examples of the input format:
   * <ul>
   * <li>"10:20, 90:0, 20:30" <=> one loop
   * <li>"10:20, 90:0, 20:30; 5.5:6.5, -90:-180, -15.2" <=> >two loops
   * <li>"" <=> the empty string results in the empty polygon (consisting of no loops)
   * <li>"full" <=> the full polygon (consisting of one full loop)
   * <li>"empty" <=> **INVALID** the empty polygon has no loops.
   * </ul>
   */
  static S2Polygon makePolygon(String str) {
    return internalMakePolygon(str, true);
  }

  /**
   * Like makePolygon(), except that it does not normalize loops (i.e., it gives you exactly what
   * you asked for).
   */
  static S2Polygon makeVerbatimPolygon(String str) {
    return internalMakePolygon(str, false);
  }

  private static final S2Polygon internalMakePolygon(String str, boolean normalizeLoops) {
    List<S2Loop> loops = Lists.newArrayList();

    for (String token : Splitter.on(';').omitEmptyStrings().split(str)) {
      S2Loop loop = makeLoop(token);
      if (normalizeLoops) {
        loop.normalize();
      }
      loops.add(loop);
    }

    return new S2Polygon(loops);
  }

  static S2Polyline makePolyline(String str) {
    List<S2Point> vertices = Lists.newArrayList();
    parseVertices(str, vertices);
    return new S2Polyline(vertices);
  }

  /**
   * Returns the result of encoding and immediately decoding the given value.
   */
  @SuppressWarnings("unchecked")
  @GwtIncompatible("ByteArrayInputStream")
  static <E> E encodeDecode(Serializable value) throws Exception {
    ByteArrayOutputStream bytes = new ByteArrayOutputStream();
    new ObjectOutputStream(bytes).writeObject(value);
    ObjectInputStream in = new ObjectInputStream(new ByteArrayInputStream(bytes.toByteArray()));
    return (E) in.readObject();
  }
}
