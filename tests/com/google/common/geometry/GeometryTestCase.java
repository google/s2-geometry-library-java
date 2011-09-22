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

import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;

import junit.framework.TestCase;

import java.util.List;
import java.util.Random;

public strictfp class GeometryTestCase extends TestCase {

  public Random rand;

  @Override
  protected void setUp() {
    rand = new Random(123456);
  }

  public void assertDoubleNear(double a, double b) {
    assertDoubleNear(a, b, 1e-9);
  }

  public void assertDoubleNear(double a, double b, double error) {
    assertTrue(a + error > b);
    assertTrue(a < b + error);
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
   * Return a right-handed coordinate frame (three orthonormal vectors). Returns
   * an array of three points: x,y,z
   */
  public ImmutableList<S2Point> getRandomFrame() {
    S2Point p0 = randomPoint();
    S2Point p1 = S2Point.normalize(S2Point.crossProd(p0, randomPoint()));
    S2Point p2 = S2Point.normalize(S2Point.crossProd(p0, p1));
    return ImmutableList.of(p0, p1, p2);
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

  int random(int n) {
    if (n == 0) {
      return 0;
    }
    return rand.nextInt(n);
  }


  // Pick "base" uniformly from range [0,maxLog] and then return
  // "base" random bits. The effect is to pick a number in the range
  // [0,2^maxLog-1] with bias towards smaller numbers.
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
   * Checks that "covering" completely covers the given region. If "check_tight"
   * is true, also checks that it does not contain any cells that do not
   * intersect the given region. ("id" is only used internally.)
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

  static void parseVertices(String str, List<S2Point> vertices) {
    if (str == null) {
      return;
    }

    for (String token : Splitter.on(',').split(str)) {
      int colon = token.indexOf(':');
      if (colon == -1) {
        throw new IllegalArgumentException(
            "Illegal string:" + token + ". Should look like '35:20'");
      }
      double lat = Double.parseDouble(token.substring(0, colon));
      double lng = Double.parseDouble(token.substring(colon + 1));
      vertices.add(S2LatLng.fromDegrees(lat, lng).toPoint());
    }
  }

  static S2Point makePoint(String str) {
    List<S2Point> vertices = Lists.newArrayList();
    parseVertices(str, vertices);
    return Iterables.getOnlyElement(vertices);
  }

  static S2Loop makeLoop(String str) {
    List<S2Point> vertices = Lists.newArrayList();
    parseVertices(str, vertices);
    return new S2Loop(vertices);
  }

  static S2Polygon makePolygon(String str) {
    List<S2Loop> loops = Lists.newArrayList();

    for (String token : Splitter.on(';').omitEmptyStrings().split(str)) {
      S2Loop loop = makeLoop(token);
      loop.normalize();
      loops.add(loop);
    }

    return new S2Polygon(loops);
  }

  static S2Polyline makePolyline(String str) {
    List<S2Point> vertices = Lists.newArrayList();
    parseVertices(str, vertices);
    return new S2Polyline(vertices);
  }
}
