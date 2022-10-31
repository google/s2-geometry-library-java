/*
 * Copyright 2014 Google Inc.
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

import static java.lang.Math.atan2;
import static java.lang.Math.ceil;
import static java.lang.Math.cos;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.pow;
import static java.lang.Math.sin;
import static java.lang.Math.sqrt;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import java.util.List;
import java.util.Random;

/**
 * A simple class that generates "Koch snowflake" fractals (see Wikipedia for an introduction).
 * There is an option to control the fractal dimension (between 1.0 and 2.0); values between 1.02
 * and 1.50 are reasonable simulations of various coastlines. The default dimension (about 1.26)
 * corresponds to the standard Koch snowflake. (The west coast of Britain has a fractal dimension of
 * approximately 1.25)
 *
 * <p>The fractal is obtained by starting with an equilateral triangle and recursively subdividing
 * each edge into four segments of equal length. Therefore the shape at level 'n' consists of 3 *
 * (4^n) edges. Multi-level fractals are also supported: if you set minLevel() to a non-negative
 * value, then the recursive subdivision has an equal probability of stopping at any of the levels
 * between the given min and max (inclusive). This yields a fractal where the perimeter of the
 * original triangle is approximately equally divided between fractals at the various possible
 * levels. If there are k distinct levels {min, ..., max}, the expected number of edges at each
 * level 'i' is approximately 3 * (4 ^ i) / k.
 */
public class S2FractalBuilder {
  private int maxLevel = -1;

  /** Value set by user. */
  private int minLevelArg = -1;

  /** Actual min level (depends on maxLevel.) */
  private int minLevel = -1;

  /** Standard Koch curve */
  private double dimension = log(4) / log(3);

  /** The ratio of the sub-edge length to the original edge length at each subdivision step. */
  private double edgeFraction = 0;

  /**
   * The distance from the original edge to the middle vertex at each subdivision step, as a
   * fraction of the original edge length.
   */
  private double offsetFraction = 0;

  private Random rand;

  /** You must call setMaxLevel() or setLevelForApproxMaxEdges() before calling makeLoop(). */
  public S2FractalBuilder(Random rand) {
    this.rand = rand;
    computeOffsets();
  }

  /** Sets the maximum subdivision level for the fractal (see above). */
  public void setMaxLevel(int maxLevel) {
    Preconditions.checkArgument(maxLevel >= 0);
    this.maxLevel = maxLevel;
    computeMinLevel();
  }

  /**
   * Sets the minimum subdivision level for the fractal (see above). The default value of -1 causes
   * the min and max levels to be the same. A minLevel of 0 should be avoided since this creates a
   * significant chance that none of the three original edges will be subdivided at all.
   */
  public void setMinLevel(int minLevelArg) {
    Preconditions.checkArgument(minLevelArg >= -1);
    this.minLevelArg = minLevelArg;
    computeMinLevel();
  }

  private void computeMinLevel() {
    if (minLevelArg >= 0 && minLevelArg <= maxLevel) {
      minLevel = minLevelArg;
    } else {
      minLevel = maxLevel;
    }
  }

  /**
   * Sets the fractal dimension. The default value of approximately 1.26 corresponds to the standard
   * Koch curve. The value must lie in the range [1.0, 2.0).
   */
  public void setFractalDimension(double dimension) {
    Preconditions.checkArgument(dimension >= 1.0);
    Preconditions.checkArgument(dimension <= 2.0);
    this.dimension = dimension;
    computeOffsets();
  }

  private void computeOffsets() {
    edgeFraction = pow(4.0, -1.0 / dimension);
    offsetFraction = sqrt(edgeFraction - 0.25);
  }

  /**
   * The following two functions set the min and/or max level to produce approximately the given
   * number of edges. (The values are rounded to a nearby value of 3 * (4 ^ n).)
   */
  public void setLevelForApproxMinEdges(int minEdges) {
    setMinLevel(levelFromEdges(minEdges));
  }

  public void setLevelForApproxMaxEdges(int maxEdges) {
    setMaxLevel(levelFromEdges(maxEdges));
  }

  /** Returns level from values in the range [1.5 * (4 ^ n), 6 * (4 ^ n)]. */
  private static int levelFromEdges(int edges) {
    return (int) ceil(0.5 * log((double) edges / 3) / log(2));
  }

  /**
   * Returns a lower bound on the ratio (Rmin / R), where 'R' is the radius passed to makeLoop(),
   * and 'Rmin' is the minimum distance from the fractal boundary to its center. This can be used to
   * inscribe another geometric figure within the fractal without intersection.
   */
  public double minRadiusFactor() {
    // The minimum radius is attained at one of the vertices created by the first subdivision
    // step as long as the dimension is not too small (at least
    // kMinDimensionForMinRadiusAtLevel1, see below).  Otherwise we fall back on the incircle
    // radius of the original triangle, which is always a lower bound (and is attained when
    // dimension = 1).
    //
    // The value below was obtained by letting AE be an original triangle edge, letting ABCDE be
    // the corresponding polyline after one subdivision step, and then letting BC be tangent to
    // the inscribed circle at the center of the fractal O.  This gives rise to a pair of
    // similar triangles whose edge length ratios can be used to solve for the corresponding
    // "edge fraction".  This method is slightly conservative because it is computed using
    // planar rather than spherical geometry.  The value below is equal to
    // -log(4)/log((2 + cbrt(2) - cbrt(4))/6).
    double kMinDimensionForMinRadiusAtLevel1 = 1.0852230903040407;
    if (dimension >= kMinDimensionForMinRadiusAtLevel1) {
      return sqrt(1 + 3 * edgeFraction * (edgeFraction - 1));
    }
    return 0.5;
  }

  /**
   * Returns the ratio (Rmax / R), where 'R' is the radius passed to makeLoop() and 'Rmax' is the
   * maximum distance from the fractal boundary to its center. This can be used to inscribe the
   * fractal within some other geometric figure without intersection.
   */
  public double maxRadiusFactor() {
    // The maximum radius is always attained at either an original triangle vertex or at a middle
    // vertex from the first subdivision step.
    return max(1.0, offsetFraction * sqrt(3) + 0.5);
  }

  private void getR2Vertices(List<R2Vector> vertices) {
    // The Koch "snowflake" consists of three Koch curves whose initial edges form an
    // equilateral triangle.
    R2Vector v0 = new R2Vector(1.0, 0.0);
    R2Vector v1 = new R2Vector(-0.5, sqrt(3) / 2);
    R2Vector v2 = new R2Vector(-0.5, -sqrt(3) / 2);
    getR2VerticesHelper(v0, v1, 0, vertices);
    getR2VerticesHelper(v1, v2, 0, vertices);
    getR2VerticesHelper(v2, v0, 0, vertices);
  }

  /**
   * Given the two endpoints (v0, v4) of an edge, recursively subdivide the edge to the desired
   * level, and insert all vertices of the resulting curve up to but not including the endpoint
   * "v4".
   */
  private void getR2VerticesHelper(R2Vector v0, R2Vector v4, int level, List<R2Vector> vertices) {
    // The second expression should return 'true' once every (maxLevel - level + 1) calls.
    if (level >= minLevel && (rand.nextInt(maxLevel - level + 1) == 0)) {
      // Stop subdivision at this level.
      vertices.add(v0);
      return;
    }
    // Otherwise compute the intermediate vertices v1, v2, and v3.
    R2Vector dir = v4.sub(v0);
    // v1 = v0 + edgeFraction * dir
    R2Vector v1 = v0.add(dir.mul(edgeFraction));
    // v2 = 0.5 * (v0 + v4) - offsetFraction * dir.ortho()
    R2Vector v2 = v0.add(v4).mul(0.5).sub(dir.ortho().mul(offsetFraction));
    // v3 = v4 - edgeFraction * dir
    R2Vector v3 = v4.sub(dir.mul(edgeFraction));

    // And recurse on the four sub-edges.
    getR2VerticesHelper(v0, v1, level + 1, vertices);
    getR2VerticesHelper(v1, v2, level + 1, vertices);
    getR2VerticesHelper(v2, v3, level + 1, vertices);
    getR2VerticesHelper(v3, v4, level + 1, vertices);
  }

  /**
   * Returns a fractal loop centered around the a-axis of the given coordinate frame, with the first
   * vertex in the direction of the positive x-axis, and the given nominal radius.
   */
  public S2Loop makeLoop(Matrix frame, S1Angle nominalRadius) {
    return new S2Loop(makeVertices(frame, nominalRadius));
  }

  /** As {@link #makeLoop(Matrix, S1Angle)} except it returns the vertices instead of loop. */
  public List<S2Point> makeVertices(Matrix frame, S1Angle nominalRadius) {
    List<R2Vector> r2Vertices = Lists.newArrayList();
    getR2Vertices(r2Vertices);
    List<S2Point> vertices = Lists.newArrayList();
    for (int i = 0; i < r2Vertices.size(); ++i) {
      // Convert each vertex to polar coordinates.
      R2Vector v = r2Vertices.get(i);
      double theta = atan2(v.y(), v.x());
      double radius = nominalRadius.radians() * v.norm();

      // We construct the loop in the given frame coordinates, with the center at (0, 0, 1).  For
      // a loop of radius 'r', the loop vertices have the form (x, y, z) where x^2 + y^2 = sin(r)
      // and z = cos(r).  The distance on the sphere (arc length) from each vertex to the center
      // is acos(cos(r)) = r.
      double z = cos(radius);
      double r = sin(radius);
      S2Point p = new S2Point(r * cos(theta), r * sin(theta), z);
      vertices.add(S2.rotate(p, frame));
    }
    return vertices;
  }
}
