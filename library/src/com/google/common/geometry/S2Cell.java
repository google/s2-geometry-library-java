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

import static com.google.common.geometry.S2.DBL_EPSILON;
import static com.google.common.geometry.S2.M_1_PI;
import static com.google.common.geometry.S2.M_PI_2;
import static com.google.common.geometry.S2.M_PI_4;
import static com.google.common.geometry.S2.area;
import static com.google.common.geometry.S2.posToIJ;
import static com.google.common.geometry.S2.posToOrientation;
import static com.google.common.geometry.S2Projections.PROJ;
import static java.lang.Math.asin;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;

import com.google.common.collect.Ordering;
import com.google.common.primitives.Doubles;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.Serializable;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;

/**
 * An S2Cell is an S2Region object that represents a cell. Unlike S2CellIds, it supports efficient
 * containment and intersection tests. However, it is also a more expensive representation.
 *
 * @author danieldanciu@google.com (Daniel Danciu) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@JsType
public final strictfp class S2Cell implements S2Region, Serializable {
  byte face;
  byte level;
  byte orientation;
  S2CellId cellId;
  double uMin;
  double uMax;
  double vMin;
  double vMax;

  /** Default constructor used only internally. */
  S2Cell() {}

  /**
   * An S2Cell always corresponds to a particular S2CellId. The other constructors are just
   * convenience methods.
   */
  @JsIgnore
  public S2Cell(S2CellId id) {
    init(id);
  }

  /** Returns the cell corresponding to the given S2 cube face. */
  public static S2Cell fromFace(int face) {
    return new S2Cell(S2CellId.fromFace(face));
  }

  /**
   * Returns a cell given its face (range 0..5), Hilbert curve position within that face (an
   * unsigned integer with {@link S2CellId#POS_BITS} bits), and level (range 0..MAX_LEVEL). The
   * given position will be modified to correspond to the Hilbert curve position at the center of
   * the returned cell. This is a static function rather than a constructor in order to indicate
   * what the arguments represent.
   */
  public static S2Cell fromFacePosLevel(int face, long pos, int level) {
    return new S2Cell(S2CellId.fromFacePosLevel(face, pos, level));
  }

  /** Convenience method to construct a leaf S2Cell containing the given point. */
  @JsIgnore
  public S2Cell(S2Point p) {
    init(S2CellId.fromPoint(p));
  }

  /** Convenience method to construct a leaf S2Cell containing the given lat,lng. */
  @JsIgnore
  public S2Cell(S2LatLng ll) {
    init(S2CellId.fromLatLng(ll));
  }

  /** Returns the S2CellId of this cell. */
  public S2CellId id() {
    return cellId;
  }

  /** Returns the face this cell is located on, in the range 0..5. */
  public int face() {
    return face;
  }

  /** Returns the level of this cell, in the range 0..S2CellId.MAX_LEVEL. */
  public byte level() {
    return level;
  }

  /** Returns the cell orientation, in the range 0..3. */
  public byte orientation() {
    return orientation;
  }

  /** Returns true if this cell is a leaf-cell, i.e. it has no children. */
  public boolean isLeaf() {
    return level == S2CellId.MAX_LEVEL;
  }

  /**
   * As {@link #getVertexRaw(int)}, except the point is normalized to unit length. For convenience,
   * the argument is reduced modulo 4 to the range [0..3].
   */
  public S2Point getVertex(int k) {
    return getVertexRaw(k).normalize();
  }

  /**
   * Returns the k<sup>th</sup> vertex of the cell (k = 0,1,2,3). Vertices are returned in CCW order
   * (lower left, lower right, upper right, upper left in the UV plane). The points are not
   * necessarily unit length. For convenience, the argument is reduced modulo 4 to the range [0..3].
   */
  public S2Point getVertexRaw(int k) {
    k &= 3;
    // Vertices are returned in the order SW, SE, NE, NW.
    return S2Projections.faceUvToXyz(
        face, ((k >> 1) ^ (k & 1)) == 0 ? uMin : uMax, (k >> 1) == 0 ? vMin : vMax);
  }

  /** As {@link #getEdgeRaw(int)}, except the point is normalized to unit length. */
  public S2Point getEdge(int k) {
    return getEdgeRaw(k).normalize();
  }

  /**
   * Returns the inward-facing normal of the great circle passing through the edge from vertex k to
   * vertex k+1 (mod 4). The normals returned by getEdgeRaw are not necessarily unit length. For
   * convenience, the argument is reduced modulo 4 to the range [0..3].
   */
  public S2Point getEdgeRaw(int k) {
    switch (k & 3) {
      case 0:
        return S2Projections.getVNorm(face, vMin); // Bottom
      case 1:
        return S2Projections.getUNorm(face, uMax); // Right
      case 2:
        return S2Projections.getVNorm(face, vMax).neg(); // Top
      default:
        return S2Projections.getUNorm(face, uMin).neg(); // Left
    }
  }

  /** As {@link S2CellId#getSizeIJ(int)}, using the level of this cell. */
  public int getSizeIJ() {
    return S2CellId.getSizeIJ(level());
  }

  /**
   * Returns true if this is not a leaf cell, in which case the array, which must contain at least
   * four non-null cells in indices 0..3, will be set to the four children of this cell in traversal
   * order. Otherwise, if this is a leaf cell, false is returned without touching the array.
   *
   * <p>This method is equivalent to the following:
   *
   * <pre>
   * for (pos = 0, id = childBegin(); !id.equals(childEnd()); id = id.next(), ++pos) {
   *   children[i].init(id);
   * }
   * </pre>
   *
   * <p>except that it is more than two times faster.
   */
  @CanIgnoreReturnValue
  public boolean subdivide(S2Cell[] children) {
    if (cellId.isLeaf()) {
      return false;
    }

    // Create four children with the appropriate bounds.
    S2CellId id = cellId.childBegin();
    R2Vector mid = getCenterUV();
    double uMid = mid.x();
    double vMid = mid.y();
    for (int pos = 0; pos < 4; ++pos, id = id.next()) {
      S2Cell child = children[pos];
      child.face = face;
      child.level = (byte) (level + 1);
      child.orientation = (byte) (orientation ^ posToOrientation(pos));
      child.cellId = id;
      // We want to split the cell in half in "u" and "v". To decide which side to set equal to the
      // midpoint value, we look at cell's (i,j) position within its parent. The index for "i" is in
      // bit 1 of ij.
      int ij = posToIJ(orientation, pos);
      // The dimension 0 index (i/u) is in bit 1 of ij.
      if ((ij & 0x2) != 0) {
        child.uMin = uMid;
        child.uMax = uMax;
      } else {
        child.uMin = uMin;
        child.uMax = uMid;
      }
      // The dimension 1 index (j/v) is in bit 0 of ij.
      if ((ij & 0x1) != 0) {
        child.vMin = vMid;
        child.vMax = vMax;
      } else {
        child.vMin = vMin;
        child.vMax = vMid;
      }
    }
    return true;
  }

  /**
   * Return the normalized direction vector corresponding to the center in (s,t)-space of the given
   * cell. This is the point at which the cell is divided into four subcells; it is not necessarily
   * the centroid of the cell in (u,v)-space or (x,y,z)-space.
   */
  public S2Point getCenter() {
    return getCenterRaw().normalize();
  }

  /**
   * Return the direction vector corresponding to the center in (s,t)-space of the given cell. This
   * is the point at which the cell is divided into four subcells; it is not necessarily the
   * centroid of the cell in (u,v)-space or (x,y,z)-space. The point returned by GetCenterRaw is not
   * necessarily unit length.
   */
  public S2Point getCenterRaw() {
    return cellId.toPointRaw();
  }

  /** Returns the bounds of this cell in (u,v)-space. */
  public R2Rect getBoundUV() {
    R2Rect rect = new R2Rect();
    setBoundUV(rect);
    return rect;
  }

  /**
   * Sets the bounds of this cell in (u,v)-space into 'bound'.
   *
   * <p>Package private to avoid leaking object mutation outside the api.
   */
  void setBoundUV(R2Rect bound) {
    bound.x().set(uMin, uMax);
    bound.y().set(vMin, vMax);
  }

  /**
   * Return the center of the cell in (u,v) coordinates (see {@code S2Projections}). Note that the
   * center of the cell is defined as the point at which it is recursively subdivided into four
   * children; in general, it is not at the midpoint of the (u,v) rectangle covered by the cell.
   */
  public R2Vector getCenterUV() {
    return cellId.getCenterUV();
  }

  /** Return the average area in steradians for cells at the given level. */
  @JsMethod(name = "averageAreaAtLevel")
  public static double averageArea(int level) {
    return PROJ.avgArea.getValue(level);
  }

  /**
   * Return the average area in steradians of cells at this level. This is accurate to within a
   * factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is extremely cheap to compute.
   */
  public double averageArea() {
    return averageArea(level);
  }

  /**
   * Return the approximate area of this cell in steradians. This method is accurate to within 3%
   * percent for all cell sizes and accurate to within 0.1% for cells at level 5 or higher (i.e.
   * 300km square or smaller). It is moderately cheap to compute.
   */
  public double approxArea() {

    // All cells at the first two levels have the same area.
    if (level < 2) {
      return averageArea(level);
    }

    // First, compute the approximate area of the cell when projected perpendicular to its normal.
    // The cross product of its diagonals gives the normal, and the length of the normal is twice
    // the projected area.
    double flatArea =
        0.5 * getVertex(2).sub(getVertex(0)).crossProd(getVertex(3).sub(getVertex(1))).norm();

    // Now, compensate for the curvature of the cell surface by pretending that the cell is shaped
    // like a spherical cap. The ratio of the area of a spherical cap to the area of its projected
    // disc turns out to be 2 / (1 + sqrt(1 - r*r)) where "r" is the radius of the disc. For
    // example, when r=0 the ratio is 1, and when r=1 the ratio is 2. Here we set Pi*r*r == flatArea
    // to find the equivalent disc.
    return flatArea * 2 / (1 + sqrt(1 - min(M_1_PI * flatArea, 1.0)));
  }

  /**
   * Return the area in steradians of this cell as accurately as possible. This method is more
   * expensive but it is accurate to 6 digits of precision even for leaf cells (whose area is
   * approximately 1e-18).
   */
  public double exactArea() {
    S2Point v0 = getVertex(0);
    S2Point v1 = getVertex(1);
    S2Point v2 = getVertex(2);
    S2Point v3 = getVertex(3);
    return area(v0, v1, v2) + area(v0, v2, v3);
  }

  // //////////////////////////////////////////////////////////////////////
  // S2Region interface (see {@code S2Region} for details):

  // NOTE: This should be marked as @Override, but clone() isn't present in GWT's version of Object,
  // so we can't mark it as such.
  @SuppressWarnings("MissingOverride")
  public S2Region clone() {
    S2Cell clone = new S2Cell();
    clone.face = this.face;
    clone.level = this.level;
    clone.orientation = this.orientation;
    clone.uMin = this.uMin;
    clone.uMax = this.uMax;
    clone.vMin = this.vMin;
    clone.vMax = this.vMax;
    return clone;
  }

  @Override
  public S2Cap getCapBound() {
    // Use the cell center in (u,v)-space as the cap axis. This vector is very close to GetCenter()
    // and faster to compute. Neither one of these vectors yields the bounding cap with minimal
    // surface area, but they are both pretty close.
    //
    // It's possible to show that the two vertices that are furthest from the (u,v)-origin never
    // determine the maximum cap size (this is a possible future optimization).
    S2Point center = S2Projections.faceUvToXyz(face, getCenterUV()).normalize();
    S2Cap cap = S2Cap.fromAxisHeight(center, 0);
    for (int k = 0; k < 4; ++k) {
      cap = cap.addPoint(getVertex(k));
    }
    return cap;
  }

  /**
   * The 4 cells around the equator extend to +/-45 degrees latitude at the midpoints of their top
   * and bottom edges. The two cells covering the poles extend down to +/-35.26 degrees at their
   * vertices. The maximum error in this calculation is 0.5 * DBL_EPSILON.
   */
  private static final double POLE_MIN_LAT = asin(sqrt(1. / 3)) - 0.5 * DBL_EPSILON;

  @Override
  public S2LatLngRect getRectBound() {
    if (level > 0) {
      // Except for cells at level 0, the latitude and longitude extremes are attained at the
      // vertices. Furthermore, the latitude range is determined by one pair of diagonally opposite
      // vertices and the longitude range is determined by the other pair.
      //
      // We first determine which corner (i,j) of the cell has the largest absolute latitude. To
      // maximize latitude, we want to find the point in the cell that has the largest absolute
      // z-coordinate and the smallest absolute x- and y-coordinates. To do this we look at each
      // coordinate (u and v), and determine whether we want to minimize or maximize that coordinate
      // based on the axis direction and the cell's (u,v) quadrant.
      double u = uMin + uMax;
      double v = vMin + vMax;
      int i = (S2Projections.getUAxis(face).z == 0 ? (u < 0) : (u > 0)) ? 1 : 0;
      int j = (S2Projections.getVAxis(face).z == 0 ? (v < 0) : (v > 0)) ? 1 : 0;
      R1Interval lat =
          R1Interval.fromPointPair(
              S2LatLng.latitude(getPoint(i, j)).radians(),
              S2LatLng.latitude(getPoint(1 - i, 1 - j)).radians());
      S1Interval lng =
          S1Interval.fromPointPair(
              S2LatLng.longitude(getPoint(i, 1 - j)).radians(),
              S2LatLng.longitude(getPoint(1 - i, j)).radians());

      // We grow the bounds slightly to make sure that the bounding rectangle contains S2LatLng(P)
      // for any point P inside the loop L defined by the four *normalized* vertices. Note that
      // normalization of a vector can change its direction by up to 0.5 * DBL_EPSILON radians, and
      // it is not enough just to add normalize() calls to the code above because the latitude/
      // longitude ranges are not necessarily determined by diagonally opposite vertex pairs after
      // normalization.
      //
      // We would like to bound the amount by which the latitude/longitude of a contained point P
      // can exceed the bounds computed above. In the case of longitude, the normalization error can
      // change the direction of rounding leading to a maximum difference in longitude of 2 *
      // DBL_EPSILON. In the case of latitude, the normalization error can shift the latitude by up
      // to 0.5 * DBL_EPSILON and the other sources of error can cause the two latitudes to differ
      // by up to another 1.5 * DBL_EPSILON, which also leads to a maximum difference of 2 *
      // DBL_EPSILON.
      return new S2LatLngRect(lat, lng)
          .expanded(S2LatLng.fromRadians(2 * DBL_EPSILON, 2 * DBL_EPSILON))
          .polarClosure();
    }

    // The face centers are the +X, +Y, +Z, -X, -Y, -Z axes in that order.
    // assert (((face < 3) ? 1 : -1) == S2Projections.getNorm(face).get(face % 3));

    S2LatLngRect bound;
    switch (face) {
      case 0:
        bound =
            new S2LatLngRect(
                new R1Interval(-M_PI_4, M_PI_4), new S1Interval(-M_PI_4, M_PI_4));
        break;
      case 1:
        bound =
            new S2LatLngRect(
                new R1Interval(-M_PI_4, M_PI_4), new S1Interval(M_PI_4, 3 * M_PI_4));
        break;
      case 2:
        bound = new S2LatLngRect(new R1Interval(POLE_MIN_LAT, M_PI_2), S1Interval.full());
        break;
      case 3:
        bound =
            new S2LatLngRect(
                new R1Interval(-M_PI_4, M_PI_4),
                new S1Interval(3 * M_PI_4, -3 * M_PI_4));
        break;
      case 4:
        bound =
            new S2LatLngRect(
                new R1Interval(-M_PI_4, M_PI_4), new S1Interval(-3 * M_PI_4, -M_PI_4));
        break;
      default:
        bound = new S2LatLngRect(new R1Interval(-M_PI_2, -POLE_MIN_LAT), S1Interval.full());
        break;
    }

    // Finally, we expand the bound to account for the error when a point P is converted to an
    // S2LatLng to test for containment. (The bound should be large enough so that it contains the
    // computed S2LatLng of any contained point, not just the infinite-precision version.) We don't
    // need to expand longitude because longitude is calculated via a single call to atan2(), which
    // is guaranteed to be semi-monotonic. (In fact the Gnu implementation is also correctly
    // rounded, but we don't even need that here.)
    return bound.expanded(S2LatLng.fromRadians(DBL_EPSILON, 0));
  }

  @Override
  public boolean mayIntersect(S2Cell cell) {
    return cellId.intersects(cell.cellId);
  }

  /**
   * Returns true if the cell contains the given point "p". Note that unlike S2Loop/S2Polygon,
   * S2Cells are considered to be closed sets. This means that points along an S2Cell edge (or at
   * a vertex) belong to the adjacent cell(s) as well.
   *
   * <p>If instead you want every point to be contained by exactly one S2Cell, you will need to
   * convert the S2Cells to S2Loops, which implement point containment this way.
   *
   * <p>The point "p" does not need to be normalized.
   */
  @Override
  @JsMethod(name = "containsPoint")
  public boolean contains(S2Point p) {
    // We can't just call XYZtoFaceUV, because for points that lie on the boundary between two faces
    // (i.e. u or v is +1/-1) we need to return true for both adjacent cells.
    //
    // We can get away with not checking the face if the point matches the face of the cell here
    // because, for the 4 faces adjacent to 'face', p will be projected outside the range of
    // ([-1,1]x[-1,1]) and thus can't intersect the cell bounds (except on the face boundary which
    // we want).
    //
    // For the face opposite 'face', the sign of the UV coordinates of P will be flipped so it will
    // automatically fall outside the cell boundary as no cells cross the origin.
    R2Vector uvPoint = S2Projections.faceXyzToUv(face, p);
    if (uvPoint == null) {
      return false;
    }

    // Expand the (u,v) bound to ensure that
    //
    //   S2Cell(S2CellId.fromPoint(p)).contains(p)
    //
    // is always true.  To do this, we need to account for the error when converting from (u,v)
    // coordinates to (s,t) coordinates. At least in the case of S2_QUADRATIC_PROJECTION, the total
    // error is at most DBL_EPSILON.
    return uvPoint.x() >= (uMin - DBL_EPSILON)
        && uvPoint.x() <= (uMax + DBL_EPSILON)
        && uvPoint.y() >= (vMin - DBL_EPSILON)
        && uvPoint.y() <= (vMax + DBL_EPSILON);
  }

  // The point 'p' does not need to be normalized.
  @Override
  public boolean contains(S2Cell cell) {
    return cellId.contains(cell.cellId);
  }

  @SuppressWarnings("AndroidJdkLibsChecker")
  private double vertexChordDist2(S2Point uvw, DoubleBinaryOperator reducer) {
    double d1 = chordDist2(uvw, uMin, vMin);
    double d2 = chordDist2(uvw, uMin, vMax);
    double d3 = chordDist2(uvw, uMax, vMin);
    double d4 = chordDist2(uvw, uMax, vMax);
    return reducer.applyAsDouble(d1, reducer.applyAsDouble(d2, reducer.applyAsDouble(d3, d4)));
  }

  /** Returns the squared chord distance from {@code uvw} to position {@code uv}. */
  private static double chordDist2(S2Point uvw, double u, double v) {
    return uvw.getDistance2(new S2Point(u, v, 1).normalize());
  }

  /**
   * Given a point {@code p} and either the lower or upper edge of the {@link S2Cell} (specified by
   * setting {@code vEnd} to false or true respectively), returns true if {@code p} is closer to the
   * interior of that edge than it is to either endpoint.
   */
  private boolean uEdgeIsClosest(S2Point p, boolean vEnd) {
    double v = vEnd ? vMax : vMin;
    // These are the normals to the planes that are perpendicular to the edge and pass through one
    // of its two endpoints.
    S2Point dir0 = new S2Point(v * v + 1, -uMin * v, -uMin);
    S2Point dir1 = new S2Point(v * v + 1, -uMax * v, -uMax);
    return p.dotProd(dir0) > 0 && p.dotProd(dir1) < 0;
  }

  /**
   * Given a point {@code p} and either the left or right edge of the {@link S2Cell} (specified by
   * setting {@code uEnd} to false or true respectively), returns true if {@code p} is closer to the
   * interior of that edge than it is to either endpoint.
   */
  private boolean vEdgeIsClosest(S2Point p, boolean uEnd) {
    double u = uEnd ? uMax : uMin;
    // See comments above.
    S2Point dir0 = new S2Point(-u * vMin, u * u + 1, -vMin);
    S2Point dir1 = new S2Point(-u * vMax, u * u + 1, -vMax);
    return p.dotProd(dir0) > 0 && p.dotProd(dir1) < 0;
  }

  /**
   * Given the dot product of a point P with the normal of a u- or v-edge at the given coordinate
   * value, returns the distance from P to that edge.
   */
  private static double edgeDistance(double dirIJ, double uv) {
    // Let P be the target point and let R be the closest point on the given edge AB. The desired
    // distance XR can be expressed as PR^2 = PQ^2 + QR^2 where Q is the point P projected onto the
    // plane through the great circle through AB. We can compute the distance PQ^2 perpendicular to
    // the plane from "dirIJ" (the dot product of the target point P with the edge normal) and the
    // squared length of the edge normal (1 + uv**2).
    double pq2 = (dirIJ * dirIJ) / (1 + uv * uv);

    // We can compute the distance QR as (1 - OQ) where O is the sphere origin, and we can compute
    // OQ^2 = 1 - PQ^2 using the Pythagorean theorem. (This calculation loses accuracy as the angle
    // approaches Pi/2.)
    double qr = 1 - sqrt(1 - pq2);
    return pq2 + qr * qr;
  }

  /**
   * Returns the distance from the given point to the cell. Returns zero if the point is inside the
   * cell.
   */
  @JsMethod(name = "getDistanceToPoint")
  public S1ChordAngle getDistance(S2Point targetXyz) {
    return S1ChordAngle.fromLength2(getDistanceInternal(targetXyz, true));
  }

  /** Returns the distance to the given cell. Returns zero if one cell contains the other. */
  public S1ChordAngle getDistance(S2Cell target) {
    // If the cells intersect, the distance is zero. We use the (u,v) ranges rather than
    // S2CellId.intersects() so that cells that share a partial edge or a corner are considered to
    // intersect.
    if (face == target.face
        && uMin <= target.uMax && uMax >= target.uMin
        && vMin <= target.vMax && vMax >= target.vMin) {
      return S1ChordAngle.ZERO;
    }

    // Otherwise, the minimum distance always occurs between a vertex of one cell and an edge of the
    // other cell (including the edge endpoints). This represents a total of 32 possible
    // (vertex, edge) pairs.
    //
    // TODO(user): This could be optimized to be at least 5x faster by pruning the set of
    // possible closest vertex/edge pairs using the faces and (u,v) ranges of both cells.
    S2Point[] va = new S2Point[4];
    S2Point[] vb = new S2Point[4];
    for (int i = 0; i < 4; i++) {
      va[i] = getVertex(i);
      vb[i] = target.getVertex(i);
    }
    S1ChordAngle minDist = S1ChordAngle.INFINITY;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        minDist = S2EdgeUtil.updateMinDistance(va[i], vb[j], vb[(j + 1) & 3], minDist);
        minDist = S2EdgeUtil.updateMinDistance(vb[i], va[j], va[(j + 1) & 3], minDist);
      }
    }
    return minDist;
  }

  /** Returns the chord distance to targetXyz, with interior distance 0 iff toInterior is true. */
  private double getDistanceInternal(S2Point targetXyz, boolean toInterior) {
    // All calculations are done in the (u,v,w) coordinates of this cell's face.
    S2Point targetUvw = S2Projections.faceXyzToUvw(face, targetXyz);

    // Compute dot products with all four upward or rightward-facing edge normals. "dirIJ" is the
    // dot product for the edge corresponding to axis I, endpoint J. For example, dir01 is the right
    // edge of the S2Cell (corresponding to the upper endpoint of the u-axis).
    double dir00 = targetUvw.x - targetUvw.z * uMin;
    double dir01 = targetUvw.x - targetUvw.z * uMax;
    double dir10 = targetUvw.y - targetUvw.z * vMin;
    double dir11 = targetUvw.y - targetUvw.z * vMax;
    boolean inside = true;
    if (dir00 < 0) {
      inside = false; // Target is to the left of the cell.
      if (vEdgeIsClosest(targetUvw, false)) {
        return edgeDistance(-dir00, uMin);
      }
    }
    if (dir01 > 0) {
      inside = false; // Target is to the right of the cell.
      if (vEdgeIsClosest(targetUvw, true)) {
        return edgeDistance(dir01, uMax);
      }
    }
    if (dir10 < 0) {
      inside = false; // Target is below the cell.
      if (uEdgeIsClosest(targetUvw, false)) {
        return edgeDistance(-dir10, vMin);
      }
    }
    if (dir11 > 0) {
      inside = false; // Target is above the cell.
      if (uEdgeIsClosest(targetUvw, true)) {
        return edgeDistance(dir11, vMax);
      }
    }
    if (inside) {
      if (toInterior) {
        return 0;
      }

      // Although you might think of S2Cells as rectangles, they are actually arbitrary
      // quadrilaterals after they are projected onto the sphere. Therefore the simplest approach is
      // just to find the minimum distance to any of the four edges.
      return Doubles.min(
          edgeDistance(-dir00, uMin),
          edgeDistance(dir01, uMax),
          edgeDistance(-dir10, vMin),
          edgeDistance(dir11, vMax));
    }

    // Otherwise, the closest point is one of the four cell vertices. Note that it is *not* trivial
    // to narrow down the candidates based on the edge sign tests above, because (1) the edges don't
    // meet at right angles and (2) there are points on the far side of the sphere that are both
    // above *and* below the cell, etc.
    return vertexChordDist2(targetUvw, Math::min);
  }

  /** Returns the maximum distance from the cell (including its interior) to the given point. */
  @JsMethod(name = "getMaxDistanceToPoint")
  public S1ChordAngle getMaxDistance(S2Point target) {
    // First check the 4 cell vertices. If all are within the hemisphere centered around target, the
    // max distance will be to one of these vertices.
    S2Point targetUvw = S2Projections.faceXyzToUvw(face, target);
    double maxDist = vertexChordDist2(targetUvw, Math::max);
    if (maxDist <= S1ChordAngle.RIGHT.getLength2()) {
      return S1ChordAngle.fromLength2(maxDist);
    }

    // Otherwise, find the minimum distance dMin to the antipodal point and the maximum distance
    // will be Pi - dMin.
    return S1ChordAngle.sub(S1ChordAngle.STRAIGHT, getDistance(target.neg()));
  }

  /** Returns the maximum distance from the cell (including its interior) to the given edge AB. */
  @JsMethod(name = "getMaxDistanceToEdge")
  public S1ChordAngle getMaxDistance(S2Point a, S2Point b) {
    // If the maximum distance from both endpoints to the cell is less than Pi/2 then the maximum
    // distance from the edge to the cell is the maximum of the two endpoint distances.
    S1ChordAngle da = getMaxDistance(a);
    S1ChordAngle db = getMaxDistance(b);
    S1ChordAngle maxDist = da.compareTo(db) < 0 ? db : da;
    // TODO(user): Use a thresholded distance via S2Predicates.
    if (maxDist.compareTo(S1ChordAngle.RIGHT) <= 0) {
      return maxDist;
    } else {
      return S1ChordAngle.sub(S1ChordAngle.STRAIGHT, getDistanceToEdge(a.neg(), b.neg()));
    }
  }

  /** Returns the maximum distance from the cell, including interior, to the given target cell. */
  public S1ChordAngle getMaxDistance(S2Cell target) {
    // If the antipodal target intersects the cell, the distance is S1ChordAngle.STRAIGHT.
    // The antipodal UV is the transpose of the original UV, interpreted within the opposite face.
    if (face == (target.face >= 3 ? target.face - 3 : target.face + 3)
        && uMin <= target.vMax && uMax >= target.vMin
        && vMin <= target.uMax && vMax >= target.uMin) {
      return S1ChordAngle.STRAIGHT;
    }

    // Otherwise, the maximum distance always occurs between a vertex of one cell and an edge of the
    // other cell (including the edge endpoints). This represents a total of 32 possible
    // (vertex, edge) pairs.
    //
    // TODO(user): When the maximum distance is at most Pi/2, the maximum is always attained
    // between a pair of vertices, and this could be made much faster by testing each vertex pair
    // once rather than the current 4 times.
    S2Point[] va = new S2Point[4];
    S2Point[] vb = new S2Point[4];
    for (int i = 0; i < 4; i++) {
      va[i] = getVertex(i);
      vb[i] = target.getVertex(i);
    }
    S1ChordAngle maxDist = S1ChordAngle.NEGATIVE;
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        maxDist = S2EdgeUtil.updateMaxDistance(va[i], vb[j], vb[(j + 1) & 3], maxDist);
        maxDist = S2EdgeUtil.updateMaxDistance(vb[i], va[j], va[(j + 1) & 3], maxDist);
      }
    }
    return maxDist;
  }

  /**
   * Returns the minimum distance from the cell to the given edge AB, or zero if the edge intersects
   * the cell interior.
   */
  public S1ChordAngle getDistanceToEdge(S2Point a, S2Point b) {
    // TODO(user): Possible optimizations:
    //  - Currently the (cell vertex, edge endpoint) distances are computed twice each, and the
    //    length of AB is computed 4 times.
    //  - To fix this, refactor GetDistance(target) so that it skips calculating the distance to
    //    each cell vertex. Instead, compute the cell vertices and distances in this function, and
    //    add a low-level getMinDistance that allows the XA, XB, and AB distances to be passed
    //    in.
    //  - It might also be more efficient to do all calculations in UVW-space, since this would
    //    involve transforming 2 points rather than 4.

    // If a and b are equal, just return the distance to the point.
    if (a.equalsPoint(b)) {
      return getDistance(a);
    }

    // First, check the minimum distance to the edge endpoints A and B. (This also detects whether
    // either endpoint is inside the cell.)
    S1ChordAngle minDist = Ordering.natural().min(getDistance(a), getDistance(b));
    if (minDist.isZero()) {
      return minDist;
    }

    // Otherwise, check whether the edge crosses the cell boundary.
    S2Point[] v = new S2Point[4];
    for (int i = 0; i < 4; i++) {
      v[i] = getVertex(i);
    }
    S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(a, b, v[3]);
    for (int i = 0; i < 4; i++) {
      if (crosser.robustCrossing(v[i]) >= 0) {
        return S1ChordAngle.ZERO;
      }
    }

    // Finally, check whether the minimum distance occurs between a cell vertex and the interior of
    // the edge AB. (Some of this work is redundant, since it also checks the distance to the
    // endpoints A and B again.)
    //
    // Note that we don't need to check the distance from the interior of AB to the interior of a
    // cell edge, because the only way that this distance can be minimal is if the two edges cross
    // (already checked above).
    for (int i = 0; i < 4; i++) {
      minDist = S2EdgeUtil.updateMinDistance(v[i], a, b, minDist);
    }

    return minDist;
  }

  /** Returns the distance from the cell boundary to the given point. */
  public S1ChordAngle getBoundaryDistance(S2Point target) {
    return S1ChordAngle.fromLength2(getDistanceInternal(target, false));
  }

  private void init(S2CellId id) {
    // Set cell properties from the ID and the FaceIJ of the ID.
    cellId = id;
    face = (byte) id.face();
    long ijo = id.toIJOrientation();
    orientation = (byte) S2CellId.getOrientation(ijo);
    level = (byte) id.level();
    int i = S2CellId.getI(ijo);
    int j = S2CellId.getJ(ijo);
    int cellSize = id.getSizeIJ();
    this.uMin = S2Projections.PROJ.ijToUV(i, cellSize);
    this.uMax = S2Projections.PROJ.ijToUV(i + cellSize, cellSize);
    this.vMin = S2Projections.PROJ.ijToUV(j, cellSize);
    this.vMax = S2Projections.PROJ.ijToUV(j + cellSize, cellSize);
  }

  private S2Point getPoint(int i, int j) {
    return S2Projections.faceUvToXyz(face, i == 0 ? uMin : uMax, j == 0 ? vMin : vMax);
  }

  @Override
  public String toString() {
    return "[" + face + ", " + level + ", " + orientation + ", " + cellId + "]";
  }

  @Override
  public int hashCode() {
    int value = 17;
    value = 37 * (37 * (37 * value + face) + orientation) + level;
    return 37 * value + id().hashCode();
  }

  @Override
  public boolean equals(Object that) {
    if (that instanceof S2Cell) {
      S2Cell thatCell = (S2Cell) that;
      return this.face == thatCell.face
          && this.level == thatCell.level
          && this.orientation == thatCell.orientation
          && this.cellId.equals(thatCell.cellId);
    }
    return false;
  }

  /* A double function of two double parameters. */
  // TODO(user): Remove this once the android JDK has this interface.
  private interface DoubleBinaryOperator {
    /** Returns the result of this function applied to {@code a} and {@code b}. */
    double applyAsDouble(double a, double b);
  }
}
