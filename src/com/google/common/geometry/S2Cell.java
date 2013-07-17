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

import static com.google.common.geometry.S2Projections.PROJ;

import com.google.common.annotations.GwtCompatible;
import com.google.common.geometry.S2CellId.FaceIJ;

import java.io.Serializable;

/**
 * An S2Cell is an S2Region object that represents a cell. Unlike S2CellIds, it
 * supports efficient containment and intersection tests. However, it is also a
 * more expensive representation.
 *
 */
@GwtCompatible(serializable = true)
public final strictfp class S2Cell implements S2Region, Serializable {
  byte face;
  byte level;
  byte orientation;
  S2CellId cellId;
  double uMin;
  double uMax;
  double vMin;
  double vMax;

  /**
   * Default constructor used only internally.
   */
  S2Cell() {
  }

  /**
   * An S2Cell always corresponds to a particular S2CellId. The other
   * constructors are just convenience methods.
   */
  public S2Cell(S2CellId id) {
    init(id);
  }

  /** Returns the cell corresponding to the given S2 cube face. */
  public static S2Cell fromFace(int face) {
    return new S2Cell(S2CellId.fromFace(face));
  }

  /**
   * Returns a cell given its face (range 0..5), Hilbert curve position within that face (an
   * unsigned integer with {@link S2CellId#POS_BITS} bits), and level (range 0..kMaxLevel).  The
   * given position will be modified to correspond to the Hilbert curve position at the center of
   * the returned cell.  This is a static function rather than a constructor in order to indicate
   * what the arguments represent.
   */
  public static S2Cell fromFacePosLevel(int face, long pos, int level) {
    return new S2Cell(S2CellId.fromFacePosLevel(face, pos, level));
  }

  // Convenience methods.
  public S2Cell(S2Point p) {
    init(S2CellId.fromPoint(p));
  }

  public S2Cell(S2LatLng ll) {
    init(S2CellId.fromLatLng(ll));
  }

  public S2CellId id() {
    return cellId;
  }

  public int face() {
    return face;
  }

  public byte level() {
    return level;
  }

  public byte orientation() {
    return orientation;
  }

  /** Returns true if this cell is a leaf-cell, i.e. it has no children. */
  public boolean isLeaf() {
    return level == S2CellId.MAX_LEVEL;
  }

  /** As {@link #getVertexRaw(int)}, except the point is normalized to unit length. */
  public S2Point getVertex(int k) {
    return S2Point.normalize(getVertexRaw(k));
  }

  /**
   * Returns the k<super>th</super> vertex of the cell (k = 0,1,2,3).  Vertices are returned in CCW
   * order (lower left, lower right, upper right, upper left in the UV plane). The points are not
   * necessarily unit length.
   */
  public S2Point getVertexRaw(int k) {
    // Vertices are returned in the order SW, SE, NE, NW.
    return S2Projections.faceUvToXyz(face,
        ((k >> 1) ^ (k & 1)) == 0 ? uMin : uMax,
        (k >> 1) == 0 ? vMin : vMax);
  }

  /** As {@link #getEdgeRaw(int)}, except the point is normalized to unit length. */
  public S2Point getEdge(int k) {
    return S2Point.normalize(getEdgeRaw(k));
  }

  /**
   * Returns the inward-facing normal of the great circle passing through the edge from vertex k to
   * vertex k+1 (mod 4).  The normals returned by getEdgeRaw are not necessarily unit length.
   */
  public S2Point getEdgeRaw(int k) {
    switch (k) {
      case 0:
        return S2Projections.getVNorm(face, vMin); // Bottom
      case 1:
        return S2Projections.getUNorm(face, uMax); // Right
      case 2:
        return S2Point.neg(S2Projections.getVNorm(face, vMax)); // Top
      default:
        return S2Point.neg(S2Projections.getUNorm(face, uMin)); // Left
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
   * for (pos=0, id=childBegin(); !id.equals(childEnd()); id = id.next(), ++pos) {
   *   children[i].init(id);
   * }
   * </pre>
   *
   * <p>except that it is more than two times faster.
   */
  public boolean subdivide(S2Cell children[]) {
    // This function is equivalent to just iterating over the child cell ids
    // and calling the S2Cell constructor, but it is about 2.5 times faster.

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
      child.orientation = (byte) (orientation ^ S2.posToOrientation(pos));
      child.cellId = id;
      // We want to split the cell in half in "u" and "v".  To decide which
      // side to set equal to the midpoint value, we look at cell's (i,j)
      // position within its parent.  The index for "i" is in bit 1 of ij.
      int ij = S2.posToIJ(orientation, pos);
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
   * Return the direction vector corresponding to the center in (s,t)-space of
   * the given cell. This is the point at which the cell is divided into four
   * subcells; it is not necessarily the centroid of the cell in (u,v)-space or
   * (x,y,z)-space. The point returned by GetCenterRaw is not necessarily unit
   * length.
   */
  public S2Point getCenter() {
    return S2Point.normalize(getCenterRaw());
  }

  public S2Point getCenterRaw() {
    return cellId.toPointRaw();
  }

  /** Returns the bounds of this cell in (u,v)-space. */
  public R2Rect getBoundUV() {
    return new R2Rect(new R1Interval(uMin, uMax), new R1Interval(vMin, vMax));
  }

  /**
   * Return the center of the cell in (u,v) coordinates (see {@code
   * S2Projections}). Note that the center of the cell is defined as the point
   * at which it is recursively subdivided into four children; in general, it is
   * not at the midpoint of the (u,v) rectangle covered by the cell
   */
  public R2Vector getCenterUV() {
    return cellId.getCenterUV();
  }

  /**
   * Return the average area for cells at the given level.
   */
  public static double averageArea(int level) {
    return PROJ.avgArea.getValue(level);
  }

  /**
   * Return the average area of cells at this level. This is accurate to within
   * a factor of 1.7 (for S2_QUADRATIC_PROJECTION) and is extremely cheap to
   * compute.
   */
  public double averageArea() {
    return averageArea(level);
  }

  /**
   * Return the approximate area of this cell. This method is accurate to within
   * 3% percent for all cell sizes and accurate to within 0.1% for cells at
   * level 5 or higher (i.e. 300km square or smaller). It is moderately cheap to
   * compute.
   */
  public double approxArea() {

    // All cells at the first two levels have the same area.
    if (level < 2) {
      return averageArea(level);
    }

    // First, compute the approximate area of the cell when projected
    // perpendicular to its normal. The cross product of its diagonals gives
    // the normal, and the length of the normal is twice the projected area.
    double flatArea = 0.5 * S2Point.crossProd(
        S2Point.sub(getVertex(2), getVertex(0)), S2Point.sub(getVertex(3), getVertex(1))).norm();

    // Now, compensate for the curvature of the cell surface by pretending
    // that the cell is shaped like a spherical cap. The ratio of the
    // area of a spherical cap to the area of its projected disc turns out
    // to be 2 / (1 + sqrt(1 - r*r)) where "r" is the radius of the disc.
    // For example, when r=0 the ratio is 1, and when r=1 the ratio is 2.
    // Here we set Pi*r*r == flat_area to find the equivalent disc.
    return flatArea * 2 / (1 + Math.sqrt(1 - Math.min(S2.M_1_PI * flatArea, 1.0)));
  }

  /**
   * Return the area of this cell as accurately as possible. This method is more
   * expensive but it is accurate to 6 digits of precision even for leaf cells
   * (whose area is approximately 1e-18).
   */
  public double exactArea() {
    S2Point v0 = getVertex(0);
    S2Point v1 = getVertex(1);
    S2Point v2 = getVertex(2);
    S2Point v3 = getVertex(3);
    return S2.area(v0, v1, v2) + S2.area(v0, v2, v3);
  }

  // //////////////////////////////////////////////////////////////////////
  // S2Region interface (see {@code S2Region} for details):

  // NOTE: This should be marked as @Override, but clone() isn't present in GWT's version of
  // Object, so we can't mark it as such.
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
    // Use the cell center in (u,v)-space as the cap axis. This vector is
    // very close to GetCenter() and faster to compute. Neither one of these
    // vectors yields the bounding cap with minimal surface area, but they
    // are both pretty close.
    //
    // It's possible to show that the two vertices that are furthest from
    // the (u,v)-origin never determine the maximum cap size (this is a
    // possible future optimization).

    S2Point center = S2Point.normalize(S2Projections.faceUvToXyz(face, getCenterUV()));
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
  private static final double POLE_MIN_LAT = Math.asin(Math.sqrt(1. / 3)) - 0.5 * S2.DBL_EPSILON;

  @Override
  public S2LatLngRect getRectBound() {
    if (level > 0) {
      // Except for cells at level 0, the latitude and longitude extremes are
      // attained at the vertices.  Furthermore, the latitude range is
      // determined by one pair of diagonally opposite vertices and the
      // longitude range is determined by the other pair.
      //
      // We first determine which corner (i,j) of the cell has the largest
      // absolute latitude.  To maximize latitude, we want to find the point in
      // the cell that has the largest absolute z-coordinate and the smallest
      // absolute x- and y-coordinates.  To do this we look at each coordinate
      // (u and v), and determine whether we want to minimize or maximize that
      // coordinate based on the axis direction and the cell's (u,v) quadrant.
      double u = uMin + uMax;
      double v = vMin + vMax;
      int i = (S2Projections.getUAxis(face).z == 0 ? (u < 0) : (u > 0)) ? 1 : 0;
      int j = (S2Projections.getVAxis(face).z == 0 ? (v < 0) : (v > 0)) ? 1 : 0;
      R1Interval lat = R1Interval.fromPointPair(
          S2LatLng.latitude(getPoint(i, j)).radians(),
          S2LatLng.latitude(getPoint(1 - i, 1 - j)).radians());
      S1Interval lng = S1Interval.fromPointPair(
          S2LatLng.longitude(getPoint(i, 1 - j)).radians(),
          S2LatLng.longitude(getPoint(1 - i, j)).radians());

      // We grow the bounds slightly to make sure that the bounding rectangle
      // contains S2LatLng(P) for any point P inside the loop L defined by the
      // four *normalized* vertices.  Note that normalization of a vector can
      // change its direction by up to 0.5 * DBL_EPSILON radians, and it is not
      // enough just to add Normalize() calls to the code above because the
      // latitude/longitude ranges are not necessarily determined by diagonally
      // opposite vertex pairs after normalization.
      //
      // We would like to bound the amount by which the latitude/longitude of a
      // contained point P can exceed the bounds computed above.  In the case of
      // longitude, the normalization error can change the direction of rounding
      // leading to a maximum difference in longitude of 2 * DBL_EPSILON.  In
      // the case of latitude, the normalization error can shift the latitude by
      // up to 0.5 * DBL_EPSILON and the other sources of error can cause the
      // two latitudes to differ by up to another 1.5 * DBL_EPSILON, which also
      // leads to a maximum difference of 2 * DBL_EPSILON.
      return new S2LatLngRect(lat, lng)
          .expanded(S2LatLng.fromRadians(2 * S2.DBL_EPSILON, 2 * S2.DBL_EPSILON))
          .polarClosure();
    }

    // The face centers are the +X, +Y, +Z, -X, -Y, -Z axes in that order.
    // assert (((face < 3) ? 1 : -1) == S2Projections.getNorm(face).get(face % 3));

    S2LatLngRect bound;
    switch (face) {
      case 0:
        bound = new S2LatLngRect(
            new R1Interval(-S2.M_PI_4, S2.M_PI_4), new S1Interval(-S2.M_PI_4, S2.M_PI_4));
        break;
      case 1:
        bound = new S2LatLngRect(
            new R1Interval(-S2.M_PI_4, S2.M_PI_4), new S1Interval(S2.M_PI_4, 3 * S2.M_PI_4));
        break;
      case 2:
        bound = new S2LatLngRect(new R1Interval(POLE_MIN_LAT, S2.M_PI_2), S1Interval.full());
        break;
      case 3:
        bound = new S2LatLngRect(
            new R1Interval(-S2.M_PI_4, S2.M_PI_4), new S1Interval(3 * S2.M_PI_4, -3 * S2.M_PI_4));
        break;
      case 4:
        bound = new S2LatLngRect(
            new R1Interval(-S2.M_PI_4, S2.M_PI_4), new S1Interval(-3 * S2.M_PI_4, -S2.M_PI_4));
        break;
      default:
        bound = new S2LatLngRect(new R1Interval(-S2.M_PI_2, -POLE_MIN_LAT), S1Interval.full());
        break;
    }

    // Finally, we expand the bound to account for the error when a point P is
    // converted to an S2LatLng to test for containment.  (The bound should be
    // large enough so that it contains the computed S2LatLng of any contained
    // point, not just the infinite-precision version.)  We don't need to expand
    // longitude because longitude is calculated via a single call to atan2(),
    // which is guaranteed to be semi-monotonic.  (In fact the Gnu implementation
    // is also correctly rounded, but we don't even need that here.)
    return bound.expanded(S2LatLng.fromRadians(S2.DBL_EPSILON, 0));
  }

  @Override
  public boolean mayIntersect(S2Cell cell) {
    return cellId.intersects(cell.cellId);
  }

  public boolean contains(S2Point p) {
    // We can't just call XYZtoFaceUV, because for points that lie on the
    // boundary between two faces (i.e. u or v is +1/-1) we need to return
    // true for both adjacent cells.
    R2Vector uvPoint = S2Projections.faceXyzToUv(face, p);
    if (uvPoint == null) {
      return false;
    }
    return (uvPoint.x() >= uMin && uvPoint.x() <= uMax
        && uvPoint.y() >= vMin && uvPoint.y() <= vMax);
  }

  // The point 'p' does not need to be normalized.
  @Override
  public boolean contains(S2Cell cell) {
    return cellId.contains(cell.cellId);
  }

  private void init(S2CellId id) {
    // Set cell properties from the ID and the FaceIJ of the ID.
    cellId = id;
    FaceIJ fij = id.toFaceIJOrientation();
    face = (byte) fij.face;
    orientation = (byte) fij.orientation;
    level = (byte) id.level();
    R2Rect uv = S2CellId.ijLevelToBoundUv(fij.i, fij.j, level);
    R1Interval u = uv.x();
    this.uMin = u.lo();
    this.uMax = u.hi();
    R1Interval v = uv.y();
    this.vMin = v.lo();
    this.vMax = v.hi();
  }

  private S2Point getPoint(int i, int j) {
    return S2Projections.faceUvToXyz(face,
        i == 0 ? uMin : uMax,
        j == 0 ? vMin : vMax);
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
      return this.face == thatCell.face && this.level == thatCell.level
          && this.orientation == thatCell.orientation && this.cellId.equals(thatCell.cellId);
    }
    return false;
  }
}
