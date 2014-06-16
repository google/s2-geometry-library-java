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
import com.google.common.base.Preconditions;

import java.io.Serializable;

/**
 * Base class for methods shared between the immutable {@link S2LatLngRect} and the mutable
 * {@link S2LatLngRect.Builder}.
 */
@GwtCompatible(serializable = true)
public abstract strictfp class S2LatLngRectBase implements S2Region, Serializable {
  protected final R1Interval lat;
  protected final S1Interval lng;

  /**
   * Constructs a rectangle from minimum and maximum latitudes and longitudes. If 
   * lo.lng() > hi.lng(), the rectangle spans the 180 degree longitude line. Both points must be
   * normalized, with lo.lat() <= hi.lat().  The rectangle contains all the points p such that
   * 'lo' <= p <= 'hi', where '<=' is defined in the obvious way.
   */
  S2LatLngRectBase(final S2LatLng lo, final S2LatLng hi) {
    lat = new R1Interval(lo.lat().radians(), hi.lat().radians());
    lng = new S1Interval(lo.lng().radians(), hi.lng().radians());
    // assert (isValid());
  }

  /**
   * Constructs a rectangle from latitude and longitude intervals. The two intervals must either
   * be both empty or both non-empty, and the latitude interval must not extend outside [-90, +90]
   * degrees.  Note that both intervals (and hence the rectangle) are closed.
   */
  S2LatLngRectBase(R1Interval lat, S1Interval lng) {
    this.lat = lat;
    this.lng = lng;
    // assert (isValid());
  }
  
  /**
   * Constructs a rectangle with lat and lng fields set to empty intervals, as defined in
   * {@link R1Interval} and {@link S1Interval}.
   */
  S2LatLngRectBase() {
    lat = R1Interval.empty();
    lng = S1Interval.empty();
  }

  /**
   * Returns true if the rectangle is valid, which essentially just means that the latitude bounds
   * do not exceed Pi/2 in absolute value and the longitude bounds do not exceed Pi in absolute
   * value.  Also, if either the latitude or longitude bound is empty then both must be.
   */
  public final boolean isValid() {
    // The lat/lng ranges must either be both empty or both non-empty.
    return (Math.abs(lat.lo()) <= S2.M_PI_2 && Math.abs(lat.hi()) <= S2.M_PI_2
        && lng.isValid() && lat.isEmpty() == lng.isEmpty());
  }

  // Accessor methods.
  public final S1Angle latLo() {
    return S1Angle.radians(lat.lo());
  }

  public final S1Angle latHi() {
    return S1Angle.radians(lat.hi());
  }

  public final S1Angle lngLo() {
    return S1Angle.radians(lng.lo());
  }

  public final S1Angle lngHi() {
    return S1Angle.radians(lng.hi());
  }

  /** Returns the latitude range of this rectangle. */
  public abstract R1Interval lat();

  /** Returns the longitude range of this rectangle. */
  public abstract S1Interval lng();

  public final S2LatLng lo() {
    return new S2LatLng(latLo(), lngLo());
  }

  public final S2LatLng hi() {
    return new S2LatLng(latHi(), lngHi());
  }

  /** Returns true if the rectangle is empty, i.e. it contains no points at all. */
  public final boolean isEmpty() {
    return lat.isEmpty();
  }

  /** Returns true if the rectangle is full, i.e. it contains all points. */
  public final boolean isFull() {
    return lat.equals(S2LatLngRect.fullLat()) && lng.isFull();
  }

  /**
   * Returns true if lng_.lo() > lng_.hi(), i.e. the rectangle crosses the 180 degree latitude line.
   */
  public final boolean isInverted() {
    return lng.isInverted();
  }

  /**
   *  Returns the k<super>th</super> vertex of the rectangle (k = 0,1,2,3) in CCW order (lower-left,
   * lower right, upper right, upper left).
   */
  public final S2LatLng getVertex(int k) {
    switch (k) {
      case 0:
        return S2LatLng.fromRadians(lat.lo(), lng.lo());
      case 1:
        return S2LatLng.fromRadians(lat.lo(), lng.hi());
      case 2:
        return S2LatLng.fromRadians(lat.hi(), lng.hi());
      case 3:
        return S2LatLng.fromRadians(lat.hi(), lng.lo());
      default:
        throw new IllegalArgumentException("Invalid vertex index.");
    }
  }

  /**
   * Returns the center of the rectangle in latitude-longitude space (in general this is not the
   * center of the region on the sphere).
   */
  public final S2LatLng getCenter() {
    return S2LatLng.fromRadians(lat.getCenter(), lng.getCenter());
  }

  /**
   * Returns the minimum distance (measured along the surface of the sphere) from a given point to 
   * the rectangle (both its boundary and its interior). The latLng must be valid.
   */
  public final S1Angle getDistance(S2LatLng p) {
    // The algorithm here is the same as in getDistance(S2LagLngRect), only with simplified
    // calculations.
    S2LatLngRectBase a = this;

    Preconditions.checkState(!a.isEmpty());
    Preconditions.checkArgument(p.isValid());

    if (a.lng().contains(p.lng().radians())) {
      return S1Angle.radians(Math.max(0.0, Math.max(p.lat().radians() - a.lat().hi(),
          a.lat().lo() - p.lat().radians())));
    }

    S1Interval interval = new S1Interval(a.lng().hi(), a.lng().complement().getCenter());
    double aLng = a.lng().lo();
    if (interval.contains(p.lng().radians())) {
      aLng = a.lng().hi();
    }

    S2Point lo = S2LatLng.fromRadians(a.lat().lo(), aLng).toPoint();
    S2Point hi = S2LatLng.fromRadians(a.lat().hi(), aLng).toPoint();
    S2Point loCrossHi =
        S2LatLng.fromRadians(0, aLng - S2.M_PI_2).normalized().toPoint();
    return S2EdgeUtil.getDistance(p.toPoint(), lo, hi, loCrossHi);
  }

  /** 
   * Returns the minimum distance (measured along the surface of the sphere) to the given
   * S2LatLngRectBase. Both S2LatLngRectBases must be non-empty.
   */
  public final S1Angle getDistance(S2LatLngRectBase other) {
    S2LatLngRectBase a = this;
    S2LatLngRectBase b = other;

    Preconditions.checkState(!a.isEmpty());
    Preconditions.checkArgument(!b.isEmpty());

    // First, handle the trivial cases where the longitude intervals overlap.
    if (a.lng().intersects(b.lng())) {
      if (a.lat().intersects(b.lat())) {
        // Intersection between a and b.
        return S1Angle.radians(0);
      }

      // We found an overlap in the longitude interval, but not in the latitude interval. This means
      // the shortest path travels along some line of longitude connecting the high-latitude of the 
      // lower rect with the low-latitude of the higher rect.
      S1Angle lo, hi;
      if (a.lat().lo() > b.lat().hi()) {
        lo = b.latHi();
        hi = a.latLo();
      } else {
        lo = a.latHi();
        hi = b.latLo();
      }
      return S1Angle.radians(hi.radians() - lo.radians());
    }

    // The longitude intervals don't overlap. In this case, the closest points occur somewhere on
    // the pair of longitudinal edges which are nearest in longitude-space.
    S1Angle aLng, bLng;
    S1Interval loHi = S1Interval.fromPointPair(a.lng().lo(), b.lng().hi());
    S1Interval hiLo = S1Interval.fromPointPair(a.lng().hi(), b.lng().lo());
    if (loHi.getLength() < hiLo.getLength()) {
      aLng = a.lngLo();
      bLng = b.lngHi();
    } else {
      aLng = a.lngHi();
      bLng = b.lngLo();
    }

    // The shortest distance between the two longitudinal segments will include at least one segment
    // endpoint. We could probably narrow this down further to a single point-edge distance by 
    // comparing the relative latitudes of the endpoints, but for the sake of clarity, we'll do all
    // four point-edge distance tests.
    S2Point aLo = new S2LatLng(a.latLo(), aLng).toPoint();
    S2Point aHi = new S2LatLng(a.latHi(), aLng).toPoint();
    S2Point aLoCrossHi =
        S2LatLng.fromRadians(0, aLng.radians() - S2.M_PI_2).normalized().toPoint();
    S2Point bLo = new S2LatLng(b.latLo(), bLng).toPoint();
    S2Point bHi = new S2LatLng(b.latHi(), bLng).toPoint();
    S2Point bLoCrossHi =
        S2LatLng.fromRadians(0, bLng.radians() - S2.M_PI_2).normalized().toPoint();

    return S1Angle.min(S2EdgeUtil.getDistance(aLo, bLo, bHi, bLoCrossHi),
        S1Angle.min(S2EdgeUtil.getDistance(aHi, bLo, bHi, bLoCrossHi),
            S1Angle.min(S2EdgeUtil.getDistance(bLo, aLo, aHi, aLoCrossHi),
                S2EdgeUtil.getDistance(bHi, aLo, aHi, aLoCrossHi))));
  }

  /**
   * Returns the width and height of this rectangle in latitude-longitude space. Empty rectangles
   * have a negative width and height.
   */
  public final S2LatLng getSize() {
    return S2LatLng.fromRadians(lat.getLength(), lng.getLength());
  }

  /** More efficient version of contains() that accepts a S2LatLng rather than an S2Point. */
  public final boolean contains(S2LatLng ll) {
    // assert (ll.isValid());
    return (lat.contains(ll.latRadians()) && lng.contains(ll.lngRadians()));
  }

  /**
   * Returns true if and only if the given point is contained in the interior of the region (i.e.
   * the region excluding its boundary). The point 'p' does not need to be normalized.
   */
  public final boolean interiorContains(S2Point p) {
    return interiorContains(new S2LatLng(p));
  }

  /**
   * More efficient version of interiorContains() that accepts a S2LatLng rather than an S2Point.
   */
  public final boolean interiorContains(S2LatLng ll) {
    // assert (ll.isValid());
    return (lat.interiorContains(ll.lat().radians()) && lng
        .interiorContains(ll.lng().radians()));
  }

  /** Returns true if and only if the rectangle contains the given other rectangle. */
  public final boolean contains(S2LatLngRectBase other) {
    return lat.contains(other.lat) && lng.contains(other.lng);
  }

  /**
   * Returns true if and only if the interior of this rectangle contains all points of the given
   * other rectangle (including its boundary).
   */
  public final boolean interiorContains(S2LatLngRectBase other) {
    return (lat.interiorContains(other.lat) && lng
        .interiorContains(other.lng));
  }

  /** Returns true if this rectangle and the given other rectangle have any points in common. */
  public final boolean intersects(S2LatLngRectBase other) {
    return lat.intersects(other.lat) && lng.intersects(other.lng);
  }

  /**
   * Returns true if this rectangle intersects the given cell. (This is an exact
   * test and may be fairly expensive, see also MayIntersect below.)
   */
  public final boolean intersects(S2Cell cell) {
    // First we eliminate the cases where one region completely contains the other. Once these are 
    // disposed of, then the regions will intersect if and only if their boundaries intersect.

    if (isEmpty()) {
      return false;
    }
    if (contains(cell.getCenterRaw())) {
      return true;
    }
    if (cell.contains(getCenter().toPoint())) {
      return true;
    }

    // Quick rejection test (not required for correctness).
    if (!intersects(cell.getRectBound())) {
      return false;
    }

    // Now check whether the boundaries intersect. Unfortunately, a latitude-longitude rectangle 
    // does not have straight edges -- two edges are curved, and at least one of them is concave.

    // Precompute the cell vertices as points and latitude-longitudes.
    S2Point[] cellV = new S2Point[4];
    S2LatLng[] cellLl = new S2LatLng[4];
    for (int i = 0; i < 4; ++i) {
      // Must be normalized.
      cellV[i] = cell.getVertex(i);
      cellLl[i] = new S2LatLng(cellV[i]);
      if (contains(cellLl[i])) {
        // Quick acceptance test.
        return true;
      }
    }

    for (int i = 0; i < 4; ++i) {
      S1Interval edgeLng = S1Interval.fromPointPair(
          cellLl[i].lng().radians(), cellLl[(i + 1) & 3].lng().radians());
      if (!lng.intersects(edgeLng)) {
        continue;
      }

      final S2Point a = cellV[i];
      final S2Point b = cellV[(i + 1) & 3];
      if (edgeLng.contains(lng.lo())) {
        if (intersectsLngEdge(a, b, lat, lng.lo())) {
          return true;
        }
      }
      if (edgeLng.contains(lng.hi())) {
        if (intersectsLngEdge(a, b, lat, lng.hi())) {
          return true;
        }
      }
      if (intersectsLatEdge(a, b, lat.lo(), lng)) {
        return true;
      }
      if (intersectsLatEdge(a, b, lat.hi(), lng)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Returns true if and only if the interior of this rectangle intersects any point (including the
   * boundary) of the given other rectangle.
   */
  public final boolean interiorIntersects(S2LatLngRectBase other) {
    return (lat.interiorIntersects(other.lat) && lng.interiorIntersects(other.lng));
  }

  /** Returns the surface area of this rectangle on the unit sphere. */
  public final double area() {
    if (isEmpty()) {
      return 0;
    }

    // This is the size difference of the two spherical caps, multiplied by the longitude ratio.
    return lng().getLength() * Math.abs(Math.sin(latHi().radians()) - Math.sin(latLo().radians()));
  }

  /** Returns true if these are the same type of rectangle and contain the same set of points. */
  @Override
  public final boolean equals(Object that) {
    if ((that == null) || this.getClass() != that.getClass()) {
      return false;
    }
    S2LatLngRectBase otherRect = (S2LatLngRectBase) that;
    return lat().equals(otherRect.lat()) && lng().equals(otherRect.lng());
  }

  /**
   * Returns true if the latitude and longitude intervals of the two rectangles are the same up to
   * the given tolerance.  See {@link R1Interval} and {@link S1Interval} for details.
   */
  public final boolean approxEquals(S2LatLngRectBase other, double maxError) {
    return lat.approxEquals(other.lat, maxError) && lng.approxEquals(other.lng, maxError);
  }

  /** Returns true if this rectangle is very nearly identical to the given other rectangle. */
  public final boolean approxEquals(S2LatLngRectBase other) {
    return approxEquals(other, 1e-15);
  }

  /**
   * As {@link #approxEquals(S2LatLngRectBase, double)}, but with separate tolerances for latitude
   * and longitude.
   */
  public final boolean approxEquals(S2LatLngRectBase other, S2LatLng maxError) {
    return lat.approxEquals(other.lat, maxError.lat().radians())
        && lng.approxEquals(other.lng, maxError.lng().radians());
  }

  @Override
  public final int hashCode() {
    int value = 17;
    value = 37 * value + lat.hashCode();
    return (37 * value + lng.hashCode());
  }

  // //////////////////////////////////////////////////////////////////////
  // S2Region interface (see {@code S2Region} for details):

  @Override
  public S2Cap getCapBound() {
    // We consider two possible bounding caps, one whose axis passes through the center of the
    // lat-lng rectangle and one whose axis is the north or south pole. We return the smaller of the
    // two caps.
    if (isEmpty()) {
      return S2Cap.empty();
    }

    double poleZ, poleAngle;
    if (lat.lo() + lat.hi() < 0) {
      // South pole axis yields smaller cap.
      poleZ = -1;
      poleAngle = S2.M_PI_2 + lat.hi();
    } else {
      poleZ = 1;
      poleAngle = S2.M_PI_2 - lat.lo();
    }
    S2Cap poleCap = S2Cap.fromAxisAngle(new S2Point(0, 0, poleZ), S1Angle.radians(poleAngle));

    // For bounding rectangles that span 180 degrees or less in longitude, the maximum cap size is
    // achieved at one of the rectangle vertices. For rectangles that are larger than 180 degrees,
    // we punt and always return a bounding cap centered at one of the two poles.
    double lngSpan = lng.hi() - lng.lo();
    if (Platform.IEEEremainder(lngSpan, 2 * S2.M_PI) >= 0) {
      if (lngSpan < 2 * S2.M_PI) {
        S2Cap midCap = S2Cap.fromAxisAngle(getCenter().toPoint(), S1Angle.radians(0));
        for (int k = 0; k < 4; ++k) {
          midCap = midCap.addPoint(getVertex(k).toPoint());
        }
        if (midCap.height() < poleCap.height()) {
          return midCap;
        }
      }
    }
    return poleCap;
  }

  /**
   * Returns true if this latitude/longitude region contains the given cell. A latitude-longitude
   * rectangle contains a cell if and only if it contains the cell's bounding rectangle, making this
   * an exact test. Note, however, that the cell must be valid; an error may result if e.g.
   * {@link S2CellId#sentinel()} is passed here.
   */
  @Override
  public final boolean contains(S2Cell cell) {
    return contains(cell.getRectBound());
  }

  /**
   * This test is cheap but is NOT exact. Use Intersects() if you want a more accurate and more
   * expensive test. Note that when this method is used by an S2RegionCoverer, the accuracy isn't
   * all that important since if a cell may intersect the region then it is subdivided, and the
   * accuracy of this method goes up as the cells get smaller.
   */
  @Override
  public final boolean mayIntersect(S2Cell cell) {
    // This test is cheap but is NOT exact (see s2latlngrect.h).
    return intersects(cell.getRectBound());
  }

  /** The point 'p' does not need to be normalized. */
  public final boolean contains(S2Point p) {
    return contains(new S2LatLng(p));
  }

  /** Returns true if the edge AB intersects the given edge of constant longitude. */
  public static final boolean intersectsLngEdge(S2Point a, S2Point b, R1Interval lat, double lng) {
    // Return true if the segment AB intersects the given edge of constant longitude. The nice thing
    // about edges of constant longitude is that they are straight lines on the sphere (geodesics).

    return S2.simpleCrossing(a, b, S2LatLng.fromRadians(lat.lo(), lng)
        .toPoint(), S2LatLng.fromRadians(lat.hi(), lng).toPoint());
  }

  /** Returns true if the edge AB intersects the given edge of constant latitude. */
  public static final boolean intersectsLatEdge(S2Point a, S2Point b, double lat, S1Interval lng) {
    // Return true if the segment AB intersects the given edge of constant latitude. Unfortunately,
    // lines of constant latitude are curves on the sphere. They can intersect a straight edge in
    // 0, 1, or 2 points.  
    // assert (S2.isUnitLength(a) && S2.isUnitLength(b));

    // First, compute the normal to the plane AB that points vaguely north.
    S2Point z = S2Point.normalize(S2.robustCrossProd(a, b));
    if (z.z < 0) {
      z = S2Point.neg(z);
    }

    // Extend this to an orthonormal frame (x,y,z) where x is the direction where the great circle
    // through AB achieves its maximum latitude.
    S2Point y = S2Point.normalize(S2.robustCrossProd(z, S2Point.Z_POS));
    S2Point x = S2Point.crossProd(y, z);
    // assert (S2.isUnitLength(x) && x.z >= 0);

    // Compute the angle "theta" from the x-axis (in the x-y plane defined above) where the great
    // circle intersects the given line of latitude.
    double sinLat = Math.sin(lat);
    if (Math.abs(sinLat) >= x.z) {
      return false; // The great circle does not reach the given latitude.
    }
    // assert (x.z > 0);
    double cosTheta = sinLat / x.z;
    double sinTheta = Math.sqrt(1 - cosTheta * cosTheta);
    double theta = Math.atan2(sinTheta, cosTheta);

    // The candidate intersection points are located +/- theta in the x-y plane. For an intersection
    // to be valid, we need to check that the intersection point is contained in the interior of the
    // edge AB and also that it is contained within the given longitude interval "lng".

    // Compute the range of theta values spanned by the edge AB.
    S1Interval abTheta = S1Interval.fromPointPair(Math.atan2(
        a.dotProd(y), a.dotProd(x)), Math.atan2(b.dotProd(y), b.dotProd(x)));

    if (abTheta.contains(theta)) {
      // Check if the intersection point is also in the given "lng" interval.
      S2Point isect = S2Point.add(S2Point.mul(x, cosTheta), S2Point.mul(y, sinTheta));
      if (lng.contains(Math.atan2(isect.y, isect.x))) {
        return true;
      }
    }
    if (abTheta.contains(-theta)) {
      // Check if the intersection point is also in the given "lng" interval.
      S2Point intersection = S2Point.sub(S2Point.mul(x, cosTheta), S2Point.mul(y, sinTheta));
      if (lng.contains(Math.atan2(intersection.y, intersection.x))) {
        return true;
      }
    }

    return false;
  }

  @Override
  public final String toString() {
    return "[Lo=" + lo() + ", Hi=" + hi() + "]";
  }

  public final String toStringDegrees() {
    return "[Lo=" + lo().toStringDegrees() + ", Hi=" + hi().toStringDegrees() + "]";
  }
}
