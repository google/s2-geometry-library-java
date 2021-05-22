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
 * Base class for methods shared between the immutable {@link S2LatLngRect} and the mutable {@link
 * S2LatLngRect.Builder}.
 */
@GwtCompatible(serializable = false)
public abstract strictfp class S2LatLngRectBase implements S2Region, Serializable {
  protected final R1Interval lat;
  protected final S1Interval lng;

  /**
   * Constructs a rectangle from minimum and maximum latitudes and longitudes. If lo.lng() >
   * hi.lng(), the rectangle spans the 180 degree longitude line. Both points must be normalized,
   * with lo.lat() <= hi.lat(). The rectangle contains all the points p such that 'lo' <= p <= 'hi',
   * where '<=' is defined in the obvious way.
   */
  S2LatLngRectBase(final S2LatLng lo, final S2LatLng hi) {
    lat = new R1Interval(lo.lat().radians(), hi.lat().radians());
    lng = new S1Interval(lo.lng().radians(), hi.lng().radians());
    // assert (isValid());
  }

  /**
   * Constructs a rectangle from latitude and longitude intervals. The two intervals must either be
   * both empty or both non-empty, and the latitude interval must not extend outside [-90, +90]
   * degrees. Note that both intervals (and hence the rectangle) are closed.
   */
  S2LatLngRectBase(R1Interval lat, S1Interval lng) {
    this.lat = lat;
    this.lng = lng;
    // assert (isValid());
  }

  /**
   * Constructs a rectangle with lat and lng fields set to empty intervals, as defined in {@link
   * R1Interval} and {@link S1Interval}.
   */
  S2LatLngRectBase() {
    lat = R1Interval.empty();
    lng = S1Interval.empty();
  }

  /**
   * Returns true if the rectangle is valid, which essentially just means that the latitude bounds
   * do not exceed Pi/2 in absolute value and the longitude bounds do not exceed Pi in absolute
   * value. Also, if either the latitude or longitude bound is empty then both must be.
   */
  public final boolean isValid() {
    // The lat/lng ranges must either be both empty or both non-empty.
    return (Math.abs(lat.lo()) <= S2.M_PI_2
        && Math.abs(lat.hi()) <= S2.M_PI_2
        && lng.isValid()
        && lat.isEmpty() == lng.isEmpty());
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

  /** Returns true if the rectangle is a point, i.e. lo() == hi() */
  public final boolean isPoint() {
    return lat().lo() == lat().hi() && lng.lo() == lng().hi();
  }

  /**
   * Returns true if lng_.lo() > lng_.hi(), i.e. the rectangle crosses the 180 degree latitude line.
   */
  public final boolean isInverted() {
    return lng.isInverted();
  }

  /**
   * Returns the k<super>th</super> vertex of the rectangle (k = 0,1,2,3) in CCW order (lower-left,
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
    // The algorithm here is the same as in getDistance(S2LatLngRect), only with simplified
    // calculations.
    S2LatLngRectBase a = this;

    Preconditions.checkState(!a.isEmpty());
    Preconditions.checkArgument(p.isValid());

    if (a.lng().contains(p.lng().radians())) {
      return S1Angle.radians(
          Math.max(
              0.0, Math.max(p.lat().radians() - a.lat().hi(), a.lat().lo() - p.lat().radians())));
    }

    S1Interval interval = new S1Interval(a.lng().hi(), a.lng().complement().getCenter());
    double aLng = a.lng().lo();
    if (interval.contains(p.lng().radians())) {
      aLng = a.lng().hi();
    }

    S2Point lo = S2LatLng.fromRadians(a.lat().lo(), aLng).toPoint();
    S2Point hi = S2LatLng.fromRadians(a.lat().hi(), aLng).toPoint();
    S2Point loCrossHi = S2LatLng.fromRadians(0, aLng - S2.M_PI_2).normalized().toPoint();
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
      S1Angle lo;
      S1Angle hi;
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
    S1Angle aLng;
    S1Angle bLng;
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
    S2Point aLoCrossHi = S2LatLng.fromRadians(0, aLng.radians() - S2.M_PI_2).normalized().toPoint();
    S2Point bLo = new S2LatLng(b.latLo(), bLng).toPoint();
    S2Point bHi = new S2LatLng(b.latHi(), bLng).toPoint();
    S2Point bLoCrossHi = S2LatLng.fromRadians(0, bLng.radians() - S2.M_PI_2).normalized().toPoint();

    return S1Angle.min(
        S2EdgeUtil.getDistance(aLo, bLo, bHi, bLoCrossHi),
        S1Angle.min(
            S2EdgeUtil.getDistance(aHi, bLo, bHi, bLoCrossHi),
            S1Angle.min(
                S2EdgeUtil.getDistance(bLo, aLo, aHi, aLoCrossHi),
                S2EdgeUtil.getDistance(bHi, aLo, aHi, aLoCrossHi))));
  }

  /**
   * Returns the undirected Hausdorff distance (measured along the surface of the sphere) to the
   * given S2LatLngRect. The directed Hausdorff distance from rectangle A to rectangle B is given by
   * {@code h(A, B) = max_{p in A} min_{q in B} d(p, q)}. The Hausdorff distance between rectangle A
   * and rectangle B is given by {@code H(A, B) = max{h(A, B), h(B, A)}}.
   */
  public final S1Angle getHausdorffDistance(S2LatLngRectBase other) {
    return S1Angle.max(
        getDirectedHausdorffDistance(other), other.getDirectedHausdorffDistance(this));
  }

  /**
   * Returns the directed Hausdorff distance (measured along the surface of the sphere) to the given
   * S2LatLngRect. The directed Hausdorff distance from rectangle A to rectangle B is given by
   * {@code h(A, B) = max_{p in A} min_{q in B} d(p, q)}. The Hausdorff distance between rectangle A
   * and rectangle B is given by {@code H(A, B) = max{h(A, B), h(B, A)}}.
   */
  public final S1Angle getDirectedHausdorffDistance(S2LatLngRectBase other) {
    if (isEmpty()) {
      return S1Angle.radians(0);
    }
    if (other.isEmpty()) {
      return S1Angle.radians(S2.M_PI); // maximum possible distance on S2
    }

    double lngDistance = lng().getDirectedHausdorffDistance(other.lng());
    // assert lngDistance >= 0;
    return getDirectedHausdorffDistance(lngDistance, lat(), other.lat());
  }

  /**
   * Return the directed Hausdorff distance from one longitudinal edge spanning latitude range
   * {@code a_lat} to the other longitudinal edge spanning latitude range {@code b_lat}, with their
   * longitudinal difference given by {@code lngDiff}.
   */
  private static S1Angle getDirectedHausdorffDistance(double lngDiff, R1Interval a, R1Interval b) {
    // By symmetry, we can assume a's longitude is 0 and b's longitude is lngDiff. Call b's two
    // endpoints bLo and bHi. Let H be the hemisphere containing a and delimited by the longitude
    // line of b. The Voronoi diagram of b on H has three edges (portions of great circles) all
    // orthogonal to b and meeting at bLo cross bHi.
    //
    // E1: (bLo, bLo cross bHi)
    // E2: (bHi, bLo cross bHi)
    // E3: (-b_mid, bLo cross bHi), where b_mid is the midpoint of b
    //
    // They subdivide H into three Voronoi regions. Depending on how longitude 0 (which contains
    // edge a) intersects these regions, we distinguish two cases:
    // Case 1: it intersects three regions. This occurs when lngDiff <= M_PI_2.
    // Case 2: it intersects only two regions. This occurs when lngDiff > M_PI_2.
    //
    // In the first case, the directed Hausdorff distance to edge b can only be realized by the
    // following points on a:
    // A1: two endpoints of a.
    // A2: intersection of a with the equator, if b also intersects the equator.
    //
    // In the second case, the directed Hausdorff distance to edge b can only be realized by the
    // following points on a:
    // B1: two endpoints of a.
    // B2: intersection of a with E3
    // B3: farthest point from bLo to the interior of D, and farthest point from bHi to the interior
    // of U, if any, where D (resp. U) is the portion of edge a below (resp. above) the intersection
    // point from B2.

    // assert lngDiff >= 0;
    // assert lngDiff <= S2.M_PI;

    if (lngDiff == 0) {
      return S1Angle.radians(a.getDirectedHausdorffDistance(b));
    }

    // Assumed longitude of b.
    double bLng = lngDiff;
    // Two endpoints of b.
    S2Point bLo = S2LatLng.fromRadians(b.lo(), bLng).toPoint();
    S2Point bHi = S2LatLng.fromRadians(b.hi(), bLng).toPoint();

    // Handling of each case outlined at the top of the function starts here.
    // This is initialized a few lines below.
    S1Angle maxDistance;

    // Cases A1 and B1.
    S2Point aLo = S2LatLng.fromRadians(a.lo(), 0).toPoint();
    S2Point aHi = S2LatLng.fromRadians(a.hi(), 0).toPoint();
    maxDistance = S2EdgeUtil.getDistance(aLo, bLo, bHi);
    maxDistance = S1Angle.max(maxDistance, S2EdgeUtil.getDistance(aHi, bLo, bHi));

    if (lngDiff <= S2.M_PI_2) {
      // Case A2.
      if (a.contains(0) && b.contains(0)) {
        maxDistance = S1Angle.max(maxDistance, S1Angle.radians(lngDiff));
      }
    } else {
      // Case B2.
      S2Point p = getBisectorIntersection(b, bLng);
      double pLat = S2LatLng.latitude(p).radians();
      if (a.contains(pLat)) {
        maxDistance = S1Angle.max(maxDistance, new S1Angle(p, bLo));
      }

      // Case B3.
      if (pLat > a.lo()) {
        maxDistance =
            S1Angle.max(
                maxDistance,
                getInteriorMaxDistance(new R1Interval(a.lo(), Math.min(pLat, a.hi())), bLo));
      }
      if (pLat < a.hi()) {
        maxDistance =
            S1Angle.max(
                maxDistance,
                getInteriorMaxDistance(new R1Interval(Math.max(pLat, a.lo()), a.hi()), bHi));
      }
    }

    return maxDistance;
  }

  // A vector orthogonal to longitude 0.
  private static final S2Point ORTHO_LNG = S2Point.Y_NEG;

  /**
   * Return the intersection of longitude 0 with the bisector of an edge on longitude 'lng' and
   * spanning latitude range 'lat'.
   */
  private static S2Point getBisectorIntersection(R1Interval lat, double lng) {
    lng = Math.abs(lng);
    double latCenter = lat.getCenter();

    // A vector orthogonal to the bisector of the given longitudinal edge.
    S2LatLng orthoBisector;
    if (latCenter >= 0) {
      orthoBisector = S2LatLng.fromRadians(latCenter - S2.M_PI_2, lng);
    } else {
      orthoBisector = S2LatLng.fromRadians(-latCenter - S2.M_PI_2, lng - S2.M_PI);
    }
    return S2.robustCrossProd(ORTHO_LNG, orthoBisector.toPoint());
  }

  /**
   * Return max distance from a point b to the segment spanning latitude range aLat on longitude 0,
   * if the max occurs in the interior of aLat. Otherwise return -1.
   */
  private static S1Angle getInteriorMaxDistance(R1Interval aLat, S2Point b) {
    // Longitude 0 is in the y=0 plane. b.x() >= 0 implies that the maximum
    // does not occur in the interior of aLat.
    if (aLat.isEmpty() || b.getX() >= 0) {
      return S1Angle.radians(-1);
    }

    // Project b to the y=0 plane. The antipodal of the normalized projection is
    // the point at which the maxium distance from b occurs, if it is contained
    // in aLat.
    S2Point intersectionPoint = new S2Point(-b.getX(), 0, -b.getZ()).normalize();
    if (aLat.interiorContains(S2LatLng.latitude(intersectionPoint).radians())) {
      return new S1Angle(b, intersectionPoint);
    } else {
      return S1Angle.radians(-1);
    }
  }

  /**
   * Returns the width and height of this rectangle in latitude-longitude space. Empty rectangles
   * have a negative width and height.
   */
  public final S2LatLng getSize() {
    return S2LatLng.fromRadians(lat.getLength(), lng.getLength());
  }

  // Returns the true centroid of the rectangle multiplied by its surface area (see s2centroids.h
  // for details on centroids).  The result is not unit length, so you may want to normalize it.
  // Note that in general the centroid is *not* at the
  // center of the rectangle, and in fact it may not even be contained by the rectangle.  (It is the
  // "center of mass" of the rectangle viewed as
  // subset of the unit sphere, i.e. it is the point in space about which this curved shape would
  // rotate.)
  //
  // The reason for multiplying the result by the rectangle area is to make it easier to compute the
  // centroid of more complicated shapes.  The centroid of a union of disjoint regions can be
  // computed simply by adding their GetCentroid() results.
  public final S2Point getCentroid() {
    // When a sphere is divided into slices of constant thickness by a set of parallel planes, all
    // slices have the same surface area.  This implies that the z-component of the centroid is
    // simply the midpoint of the z-interval spanned by the S2LatLngRect.
    //
    // Similarly, it is easy to see that the (x,y) of the centroid lies in the plane through the
    // midpoint of the rectangle's longitude interval.  We only need to determine the distance "d"
    // of this point from the z-axis.
    //
    // Let's restrict our attention to a particular z-value.  In this z-plane, the S2LatLngRect is a
    // circular arc.  The centroid of this arc lies on a radial line through the midpoint of the
    // arc, and at a distance from the z-axis of
    //
    //     r * (sin(alpha) / alpha)
    //
    // where r = sqrt(1-z^2) is the radius of the arc, and "alpha" is half of the arc length (i.e.,
    // the arc covers longitudes [-alpha, alpha]).
    //
    // To find the centroid distance from the z-axis for the entire rectangle, we just need to
    // integrate over the z-interval.  This gives
    //
    //    d = Integrate[sqrt(1-z^2)*sin(alpha)/alpha, z1..z2] / (z2 - z1)
    //
    // where [z1, z2] is the range of z-values covered by the rectangle.  This simplifies to
    //
    //    d = sin(alpha)/(2*alpha*(z2-z1))*(z2*r2 - z1*r1 + theta2 - theta1)
    //
    // where [theta1, theta2] is the latitude interval, z1=sin(theta1), z2=sin(theta2),
    // r1=cos(theta1), and r2=cos(theta2).
    //
    // Finally, we want to return not the centroid itself, but the centroid scaled by the area of
    // the rectangle.  The area of the rectangle is
    //
    //    A = 2 * alpha * (z2 - z1)
    //
    // which fortunately appears in the denominator of "d".

    if (isEmpty()) {
      return new S2Point();
    }
    double z1 = Math.sin(latLo().radians());
    double z2 = Math.sin(latHi().radians());
    double r1 = Math.cos(latLo().radians());
    double r2 = Math.cos(latHi().radians());
    double alpha = 0.5 * lng.getLength();
    double r = Math.sin(alpha) * (r2 * z2 - r1 * z1 + lat.getLength());
    double lngCenter = lng.getCenter();
    double z = alpha * (z2 + z1) * (z2 - z1); // scaled by the area
    return new S2Point(r * Math.cos(lngCenter), r * Math.sin(lngCenter), z);
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
    return (lat.interiorContains(ll.lat().radians()) && lng.interiorContains(ll.lng().radians()));
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
    return (lat.interiorContains(other.lat) && lng.interiorContains(other.lng));
  }

  /** Returns true if this rectangle and the given other rectangle have any points in common. */
  public final boolean intersects(S2LatLngRectBase other) {
    return lat.intersects(other.lat) && lng.intersects(other.lng);
  }

  /**
   * Returns true if this rectangle intersects the given cell. (This is an exact test and may be
   * fairly expensive, see also MayIntersect below.)
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
      S1Interval edgeLng =
          S1Interval.fromPointPair(cellLl[i].lng().radians(), cellLl[(i + 1) & 3].lng().radians());
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

  /** Returns true if the boundary of this rectangle intersects the given geodesic edge (v0, v1). */
  public final boolean boundaryIntersects(S2Point v0, S2Point v1) {
    if (isEmpty()) {
      return false;
    }
    if (!lng.isFull()) {
      if (intersectsLngEdge(v0, v1, lat, lng.lo())) {
        return true;
      }
      if (intersectsLngEdge(v0, v1, lat, lng.hi())) {
        return true;
      }
    }
    if (lat.lo() != -S2.M_PI_2 && intersectsLatEdge(v0, v1, lat.lo(), lng)) {
      return true;
    }
    if (lat.hi() != S2.M_PI_2 && intersectsLatEdge(v0, v1, lat.hi(), lng)) {
      return true;
    }
    return false;
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
   * the given tolerance. See {@link R1Interval} and {@link S1Interval} for details.
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

    double poleZ;
    double poleAngle;
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
   * an exact test. Note, however, that the cell must be valid; an error may result if e.g. {@link
   * S2CellId#sentinel()} is passed here.
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
  @Override
  public final boolean contains(S2Point p) {
    return contains(new S2LatLng(p));
  }

  /** Returns true if the edge AB intersects the given edge of constant longitude. */
  public static final boolean intersectsLngEdge(S2Point a, S2Point b, R1Interval lat, double lng) {
    // Return true if the segment AB intersects the given edge of constant longitude. The nice thing
    // about edges of constant longitude is that they are straight lines on the sphere (geodesics).

    return S2.simpleCrossing(
        a,
        b,
        S2LatLng.fromRadians(lat.lo(), lng).toPoint(),
        S2LatLng.fromRadians(lat.hi(), lng).toPoint());
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
    S1Interval abTheta =
        S1Interval.fromPointPair(
            Math.atan2(a.dotProd(y), a.dotProd(x)), Math.atan2(b.dotProd(y), b.dotProd(x)));

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
