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

import static com.google.common.geometry.S2Predicates.sign;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * S2ConvexHullQuery builds the convex hull of any collection of points, polylines, loops, and
 * polygons. It returns a single convex loop.
 *
 * <p>The convex hull is defined as the smallest convex region on the sphere that contains all of
 * the input geometry. Recall that a region is "convex" if for every pair of points inside the
 * region, the straight edge between them is also inside the region. In our case, a "straight" edge
 * is a geodesic, i.e. the shortest path on the sphere between two points.
 *
 * <p>Containment of input geometry is defined as follows:
 *
 * <ul>
 *   <li>Each input loop and polygon is contained by the convex hull exactly (i.e., according to
 *       S2Polygon.contains(S2Polygon)).
 *   <li>Each input point is either contained by the convex hull or is a vertex of the convex hull.
 *       (Recall that S2Loops do not necessarily contain their vertices.)
 *   <li>For each input polyline, the convex hull contains all of its vertices according to the rule
 *       for points above. (The definition of convexity then ensures that the convex hull also
 *       contains the polyline edges.)
 * </ul>
 *
 * <p>To use this class, call the add*() methods to add your input geometry, and then call
 * getConvexHull(). Note that getConvexHull() does *not* reset the state; you can continue adding
 * geometry if desired and compute the convex hull again. If you want to start from scratch, simply
 * declare a new S2ConvexHullQuery object (they are cheap to create).
 *
 * <p>This class is not thread-safe.
 *
 * <p>This class implements Andrew's monotone chain algorithm, which is a variant of the Graham scan
 * (see https://en.wikipedia.org/wiki/Graham_scan). The time complexity is O(n log n), and the space
 * required is O(n). In fact only the call to "sort" takes O(n log n) time; the rest of the
 * algorithm is linear.
 *
 * <p>Demonstration of the algorithm and code:
 * en.wikibooks.org/wiki/Algorithm_Implementation/Geometry/Convex_hull/Monotone_chain
 */
public final class S2ConvexHullQuery {
  /** The length of edges to expand away from degenerate points to form a polygon. */
  private static final double OFFSET_FOR_SINGLE_POINT_LOOP = 1e-15;

  private final S2LatLngRect.Builder bound = S2LatLngRect.Builder.empty();
  private final List<S2Point> points = new ArrayList<>();

  /** Adds a point to the input geometry. */
  public void addPoint(S2Point point) {
    bound.addPoint(point);
    points.add(point);
  }

  /** Adds a polyline to the input geometry. */
  public void addPolyline(S2Polyline polyline) {
    bound.union(polyline.getRectBound());
    points.addAll(polyline.vertices());
  }

  /** Adds a loop to the input geometry. */
  public void addLoop(S2Loop loop) {
    // Only loops at depth 0 can contribute to the convex hull.
    if (loop.depth() != 0) {
      return;
    }
    bound.union(loop.getRectBound());
    if (loop.isEmptyOrFull()) {
      // The empty and full loops consist of a single fake "vertex" that should not be added to our
      // point collection.
      return;
    }
    for (int i = 0; i < loop.numVertices(); ++i) {
      points.add(loop.vertex(i));
    }
  }

  /** Adds a polygon to the input geometry. */
  public void addPolygon(S2Polygon polygon) {
    for (int i = 0; i < polygon.numLoops(); ++i) {
      addLoop(polygon.loop(i));
    }
  }

  /**
   * Computes a bounding cap for the input geometry provided.
   *
   * <p>Note that this method does not clear the geometry; you can continue adding to it and call
   * this method again if desired.
   */
  public S2Cap getCapBound() {
    // We keep track of a rectangular bound rather than a spherical cap because it is easy to
    // compute a tight bound for a union of rectangles, whereas it is quite difficult to compute a
    // tight bound around a union of caps. Also, polygons and polylines implement getCapBound() in
    // terms of getRectBound() for this same reason, so it is much better to keep track of a
    // rectangular bound as we go along and convert it at the end.
    //
    // TODO(user): We could compute an optimal bound by implementing Welzl's algorithm.
    // However we would still need to have special handling of loops and polygons, since if a loop
    // spans more than 180 degrees in any direction (i.e., if it contains two antipodal points),
    // then it is not enough just to bound its vertices. In this case the only convex bounding cap
    // is S2Cap.full(), and the only convex bounding loop is the full loop.
    return bound.getCapBound();
  }

  /**
   * Computes the convex hull of the input geometry provided.
   *
   * <p>If there is no geometry, this method returns an empty loop containing no points (see
   * S2Loop.isEmpty()).
   *
   * <p>If the geometry spans more than half of the sphere, this method returns a full loop
   * containing the entire sphere (see S2Loop.isFull()).
   *
   * <p>If the geometry contains 1 or 2 points, or a single edge, this method returns a very small
   * loop consisting of three vertices (which are a superset of the input vertices).
   *
   * <p>Note that this method does not clear the geometry; you can continue adding to it and call
   * this method again if desired.
   */
  public S2Loop getConvexHull() {
    S2Cap cap = getCapBound();
    if (cap.height() >= 1) {
      // The bounding cap is not convex. The current bounding cap implementation is not optimal,
      // but nevertheless it is likely that the input geometry itself is not contained by any convex
      // polygon. In any case, we need a convex bounding cap to proceed with the algorithm below
      // (in order to construct a point "origin" that is definitely outside the convex hull).
      return S2Loop.full();
    }
    // This code implements Andrew's monotone chain algorithm, which is a simple variant of the
    // Graham scan. Rather than sorting by x-coordinate, instead we sort the points in CCW order
    // around an origin O such that all points are guaranteed to be on one side of some geodesic
    // through O. This ensures that as we scan through the points, each new point can only belong at
    // the end of the chain (i.e., the chain is monotone in terms of the angle around O from the
    // starting point).
    S2Point origin = cap.axis().ortho();
    Collections.sort(points, new OrderedCcwAround(origin));

    // Remove duplicates. We need to do this before checking whether there are fewer than 3 points.
    ImmutableSet<S2Point> uniquePoints = ImmutableSet.copyOf(points);
    points.clear();
    points.addAll(uniquePoints);

    // Special cases for fewer than 3 points.
    if (points.size() < 3) {
      if (points.isEmpty()) {
        return S2Loop.empty();
      } else if (points.size() == 1) {
        return getSinglePointLoop(points.get(0));
      } else {
        return getSingleEdgeLoop(points.get(0), points.get(1));
      }
    }

    // Verify that all points lie within a 180 degree span around the origin.
    Preconditions.checkState(sign(origin, points.get(0), Iterables.getLast(points)) >= 0);

    // Generate the lower and upper halves of the convex hull. Each half consists of the maximal
    // subset of vertices such that the edge chain makes only left (CCW) turns.
    List<S2Point> lower = getMonotoneChain(points);
    List<S2Point> upper = getMonotoneChain(Lists.reverse(points));

    // Remove the duplicate vertices and combine the chains.
    Preconditions.checkState(lower.get(0).equals(Iterables.getLast(upper)));
    Preconditions.checkState(Iterables.getLast(lower).equals(upper.get(0)));
    lower.remove(lower.size() - 1);
    upper.remove(upper.size() - 1);
    lower.addAll(upper);
    return new S2Loop(lower);
  }

  /** A comparator for sorting points in CCW around a central point "center". */
  private static final class OrderedCcwAround implements Comparator<S2Point> {
    private final S2Point center;

    OrderedCcwAround(S2Point center) {
      this.center = center;
    }

    @Override
    public int compare(S2Point x, S2Point y) {
      if (lessThan(x, y)) {
        return -1;
      } else if (lessThan(y, x)) {
        return 1;
      } else {
        return 0;
      }
    }

    private boolean lessThan(S2Point x, S2Point y) {
      // If X and Y are equal, this will return false (as desired).
      return sign(center, x, y) > 0;
    }
  }

  /**
   * Iterate through the given points, selecting the maximal subset of points such that the edge
   * chain makes only left (CCW) turns.
   */
  private static List<S2Point> getMonotoneChain(List<S2Point> points) {
    List<S2Point> output = new ArrayList<>();
    for (S2Point p : points) {
      // Remove any points that would cause the chain to make a clockwise turn.
      while (output.size() >= 2
          && sign(output.get(output.size() - 2), Iterables.getLast(output), p) <= 0) {
        output.remove(output.size() - 1);
      }
      output.add(p);
    }
    return output;
  }

  /**
   * Constructs a 3-vertex polygon consisting of "p" and two nearby vertices. Note that contains(p)
   * may be false for the resulting loop (see comments at top of file).
   */
  private static S2Loop getSinglePointLoop(S2Point p) {
    S2Point d0 = S2.ortho(p);
    S2Point d1 = S2Point.crossProd(p, d0);
    return new S2Loop(
        ImmutableList.of(
            p,
            S2Point.normalize(S2Point.add(p, S2Point.mul(d0, OFFSET_FOR_SINGLE_POINT_LOOP))),
            S2Point.normalize(S2Point.add(p, S2Point.mul(d1, OFFSET_FOR_SINGLE_POINT_LOOP)))));
  }

  /** Construct a loop consisting of the two vertices and their midpoint. */
  private static S2Loop getSingleEdgeLoop(S2Point a, S2Point b) {
    // If the points are exactly antipodal we return the full loop.
    //
    // Note that we could use the code below even in this case (which would return a zero-area loop
    // that follows the edge AB), except that (1) the direction of AB is defined using symbolic
    // perturbations and therefore is not predictable by ordinary users, and (2) S2Loop disallows
    // antipodal adjacent vertices and so we would need to use 4 vertices to define the degenerate
    // loop. Note that the S2Loop antipodal vertex restriction is historical and now could easily
    // be removed, however it would still have the problem that the edge direction is not easily
    // predictable.
    if (a.equalsPoint(b.neg())) {
      return S2Loop.full();
    }

    // Construct a loop consisting of the two vertices and their midpoint. We use
    // S2EdgeUtil.interpolate() to ensure that the midpoint is very close to the edge even when its
    // endpoints are nearly antipodal.
    S2Loop loop = new S2Loop(ImmutableList.of(a, b, S2EdgeUtil.interpolate(0.5, a, b)));
    // The resulting loop may be clockwise, so invert it if necessary.
    loop.normalize();
    return loop;
  }
}
