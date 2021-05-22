/*
 * Copyright 2019 Google Inc.
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
import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.geometry.S2ShapeUtil.CentroidMeasure;
import java.util.AbstractList;
import java.util.Comparator;
import java.util.List;

/**
 * Defines various angle and area measures for {@link S2Shape} objects. Unlike the built-in {@link
 * S2Polygon} and {@link S2Polyline} methods, these methods allow the underlying data to be
 * represented arbitrarily.
 */
@GwtCompatible
final strictfp class S2ShapeMeasures {

  private S2ShapeMeasures() {}

  /**
   * Returns the sum of all polyline lengths on the unit sphere for shapes of dimension 1, or {@link
   * S1Angle#ZERO} otherwise. See {@link #perimeter(S2Shape)} for shapes of dimension 2.
   *
   * <p>See {@link S2ShapeIndexMeasures#length(S2ShapeIndex)} for more info.
   */
  public static S1Angle length(S2Shape shape) {
    if (shape.dimension() != 1) {
      return S1Angle.ZERO;
    }
    S1Angle.Builder builder = new S1Angle.Builder();
    for (int chainId = 0; chainId < shape.numChains(); chainId++) {
      builder.add(polylineLength(shape, chainId));
    }
    return builder.build();
  }

  /**
   * Returns the length of the polyline, or {@link S1Angle#ZERO} if the polyline has fewer than two
   * vertices.
   */
  @VisibleForTesting
  static S1Angle polylineLength(S2Shape shape, int chainId) {
    S1Angle.Builder builder = new S1Angle.Builder();
    forEachChainEdge(shape, chainId, (a, b) -> builder.add(a.angle(b)));
    return builder.build();
  }

  /**
   * Returns the sum of all loop perimeters on the unit sphere for shapes of dimension 2, or {@link
   * S1Angle#ZERO} otherwise. See {@link #length(S2Shape)} for shapes of dimension 1.
   */
  public static S1Angle perimeter(S2Shape shape) {
    if (shape.dimension() != 2) {
      return S1Angle.ZERO;
    }
    S1Angle.Builder builder = new S1Angle.Builder();
    for (int chainId = 0; chainId < shape.numChains(); chainId++) {
      builder.add(loopPerimeter(shape, chainId));
    }
    return builder.build();
  }

  /** Returns the perimeter of the loop, or {@link S1Angle#ZERO} if the loop has 0 or 1 vertices. */
  @VisibleForTesting
  static S1Angle loopPerimeter(S2Shape shape, int chainId) {
    if (shape.getChainLength(chainId) <= 1) {
      return S1Angle.ZERO;
    }
    S1Angle.Builder builder = new S1Angle.Builder();
    forEachChainEdge(shape, chainId, (a, b) -> builder.add(a.angle(b)));
    return builder.build();
  }

  /**
   * For shapes of dimension 2, returns the area of the shape on the unit sphere. The result is
   * between 0 and 4*Pi steradians. Otherwise returns zero. This method has good relative accuracy
   * for both very large and very small regions.
   */
  public static double area(S2Shape shape) {
    if (shape.dimension() != 2) {
      return 0;
    }
    double area = 0;
    for (int chainId = 0; chainId < shape.numChains(); chainId++) {
      area += signedLoopArea(shape, chainId);
    }
    // Note that signedLoopArea() guarantees that the full loop (containing all points on the
    // sphere) has a very small negative area.
    if (area < 0.0) {
      area += 4 * S2.M_PI;
    }
    return area;
  }

  /**
   * Returns the area of the loop interior, i.e. the region on the left side of the loop. The result
   * is between 0 and 4*Pi steradians. The implementation ensures that nearly-degenerate clockwise
   * loops have areas close to zero, while nearly-degenerate counter-clockwise loops have areas
   * close to 4*Pi.
   */
  @VisibleForTesting
  static double loopArea(S2Shape shape, int chainId) {
    return loopArea(vertices(shape, chainId));
  }

  /** Same as {@link #loopArea(S2Shape, int)}, but takes a loop as a list of vertices. */
  static double loopArea(List<S2Point> loop) {
    double area = signedLoopArea(loop);
    assert (Math.abs(area) <= 2 * S2.M_PI);
    if (area < 0.0) {
      area += 4 * S2.M_PI;
    }
    return area;
  }

  /**
   * Returns the area of the loop interior, i.e. the region on the left side of the loop. The result
   * is between 0 and 4*Pi steradians. The implementation ensures that nearly-degenerate clockwise
   * loops have areas close to zero, while nearly-degenerate counter-clockwise loops have areas
   * close to 4*Pi.
   */
  private static double signedLoopArea(S2Shape shape, int chainId) {
    return signedLoopArea(vertices(shape, chainId));
  }

  /** Same as {@link #signedLoopArea(S2Shape, int)}, but takes a loop as a list of vertices. */
  private static double signedLoopArea(List<S2Point> loop) {
    MutableDouble mutableArea = new MutableDouble();
    S2ShapeUtil.visitSurfaceIntegral(loop, (a, b, c) -> mutableArea.d += S2.signedArea(a, b, c));
    double maxError = S2.getTurningAngleMaxError(loop.size());
    assert (Math.abs(mutableArea.d) <= 4 * S2.M_PI + maxError);
    double area = mutableArea.d % (4 * S2.M_PI);

    if (area == -2 * S2.M_PI) {
      area = 2 * S2.M_PI;
    }

    if (Math.abs(area) <= maxError) {
      double curvature = turningAngle(loop);
      // Zero-area loops should have a curvature of approximately +/- 2*Pi.
      assert (!(area == 0 && curvature == 0));
      if (curvature == 2 * S2.M_PI) {
        // Degenerate
        return 0.0;
      }
      if (area <= 0 && curvature > 0) {
        return Double.MIN_VALUE;
      }
      // Full loops are handled by the case below.
      if (area >= 0 && curvature < 0) {
        return -Double.MIN_VALUE;
      }
    }
    return area;
  }

  static double turningAngle(S2Shape shape, int chainId) {
    return turningAngle(vertices(shape, chainId));
  }

  /**
   * Returns the geodesic curvature of the loop, defined as the sum of the turn angles at each
   * vertex (see {@link S2#turnAngle(S2Point, S2Point, S2Point)}). The result is positive if the
   * loop is counter-clockwise, negative if the loop is clockwise, and zero if the loop is a great
   * circle. The geodesic curvature is equal to 2*Pi minus the area of the loop.
   *
   * <p>The following cases are handled specially:
   *
   * <ul>
   *   <li>Degenerate loops (consisting of an isolated vertex or composed entirely of sibling edge
   *       pairs) have a curvature of 2*Pi exactly.
   *   <li>The full loop (containing all points, and represented as a loop with no vertices) has a
   *       curvature of -2*Pi exactly.
   *   <li>All other loops have a non-zero curvature in the range (-2*Pi, 2*Pi). For any such loop,
   *       reversing the order of the vertices is guaranteed to negate the curvature. This property
   *       can be used to define a unique normalized orientation for every loop.
   * </ul>
   */
  static double turningAngle(List<S2Point> loop) {
    // By convention, a loop with no vertices contains all points on the sphere.
    if (loop.isEmpty()) {
      return -2 * S2.M_PI;
    }

    // Remove any degeneracies from the loop.
    loop = pruneDegeneracies(loop);

    // If the entire loop was degenerate, it's turning angle is defined as 2*Pi.
    if (loop.isEmpty()) {
      return 2 * S2.M_PI;
    }

    // To ensure that we get the same result when the vertex order is rotated, and that the result
    // is negated when the vertex order is reversed, we need to add up the individual turn angles in
    // a consistent order. (In general, adding up a set of numbers in a different order can change
    // the sum due to rounding errors).
    //
    // Furthermore, if we just accumulate an ordinary sum then the worst-case error is quadratic in
    // the number of vertices. (This can happen with spiral shapes, where the partial sum of the
    // turning angles can be linear in the number of vertices). To avoid this we use the Kahan
    // summation algorithm (http://en.wikipedia.org/wiki/Kahan_summation_algorithm).
    LoopOrder loopOrder = canonicalLoopOrder(loop);
    int i = loopOrder.first;
    int dir = loopOrder.dir;
    int n = loop.size();
    double sum =
        S2.turnAngle(loop.get((i + n - dir) % n), loop.get(i % n), loop.get((i + dir) % n));
    double compensation = 0; // Kahan summation algorithm
    for (int x = 0; x < n - 1; x++) {
      i += dir;
      double angle =
          S2.turnAngle(loop.get((i - dir) % n), loop.get(i % n), loop.get((i + dir) % n));
      double oldSum = sum;
      angle += compensation;
      sum += angle;
      compensation = (oldSum - sum) + angle;
    }
    double maxCurvature = 2 * S2.M_PI - 4 * S2.DBL_EPSILON;
    sum += compensation;
    return Math.max(-maxCurvature, Math.min(maxCurvature, dir * sum));
  }

  /**
   * Returns an index "first" and a direction "dir" such that the vertex sequence (first, first +
   * dir, ..., first + (n - 1) * dir) does not change when the loop vertex order is rotated or
   * reversed. This allows the loop vertices to be traversed in a canonical order.
   */
  @VisibleForTesting
  static LoopOrder canonicalLoopOrder(List<S2Point> loop) {
    // In order to handle loops with duplicate vertices and/or degeneracies, we return the LoopOrder
    // that minimizes the entire corresponding vertex *sequence*. For example, suppose that vertices
    // are sorted alphabetically, and consider the loop CADBAB. The canonical loop order would be
    // (4, 1), corresponding to the vertex sequence ABCADB. (For comparison, loop order (4, -1)
    // yields the sequence ABDACB).
    //
    // If two or more loop orders yield identical minimal vertex sequences, then it doesn't matter
    // which one we return (since they yield the same result).

    // For efficiency, we divide the process into two steps. First we find the smallest vertex, and
    // the set of vertex indices where that vertex occurs (noting that the loop may contain
    // duplicate vertices). Then we consider both possible directions starting from each such vertex
    // index, and return the LoopOrder corresponding to the smallest vertex sequence.
    if (loop.isEmpty()) {
      return new LoopOrder(0, 1);
    }

    List<Integer> minIndexes = Lists.newArrayList(0);
    for (int i = 1; i < loop.size(); i++) {
      if (loop.get(i).compareTo(loop.get(minIndexes.get(0))) <= 0) {
        if (loop.get(i).compareTo(loop.get(minIndexes.get(0))) < 0) {
          minIndexes.clear();
        }
        minIndexes.add(i);
      }
    }
    LoopOrder minOrder = new LoopOrder(minIndexes.get(0), 1);
    Comparator<LoopOrder> loopOrderComparator = new LoopOrderComparator(loop);
    for (int minIndex : minIndexes) {
      LoopOrder loopOrder1 = new LoopOrder(minIndex, 1);
      LoopOrder loopOrder2 = new LoopOrder(minIndex + loop.size(), -1);
      if (loopOrderComparator.compare(loopOrder1, minOrder) < 0) {
        minOrder = loopOrder1;
      }
      if (loopOrderComparator.compare(loopOrder2, minOrder) < 0) {
        minOrder = loopOrder2;
      }
    }
    return minOrder;
  }

  /**
   * Returns a new loop obtained by removing all degeneracies from "input". In particular, the
   * result will not contain any adjacent duplicate vertices or sibling edge pairs, i.e. vertex
   * sequences of the form (A, A) or (A, B, A).
   */
  @VisibleForTesting
  static List<S2Point> pruneDegeneracies(List<S2Point> input) {
    List<S2Point> loop = Lists.newArrayListWithCapacity(input.size());
    for (S2Point p : input) {
      if (loop.isEmpty() || !p.equalsPoint(Iterables.getLast(loop))) {
        if (loop.size() >= 2 && p.equalsPoint(loop.get(loop.size() - 2))) {
          loop.remove(loop.size() - 1);
        } else {
          loop.add(p);
        }
      }
    }

    // Check whether the loop was completely degenerate.
    if (loop.size() < 3) {
      return ImmutableList.of();
    }

    // Otherwise some portion of the loop is guaranteed to be non-degenerate.
    // However there may still be some degenerate portions to remove.
    if (loop.get(0).equalsPoint(Iterables.getLast(loop))) {
      loop.remove(loop.size() - 1);
    }

    // If the loop begins with BA and ends with A, then there is an edge pair of the form ABA at the
    // end/start of the loop. Remove all such pairs. As noted above, this is guaranteed to leave a
    // non-degenerate loop.
    int i = 0;
    while (loop.get(i + 1).equalsPoint(loop.get(loop.size() - i - 1))) {
      i++;
    }
    return loop.subList(i, loop.size() - i);
  }

  /**
   * Returns the centroid of shape multiplied by the measure of shape.
   *
   * <p>See {@link S2ShapeIndexMeasures#centroid(S2ShapeIndex)} for more info.
   */
  public static S2Point centroid(S2Shape shape) {
    S2Point.Builder builder = new S2Point.Builder();
    int dimension = shape.dimension();
    int numChains = shape.numChains();
    switch (dimension) {
      case 0:
        for (int chainId = 0; chainId < numChains; chainId++) {
          builder.add(shape.getChainVertex(chainId, 0));
        }
        break;
      case 1:
        for (int chainId = 0; chainId < numChains; chainId++) {
          builder.add(polylineCentroid(shape, chainId));
        }
        break;
      case 2:
        for (int chainId = 0; chainId < numChains; chainId++) {
          builder.add(loopCentroid(shape, chainId));
        }
        break;
      default:
        throw new IllegalArgumentException("Unexpected S2Shape dimension: " + shape.dimension());
    }
    return builder.build();
  }

  /**
   * Returns the true centroid of the polyline multiplied by the length of the polyline.
   *
   * <p>Scaling by the polyline length makes it easy to compute the centroid of several polylines
   * (by simply adding up their centroids).
   *
   * <p>CAVEAT: Returns {@link S2Point#ORIGIN} for degenerate polylines (e.g., AA). [Note that this
   * answer is correct; the result of this function is a line integral over the polyline, whose
   * value is always zero if the polyline is degenerate].
   */
  @VisibleForTesting
  static S2Point polylineCentroid(S2Shape shape, int chainId) {
    S2Point.Builder builder = new S2Point.Builder();
    forEachChainEdge(shape, chainId, (a, b) -> builder.add(S2.trueCentroid(a, b)));
    return builder.build();
  }

  /**
   * Returns the true centroid of the loop multiplied by the area of the loop.
   *
   * <p>See {@link S2ShapeIndexMeasures#centroid(S2ShapeIndex)} for more info.
   */
  @VisibleForTesting
  static S2Point loopCentroid(S2Shape shape, int chainId) {
    CentroidMeasure centroidMeasure = new CentroidMeasure();
    S2ShapeUtil.visitSurfaceIntegral(
        new AbstractList<S2Point>() {
          @Override
          public S2Point get(int i) {
            return shape.getChainVertex(chainId, i);
          }

          @Override
          public int size() {
            return shape.getChainLength(chainId);
          }
        },
        centroidMeasure);
    return centroidMeasure.value();
  }

  private static List<S2Point> vertices(S2Shape shape, int chainId) {
    return new AbstractList<S2Point>() {
      int size = shape.getChainLength(chainId);

      @Override
      public S2Point get(int i) {
        return shape.getChainVertex(chainId, i);
      }

      @Override
      public int size() {
        return size;
      }
    };
  }

  /** Passes each edge (a, b) in the chain of shape at index chainId to edgeConsumer. */
  private static void forEachChainEdge(
      S2Shape shape, int chainId, BiConsumer<S2Point, S2Point> edgeConsumer) {
    int chainLength = shape.getChainLength(chainId);
    if (chainLength == 0) {
      return;
    }
    S2Point prev = shape.getChainVertex(chainId, 0);
    for (int edgeOffset = 1; edgeOffset <= chainLength; edgeOffset++) {
      S2Point next = shape.getChainVertex(chainId, edgeOffset);
      edgeConsumer.accept(prev, next);
      prev = next;
    }
  }

  /**
   * Represents a cyclic ordering of the loop vertices, starting at the index "first" and proceeding
   * in direction "dir" (either +1 or -1). "first" and "dir" must be chosen such that (first, ...,
   * first + n * dir) are all in the range [0, 2*n-1].
   */
  @VisibleForTesting
  static class LoopOrder {
    final int first;
    final int dir;

    LoopOrder(int first, int dir) {
      this.first = first;
      this.dir = dir;
    }

    boolean equalsLoopOrder(LoopOrder other) {
      return first == other.first && dir == other.dir;
    }

    @Override
    public boolean equals(Object other) {
      return other instanceof LoopOrder && equalsLoopOrder((LoopOrder) other);
    }

    @Override
    public int hashCode() {
      return first + dir;
    }
  }

  private static class LoopOrderComparator implements Comparator<LoopOrder> {
    private final List<S2Point> loop;
    private final IntFunction<S2Point> vertex;

    LoopOrderComparator(List<S2Point> loop) {
      this.loop = loop;
      vertex = i -> loop.get(i % loop.size());
    }

    @Override
    public int compare(LoopOrder loopOrder1, LoopOrder loopOrder2) {
      if (loopOrder1.equalsLoopOrder(loopOrder2)) {
        return 0;
      }
      assert (vertex.apply(loopOrder1.first).equalsPoint(vertex.apply(loopOrder2.first)));
      int i1 = loopOrder1.first;
      int i2 = loopOrder2.first;
      for (int n = loop.size(); --n > 0; ) {
        i1 += loopOrder1.dir;
        i2 += loopOrder2.dir;
        int compare = vertex.apply(i1).compareTo(vertex.apply(i2));
        if (compare != 0) {
          return compare;
        }
      }
      return 0;
    }
  }

  /** A consumer which accepts two arguments. */
  private interface BiConsumer<T, U> {
    void accept(T t, U u);
  }

  /** A function which accepts an int. */
  private interface IntFunction<T> {
    T apply(int i);
  }

  /** Wraps a mutable primitive double. */
  private static class MutableDouble {
    double d = 0;
  }
}
