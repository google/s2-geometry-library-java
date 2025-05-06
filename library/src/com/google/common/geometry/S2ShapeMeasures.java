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

import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import java.util.AbstractList;
import java.util.Comparator;
import java.util.List;

/**
 * Defines various angle and area measures for {@link S2Shape} objects. Unlike the built-in {@link
 * S2Polygon} and {@link S2Polyline} methods, these methods allow the underlying data to be
 * represented arbitrarily.
 */
public final class S2ShapeMeasures {

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
      area += 4 * PI;
    }
    return area;
  }

  /**
   * Like area(), except that this method is faster and has more error. The additional error is at
   * most 2.22e-15 steradians per vertex, which works out to about 0.09 square meters per vertex on
   * the Earth's surface. For example, a loop with 100 vertices has a maximum error of about 9
   * square meters. (The actual error is typically much smaller than this.)
   */
  public static double approxArea(S2Shape shape) {
    if (shape.dimension() != 2) {
      return 0;
    }
    double area = 0;
    for (int chainId = 0; chainId < shape.numChains(); ++chainId) {
      area += approxLoopArea(shape, chainId);
    }

    // Special case to ensure that full polygons are handled correctly.
    if (area <= 4 * PI) {
      return area;
    }
    return area % (4 * PI);
  }

  /**
   * Like loopArea(), except that this method is faster and has more error. The result is between 0
   * and 4*Pi steradians. The maximum error is 2.22e-15 steradians per loop vertex, which works out
   * to about 0.09 square meters per vertex on the Earth's surface. For example, a loop with 100
   * vertices has a maximum error of about 9 square meters. (The actual error is typically much
   * smaller than this.) The error bound can be computed using getTurningAngleMaxError(), which
   * returns the maximum error in steradians.
   */
  public static double approxLoopArea(S2Shape shape, int chainId) {
    return 2 * PI - turningAngle(shape, chainId);
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

  /** Same as {@link #loopArea(S2Shape, int)}, but takes a list of vertices. */
  static double loopArea(List<S2Point> loop) {
    double area = signedLoopArea(loop);
    assert abs(area) <= 2 * PI;
    if (area < 0.0) {
      area += 4 * PI;
    }
    return area;
  }

  /**
   * Returns either the positive area of the region on the left side of the loop, or the negative
   * area of the region on the right side of the loop, whichever is smaller in magnitude. The result
   * is between -2*Pi and 2*Pi steradians. This method is used to accurately compute the area of
   * polygons consisting of multiple loops.
   *
   * <p>The following cases are handled specially:
   *
   * <ul>
   *   <li>Counter-clockwise loops are guaranteed to have positive area, and clockwise loops are
   *       guaranteed to have negative area.
   *   <li>Degenerate loops (consisting of an isolated vertex or composed entirely of sibling edge
   *       pairs) have an area of exactly zero.
   *   <li>The full loop (containing all points, and represented as a loop with no vertices) has a
   *       negative area with the minimum possible magnitude. (This is the "signed equivalent" of
   *       having an area of 4*Pi.)
   * </ul>
   */
  private static double signedLoopArea(S2Shape shape, int chainId) {
    return signedLoopArea(vertices(shape, chainId));
  }

  /** Same as {@link #signedLoopArea(S2Shape, int)}, but takes a list of vertices. */
  private static double signedLoopArea(List<S2Point> loop) {
    // It is surprisingly difficult to compute the area of a loop robustly. The main issues are (1)
    // whether degenerate loops are considered to be CCW or not (i.e., whether their area is close
    // to 0 or 4*Pi), and (2) computing the areas of small loops with good relative accuracy.
    //
    // With respect to degeneracies, we would like getArea() to be consistent with
    // S2Loop.contains(S2Point) in that loops that contain many points should have large areas, and
    // loops that contain few points should have small areas. For example, if a degenerate triangle
    // is considered CCW according to S2Predicates.sign(), then it will contain very few points and
    // its area should be approximately zero. On the other hand if it is considered clockwise, then
    // it will contain virtually all points and so its area should be approximately 4*Pi.
    //
    // More precisely, let U be the set of S2Points for which S2.isUnitLength() is true, let P(U) be
    // the projection of those points onto the mathematical unit sphere, and let V(P(U)) be the
    // Voronoi diagram of the projected points. Then for every loop x, we would like GetArea() to
    // approximately equal the sum of the areas of the Voronoi regions of the points p for which
    // x.Contains(p) is true.
    //
    // The second issue is that we want to compute the area of small loops accurately. This
    // requires having good relative precision rather than good absolute precision. For example, if
    // the area of a loop is 1e-12 and the error is 1e-15, then the area only has 3 digits of
    // accuracy. (For reference, 1e-12 is about 40 square meters on the surface of the earth.) We
    // would like to have good relative accuracy even for small loops.
    //
    // To achieve these goals, we combine two different methods of computing the area. This first
    // method is based on the Gauss-Bonnet theorem, which says that the area enclosed by the loop
    // equals 2*Pi minus the total geodesic curvature of the loop (i.e., the sum of the "turning
    // angles" at all the loop vertices). The big advantage of this method is that as long as we
    // use S2Predicates.sign() to compute the turning angle at each vertex, then degeneracies are
    // always handled correctly. In other words, if a degenerate loop is CCW according to the
    // symbolic perturbations used by S2Predicates.sign(), then its turning angle will be
    // approximately 2*Pi.
    //
    // The disadvantage of the Gauss-Bonnet method is that its absolute error is about 2e-15 times
    // the number of vertices (see {@link getTurningAngleMaxError(S2Shape, int)}. So, it cannot
    // compute the area of small loops accurately.
    //
    // The second method is based on splitting the loop into triangles and summing the area of each
    // triangle. To avoid the difficulty and expense of decomposing the loop into a union of non-
    // overlapping triangles, instead we compute a signed sum over triangles that may overlap (see
    // the comments for S2ShapeUtil.visitSurfaceIntegral). The advantage of this method is that the
    // area of each triangle can be computed with much better relative accuracy (using l'Huilier's
    // theorem). The disadvantage is that the result is a signed area: CCW loops may yield a small
    // positive value, while CW loops may yield a small negative value (which is converted to a
    // positive area by adding 4*Pi). This means that small errors in computing the signed area may
    // translate into a very large error in the result (if the sign of the sum is incorrect).
    //
    // So, our strategy is to combine these two methods as follows. First we compute the area using
    // the "signed sum over triangles" approach (since it is generally more accurate). We also
    // estimate the maximum error in this result. If the signed area is too close to zero (i.e.,
    // zero is within the error bounds), then we double-check the sign of the result using the
    // Gauss-Bonnet method. If the two methods disagree, we return the smallest possible positive
    // or negative area based on the result of {@link turningAngle()}. Otherwise we return the area
    // that we computed originally.
    //
    // The signed area should be between approximately -4*Pi and 4*Pi.
    // Normalize it to be in the range [-2*Pi, 2*Pi].
    MutableDouble mutableArea = new MutableDouble();
    S2ShapeUtil.visitSurfaceIntegral(loop, (a, b, c) -> mutableArea.d += S2.signedArea(a, b, c));
    // TODO(user): This error estimate is approximate.
    double maxError = S2.getTurningAngleMaxError(loop.size());

    // Normalize the area to be in the range (-2*Pi, 2*Pi]. Effectively this means that hemispheres
    // are always interpreted as having positive area.
    assert abs(mutableArea.d) <= 4 * PI + maxError;
    double area = mutableArea.d % (4 * PI);
    if (area == -2 * PI) {
      area = 2 * PI;
    }

    // If the area is a small negative or positive number, verify that the sign of the result is
    // consistent with the loop orientation.
    if (abs(area) <= maxError) {
      double curvature = turningAngle(loop);
      // Zero-area loops should have a curvature of approximately +/- 2*Pi.
      assert !(area == 0 && curvature == 0);
      if (curvature == 2 * PI) {
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
      return -2 * PI;
    }

    // Remove any degeneracies from the loop.
    loop = pruneDegeneracies(loop);

    // If the entire loop was degenerate, it's turning angle is defined as 2*Pi.
    if (loop.isEmpty()) {
      return 2 * PI;
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
    double maxCurvature = 2 * PI - 4 * S2.DBL_EPSILON;
    sum += compensation;
    return max(-maxCurvature, min(maxCurvature, dir * sum));
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
      // Remove duplicate vertices.
      if (loop.isEmpty() || !p.equalsPoint(Iterables.getLast(loop))) {
        // Remove edge pairs of the form ABA.
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

    // Otherwise some portion of the loop is guaranteed to be non-degenerate. However there may
    // still be some degenerate portions to remove.
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
  public static S2Point polylineCentroid(S2Shape shape, int chainId) {
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
    return loopCentroid(vertices(shape, chainId));
  }

  /** Same as loopCentroid(S2Shape shape, int chainId) but takes a list of vertices. */
  @VisibleForTesting
  static S2Point loopCentroid(List<S2Point> loop) {
    double[] sum = new double[3];
    S2ShapeUtil.visitSurfaceIntegral(
        loop,
        (a, b, c) -> {
          S2Point centroid = S2.trueCentroid(a, b, c);
          sum[0] += centroid.x;
          sum[1] += centroid.y;
          sum[2] += centroid.z;
        });
    return new S2Point(sum[0], sum[1], sum[2]);
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
      assert vertex.apply(loopOrder1.first).equalsPoint(vertex.apply(loopOrder2.first));
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
