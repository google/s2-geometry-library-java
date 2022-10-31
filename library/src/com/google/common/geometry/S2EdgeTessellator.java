/*
 * Copyright 2022 Google Inc.
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

import com.google.common.collect.Iterables;
import com.google.errorprone.annotations.CheckReturnValue;
import java.util.List;
import java.util.function.Consumer;

/**
 * Given an edge in some 2D projection (e.g., Mercator), S2EdgeTessellator converts the edge into a
 * chain of spherical geodesic edges such that the maximum distance between the original edge and
 * the geodesic edge chain is at most "tolerance". Similarly, it can convert a spherical geodesic
 * edge into a chain of edges in a given 2D projection such that the maximum distance between the
 * geodesic edge and the chain of projected edges is at most "tolerance".
 *
 * <p>Tessellation is implemented by subdividing the edge until the estimated maximum error is
 * below the given tolerance. Estimating error is a hard problem, especially when the only methods
 * available are point evaluation of the projection and its inverse. (These are the only methods
 * that {@link Projection} provides, which makes it easier and less error-prone to implement new
 * projections.)
 *
 * <p>One technique that significantly increases robustness is to treat the geodesic and projected
 * edges as parametric curves rather than geometric ones. Given a spherical edge AB and a projection
 * p:S2->R2, let f(t) be the normalized arc length parametrization of AB and let g(t) be the
 * normalized arc length parameterization of the projected edge p(A)p(B). (In other words, f(0)=A,
 * f(1)=B, g(0)=p(A), g(1)=p(B).) We now define the geometric error as the maximum distance from the
 * point p^-1(g(t)) to the geodesic edge AB for any t in [0,1], where p^-1 denotes the inverse
 * projection. In other words, the geometric error is the maximum distance from any point on the
 * projected edge (mapped back onto the sphere) to the geodesic edge AB. On the other hand we define
 * the parametric error as the maximum distance between the points f(t) and p^-1(g(t)) for any t in
 * [0,1], i.e. the maximum distance (measured on the sphere) between the geodesic and projected
 * points at the same interpolation fraction t.
 *
 * <p>The easiest way to estimate the parametric error is to simply evaluate both edges at their
 * midpoints and measure the distance between them (the "midpoint method"). This is very fast and
 * works quite well for most edges, however it has one major drawback: it doesn't handle points of
 * inflection (i.e., points where the curvature changes sign). For example, edges in the Mercator
 * and Plate Carree projections always curve towards the equator relative to the corresponding
 * geodesic edge, so in these projections there is a point of inflection whenever the projected edge
 * crosses the equator. The worst case occurs when the edge endpoints have different longitudes but
 * the same absolute latitude, since in that case the error is non-zero but the edges have exactly
 * the same midpoint (on the equator).
 *
 * <p>One solution to this problem is to split the input edges at all inflection points (i.e., along
 * the equator in the case of the Mercator and Plate Carree projections). However for general
 * projections these inflection points can occur anywhere on the sphere (e.g., consider the
 * Transverse Mercator projection). This could be addressed by adding methods to the S2Projection
 * interface to split edges at inflection points but this would make it harder and more error-prone
 * to implement new projections.
 *
 * <p>Another problem with this approach is that the midpoint method sometimes underestimates the
 * true error even when edges do not cross the equator. For the Plate Carree and Mercator
 * projections, the midpoint method can underestimate the error by up to 3%.
 *
 * <p>Both of these problems can be solved as follows. We assume that the error can be modeled as a
 * convex combination of two worst-case functions, one where the error is maximized at the edge
 * midpoint and another where the error is *minimized* (i.e., zero) at the edge midpoint. For
 * example, we could choose these functions as:
 *
 * <p>E1(x) = 1 - x^2 E2(x) = x * (1 - x^2)
 *
 * <p>where for convenience we use an interpolation parameter "x" in the range [-1, 1] rather than
 * the original "t" in the range [0, 1]. Note that both error functions must have roots at x = {-1,
 * 1} since the error must be zero at the edge endpoints. E1 is simply a parabola whose maximum
 * value is 1 attained at x = 0, while E2 is a cubic with an additional root at x = 0, and whose
 * maximum value is 2 * sqrt(3) / 9 attained at x = 1 / sqrt(3).
 *
 * <p>Next, it is convenient to scale these functions so that the both have a maximum value of 1. E1
 * already satisfies this requirement, and we simply redefine E2 as
 *
 * <p>E2(x) = x * (1 - x^2) / (2 * sqrt(3) / 9)
 *
 * <p>Now define x0 to be the point where these two functions intersect, i.e. the point in the range
 * (-1, 1) where E1(x0) = E2(x0). This value has the very convenient property that if we evaluate
 * the actual error E(x0), then the maximum error on the entire interval [-1, 1] is bounded by
 *
 * <p>E(x) <= E(x0) / E1(x0)
 *
 * <p>since whether the error is modeled using E1 or E2, the resulting function has the same maximum
 * value (namely E(x0) / E1(x0)). If it is modeled as some other convex combination of E1 and E2,
 * the maximum value can only decrease.
 *
 * <p>Finally, since E2 is not symmetric about the y-axis, we must also allow for the possibility
 * that the error is a convex combination of E1 and -E2. This can be handled by evaluating the error
 * at E(-x0) as well, and then computing the final error bound as
 *
 * <p>E(x) <= max(E(x0), E(-x0)) / E1(x0) .
 *
 * <p>Effectively, this method is simply evaluating the error at two points about 1/3 and 2/3 of the
 * way along the edges, and then scaling the maximum of these two errors by a constant factor.
 * Intuitively, the reason this works is that if the two edges cross somewhere in the interior, then
 * at least one of these points will be far from the crossing.
 *
 * <p>The actual algorithm implemented below has some additional refinements. First, edges longer
 * than 90 degrees are always subdivided; this avoids various unusual situations that can happen
 * with very long edges, and there is really no reason to avoid adding vertices to edges that are so
 * long.
 *
 * <p>Second, the error function E1 above needs to be modified to take into account spherical
 * distortions. (It turns out that spherical distortions are beneficial in the case of E2, i.e. they
 * only make its error estimates slightly more conservative.) To do this, we model E1 as the maximum
 * error in a Plate Carree edge of length 90 degrees or less. This turns out to be an edge from
 * 45:-90 to 45:90 (in lat:lng format). The corresponding error as a function of "x" in the range
 * [-1, 1] can be computed as the distance between the Plate Caree edge point (45, 90 * x) and the
 * geodesic edge point (90 - 45 * abs(x), 90 * sgn(x)). Using the Haversine formula, the
 * corresponding function E1 (normalized to have a maximum value of 1) is:
 *
 * <pre>{@code
 *   E1(x) =
 *     asin(sqrt(sin(Pi / 8 * (1 - x)) ^ 2 +
 *               sin(Pi / 4 * (1 - x)) ^ 2 * cos(Pi / 4) * sin(Pi / 4 * x))) /
 *     asin(sqrt((1 - 1 / sqrt(2)) / 2))
 * }</pre>
 *
 * <p>Note that this function does not need to be evaluated at runtime, it simply affects the
 * calculation of the value x0 where E1(x0) = E2(x0) and the corresponding scaling factor C = 1 /
 * E1(x0).
 *
 * <p>In the case of the Mercator and Plate Carree projections this strategy produces a conservative
 * upper bound (verified using 10 million random edges). Furthermore the bound is nearly tight; the
 * scaling constant is C = 1.19289, whereas the maximum observed value was 1.19254.
 *
 * <p>Compared to the simpler midpoint evaluation method, this strategy requires more function
 * evaluations (currently twice as many, but with a smarter tessellation algorithm it will only be
 * 50% more). It also results in a small amount of additional tessellation (about 1.5%) compared to
 * the midpoint method, but this is due almost entirely to the fact that the midpoint method does
 * not yield conservative error estimates.
 *
 * <p>For random edges with a tolerance of 1 meter, the expected amount of overtessellation is as
 * follows:
 * <ul>
 * <li>Plate Caree Midpoint: 1.8%
 * <li>Plate Caree Cubic: 3.0%
 * <li>Mercator Midpoint: 15.8%
 * <li>Mercator Cubic: 17.4%
 * </ul>
 */
@CheckReturnValue
final class S2EdgeTessellator {

  /**
   * The interpolation fraction at which the two edges are evaluated in order to measure the error
   * between them. (Edges are evaluated at two points measured this fraction from either end.) With
   * respect to the algorithm description above, this value is t0 = (1 - x0) / 2 in the range [0, 1]
   * that corresponds to x0 in the range [-1, 1] chosen such that E1(x0) == E2(x0).
   */
  private static final double INTERPOLATION_FRACTION = 0.31215691082248315;

  /** The following is the value of E1(x0) == E2(x0). */
  private static final double SCALE_FACTOR = 0.8382999256988851;

  /**
   * Returns the minimum supported tolerance (which corresponds to a distance less than one
   * micrometer on the Earth's surface). This is part of the API of S2EdgeTessellator.
   */
  public static final S1Angle MIN_TOLERANCE = S1Angle.radians(1e-13);

  private final Projection projection;
  private final S1ChordAngle scaledTolerance;

  /**
   * Constructs an S2EdgeTessellator using the given projection and error tolerance.
   *
   * <p>Method:
   * <li>AppendProjected- (input: S2 geodesics, output: Planar projected edges)
   * <li>AppendUnprojected- (input: Planar projected edges, output: S2 geodesics)
   */
  public S2EdgeTessellator(Projection projection, S1Angle tolerance) {
    this.projection = projection;
    this.scaledTolerance =
        S1ChordAngle.fromS1Angle(S1Angle.max(MIN_TOLERANCE, tolerance).mul(SCALE_FACTOR));
  }

  /**
   * Converts the spherical geodesic edge AB to a chain of planar edges in the given projection and
   * appends the corresponding vertices to "vertices".
   *
   * <p>This method can be called multiple times with the same output R2Vector list to convert an
   * entire polyline or loop. All vertices of the first edge are appended, but the first vertex of
   * each subsequent edge is omitted (and must match the last vertex of the previous edge).
   *
   * <p>If the given projection has one or more coordinate axes that "wrap", then every vertex's
   * coordinates will be as close as possible to the previous vertex's coordinates. Note that this
   * may yield vertices whose coordinates are outside the usual range. For example, tessellating the
   * edge (0:170, 0:-170) (in lat:lng notation) yields (0:170, 0:190).
   */
  public void appendProjected(S2Point a, S2Point b, List<R2Vector> vertices) {
    R2Vector pa = projection.project(a);
    if (vertices.isEmpty()) {
      vertices.add(pa);
    } else {
      pa = projection.wrapDestination(Iterables.getLast(vertices), pa);
    }
    R2Vector pb = projection.project(b);
    appendProjectedHelper(pa, a, pb, b, vertices::add);
  }

  /**
   * Converts the planar edge AB in the given projection to a chain of spherical geodesic edges and
   * appends the vertices to "vertices".
   *
   * <p>This method can be called multiple times with the same output S2Point list to convert an
   * entire polyline or loop. All vertices of the first edge are appended, but the first vertex of
   * each subsequent edge is omitted (and is required to match that last vertex of the previous
   * edge).
   *
   * <p>Note that to construct an S2Loop, you must remove the last S2Point at the very end to
   * eliminate the duplicate first and last vertex. Note also that if the given projection involves
   * coordinate "wrapping" (e.g. across the 180 degree meridian) then the first and last vertices
   * may not be exactly the same.
   */
  public void appendUnprojected(R2Vector pa, R2Vector pb, List<S2Point> vertices) {
    S2Point a = projection.unproject(pa);
    S2Point b = projection.unproject(pb);
    if (vertices.isEmpty()) {
      vertices.add(a);
    }
    appendUnprojectedHelper(pa, a, pb, b, vertices::add);
  }

  private void appendProjectedHelper(
      R2Vector pa, S2Point a, R2Vector pbIn, S2Point b, Consumer<R2Vector> vertexAdder) {
    R2Vector pb = projection.wrapDestination(pa, pbIn);
    if (estimateMaxError(pa, a, pb, b).lessOrEquals(scaledTolerance)) {
      vertexAdder.accept(pb);
    } else {
      S2Point mid = a.add(b).normalize();
      R2Vector projectedMid = projection.wrapDestination(pa, projection.project(mid));
      appendProjectedHelper(pa, a, projectedMid, mid, vertexAdder);
      appendProjectedHelper(projectedMid, mid, pb, b, vertexAdder);
    }
  }

  private void appendUnprojectedHelper(
      R2Vector pa, S2Point a, R2Vector pbIn, S2Point b, Consumer<S2Point> vertexAdder) {
    R2Vector pb = projection.wrapDestination(pa, pbIn);
    if (estimateMaxError(pa, a, pb, b).lessOrEquals(scaledTolerance)) {
      vertexAdder.accept(b);
    } else {
      R2Vector projectedMid = Projection.interpolate(0.5, pa, pb);
      S2Point mid = projection.unproject(projectedMid);
      appendUnprojectedHelper(pa, a, projectedMid, mid, vertexAdder);
      appendUnprojectedHelper(projectedMid, mid, pb, b, vertexAdder);
    }
  }

  private S1ChordAngle estimateMaxError(R2Vector pa, S2Point a, R2Vector pb, S2Point b) {
    if (a.dotProd(b) < -1e-14) {
      return S1ChordAngle.INFINITY;
    }

    double t1 = INTERPOLATION_FRACTION;
    double t2 = 1 - INTERPOLATION_FRACTION;
    S2Point mid1 = S2EdgeUtil.interpolate(t1, a, b);
    S2Point mid2 = S2EdgeUtil.interpolate(t2, a, b);
    S2Point projectedMid1 = projection.unproject(Projection.interpolate(t1, pa, pb));
    S2Point projectedMid2 = projection.unproject(Projection.interpolate(t2, pa, pb));
    S1ChordAngle mid1Angle = new S1ChordAngle(mid1, projectedMid1);
    S1ChordAngle mid2Angle = new S1ChordAngle(mid2, projectedMid2);
    return S1ChordAngle.max(mid1Angle, mid2Angle);
  }
}
