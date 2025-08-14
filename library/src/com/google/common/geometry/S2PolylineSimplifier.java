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

import static com.google.common.geometry.S2.M_PI_2;
import static com.google.common.geometry.S2Predicates.DBL_ERR;
import static java.lang.Math.PI;
import static java.lang.Math.asin;
import static java.lang.Math.atan2;
import static java.lang.Math.sqrt;

import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;

/**
 * This is a helper class for simplifying polylines. It allows you to compute a maximal edge that
 * intersects a sequence of discs, and that optionally avoids a different sequence of discs. The
 * results are conservative in that the edge is guaranteed to intersect or avoid the specified discs
 * using exact arithmetic (see S2Predicates.java).
 *
 * <p>Note that S2Builder can also simplify polylines and supports more features (e.g., snapping to
 * S2CellId centers), so it is only recommended to use this class if S2Builder does not meet your
 * needs.
 *
 * <p>Here is a simple example showing how to simplify a polyline into a sequence of edges that stay
 * within "maxError" of the original edges:
 *
 * {@snippet :
 * List<S2Point> v = { ... };
 * S2PolylineSimplifier simplifier = new S2PolylineSimplifier();
 * simplifier.init(v.get(0));
 * for (int i = 1; i < v.size(); ++i) {
 *   if (!simplifier.extend(v.get(i))) {
 *     outputEdge(simplifier.src(), v.get(i - 1));
 *     simplifier.init(v.get(i - 1));
 *   }
 *   simplifier.targetDisc(v.get(i), maxError);
 * }
 * outputEdge(simplifer.src(), v.back());
 *
 * }
 *
 * <p>Note that the points targeted by targetDisc do not need to be the same as the candidate
 * endpoints passed to extend(). So for example, you could target the original vertices of a
 * polyline, but only consider endpoints that are snapped to E7 coordinates or S2CellId centers.
 *
 * <p>Please be aware that this class works by maintaining a range of acceptable angles (bearings)
 * from the start vertex to the hypothetical destination vertex. It does not keep track of distances
 * to any of the discs to be targeted or avoided. Therefore to use this class correctly, constraints
 * should be added in increasing order of distance. (The actual requirement is slightly weaker than
 * this, which is why it is not enforced, but basically you should only call targetDisc() and
 * avoidDisc() with arguments that you want to constrain the immediately following call to
 * extend().)
 */
public class S2PolylineSimplifier {
  /** Output edge source vertex. */
  private S2Point src;

  /** First vector of an orthonormal frame for mapping vectors to angles. */
  private S2Point xDir;

  /** Second vector of an orthonormal frame for mapping vectors to angles. */
  private S2Point yDir;

  /** Allowable range of angles for the output edge. */
  private S1Interval window;

  /**
   * We store the discs to avoid individually until targetDisc() is first called with a disc that
   * does not contain the source vertex. At that time all such discs are processed by using them to
   * constrain "window", and this list is cleared.
   */
  private final ArrayList<RangeToAvoid> rangesToAvoid = new ArrayList<>();

  /** Constructs a new S2PolylineSimplifier. */
  public S2PolylineSimplifier() {}

  /** Starts a new simplified edge at "src". */
  public void init(S2Point src) {
    this.src = src;
    window = S1Interval.full();
    rangesToAvoid.clear();

    // Precompute basis vectors for the tangent space at "src". This is similar to getFrame() except
    // that we don't normalize the vectors. As it turns out, the two basis vectors below have the
    // same magnitude, up to the length error in S2Point.normalize().

    // Find the index of the component whose magnitude is smallest.
    S2Point tmp = src.fabs();
    int i = tmp.x < tmp.y ? (tmp.x < tmp.z ? 0 : 2) : (tmp.y < tmp.z ? 1 : 2);

    // Access the elements of 'src' by index.
    double[] s = new double[] {src.x, src.y, src.z};

    // Vectors for x and y tangent directions.
    double[] xVector = new double[3];
    double[] yVector = new double[3];

    // We define the "y" basis vector as the cross product of "src" and the basis vector for axis
    // "i". Let "j" and "k" be the indices of the other two components in cyclic order.
    int j = (i == 2 ? 0 : i + 1);
    int k = (i == 0 ? 2 : i - 1);
    yVector[i] = 0;
    yVector[j] = s[k];
    yVector[k] = -s[j];

    // Compute the cross product of "yDir" and "src". We write out the cross product here mainly for
    // documentation purposes; it also happens to save a few multiplies because unfortunately the
    // optimizer does *not* get rid of multiplies by zero (since these multiplies propagate NaN, for
    // example).
    xVector[i] = s[j] * s[j] + s[k] * s[k];
    xVector[j] = -s[j] * s[i];
    xVector[k] = -s[k] * s[i];

    xDir = new S2Point(xVector);
    yDir = new S2Point(yVector);
  }

  /**
   * Returns true if the edge (src, dst) satisfies all of the targeting requirements so far. Returns
   * false if the edge would be longer than 90 degrees (such edges are not supported) or cannot
   * satisfy the constraints.
   */
  public boolean extend(S2Point dst) {
    // We limit the maximum edge length to 90 degrees in order to simplify the error bounds. (The
    // error gets arbitrarily large as the edge length approaches 180 degrees.)
    S1ChordAngle edgeLength = new S1ChordAngle(src, dst);
    if (edgeLength.greaterThan(S1ChordAngle.RIGHT)) {
      return false;
    }

    // Otherwise check whether this vertex is in the acceptable angle range.
    double dir = getDirection(dst);
    if (!window.contains(dir)) {
      return false;
    }

    // Also check any angles ranges to avoid that have not been processed yet.
    for (RangeToAvoid range : rangesToAvoid) {
      if (range.interval.contains(dir)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Requires that the output edge must pass through the disc specified by the given point and
   * radius.
   *
   * <p>Returns true if it is possible to intersect the target disc, given previous constraints.
   */
  @CanIgnoreReturnValue
  public boolean targetDisc(S2Point point, S1ChordAngle radius) {
    // Shrink the target interval by the maximum error from all sources. This guarantees that the
    // output edge will intersect the given disc.
    double semiwidth = getSemiwidth(point, radius, -1 /* round down */);
    if (semiwidth >= PI) {
      // The target disc contains "src", so there is nothing to do.
      return true;
    }
    if (semiwidth < 0) {
      window = S1Interval.empty();
      return false;
    }

    // Otherwise compute the angle interval corresponding to the target disc and intersect it with
    // the current window.
    double center = getDirection(point);
    S1Interval target = S1Interval.fromPoint(center).expanded(semiwidth);
    window = window.intersection(target);

    // If there are any angle ranges to avoid, they can be processed now.
    for (RangeToAvoid range : rangesToAvoid) {
      avoidRange(range.interval, range.onLeft);
    }
    rangesToAvoid.clear();

    return !window.isEmpty();
  }

  /**
   * Requires that the output edge must avoid the given disc. "discOnLeft" specifies whether the
   * disc must be to the left or right of the output edge AB. (This feature allows the simplified
   * edge to preserve the topology of the original polyline with respect to other nearby points.)
   *
   * <p>More precisely, let AB be the output edge, P be the center of the disc, and r be its radius.
   *
   * <pre>Then this method ensures that {@code
   *   (1) distance(AB, P) > r, and
   *   (2) if dotProd(AB, AP) > 0, then sign(ABP) > 0 iff discOnLeft is true.
   * }
   *
   * <p>The second condition says that "discOnLeft" has an effect if and only if P is not behind the
   * source vertex A with respect to the direction AB.
   *
   * <p>If your input is a polyline, you can compute "discOnLeft" as follows. Let the polyline be
   * ABCDE and assume that it already avoids a set of points X_i. Suppose that you have already
   * added ABC to the simplifier, and now want to extend the edge chain to D. First find the X_i
   * that are near the edge CD, then discard the ones such that {@code AX_i <= AC} or {@code AX_i >= AD}
   * (since these points have either already been considered or aren't relevant yet). Now X_i is to
   * the left of the polyline if and only if S2Predicates.orderedCCW(A, D, X_i, C) (in other words,
   * if X_i is to the left of the angle wedge ACD). Note that simply testing S2Predicates.sign(C, D,
   * X_i) or S2Predicates.sign(A, D, X_i) does not handle all cases correctly.
   *
   * <p>Returns true if the disc can be avoided given previous constraints, or if the discs to avoid
   * have not been processed yet. Returns false if the disc cannot be avoided.
   */
  @CanIgnoreReturnValue
  public boolean avoidDisc(S2Point point, S1ChordAngle radius, boolean discOnLeft) {
    // Expand the interval by the maximum error from all sources. This guarantees that the final
    // output edge will avoid the given disc.
    double semiwidth = getSemiwidth(point, radius, 1 /*round up*/);
    if (semiwidth >= PI) {
      // The disc to avoid contains "src", so it can't be avoided.
      window = S1Interval.empty();
      return false;
    }
    // Compute the disallowed range of angles: the angle subtended by the disc on one side, and 90
    // degrees on the other (to satisfy "discOnLeft").
    double center = getDirection(point);
    double dLeft = discOnLeft ? M_PI_2 : semiwidth;
    double dRight = discOnLeft ? semiwidth : M_PI_2;
    S1Interval avoidInterval =
        new S1Interval(
            Platform.IEEEremainder(center - dRight, 2 * PI),
            Platform.IEEEremainder(center + dLeft, 2 * PI));

    if (window.isFull()) {
      // Discs to avoid can't be processed until window is reduced to at most 180 degrees by a call
      // to TargetDisc(). Save it for later.
      rangesToAvoid.add(new RangeToAvoid(avoidInterval, discOnLeft));
      return true;
    }
    avoidRange(avoidInterval, discOnLeft);
    return !window.isEmpty();
  }

  private double getDirection(S2Point p) {
    return atan2(p.dotProd(yDir), p.dotProd(xDir));
  }

  /**
   * Computes half the angle in radians subtended from the source vertex by a disc of radius "r"
   * centered at "p", rounding the result conservatively up or down according to whether
   * roundDirection is +1 or -1. (So for example, if roundDirection == +1 then the return value is
   * an upper bound on the true result.)
   */
  private double getSemiwidth(S2Point p, S1ChordAngle r, int roundDirection) {
    // Using spherical trigonometry,
    //
    //   sin(semiwidth) = sin(r) / sin(a)
    //
    // where "a" is the angle between "src" and "p". Rather than measuring these angles, instead we
    // measure the squared chord lengths through the interior of the sphere (i.e., Cartersian
    // distance). Letting "r2" be the squared chord distance corresponding to "r", and "a2" be the
    // squared chord distance corresponding to "a", we use the relationships
    //
    //    sin^2(r) = r2 (1 - r2 / 4)
    //    sin^2(a) = d2 (1 - d2 / 4)
    //
    // which follow from the fact that r2 = (2 * sin(r / 2)) ^ 2, etc.

    // "a2" has a relative error up to 5 * DBL_ERR, plus an absolute error of up to 64 * DBL_ERR^2
    // (because "src" and "p" may differ from unit length by up to 4 * DBL_ERR). We can correct for
    // the relative error later, but for the absolute error we use "roundDirection" to account for
    // it now.
    double r2 = r.getLength2();
    double a2 = new S1ChordAngle(src, p).getLength2();
    a2 -= 64 * DBL_ERR * DBL_ERR * roundDirection;
    if (a2 <= r2) {
      return PI; // The given disc contains "src".
    }

    double sin2R = r2 * (1 - 0.25 * r2);
    double sin2A = a2 * (1 - 0.25 * a2);
    double semiwidth = asin(sqrt(sin2R / sin2A));

    // We compute bounds on the errors from all sources:
    //
    //   - The call to getSemiwidth (this call).
    //   - The call to getDirection that computes the center of the interval.
    //   - The call to getDirection in Extend that tests whether a given point
    //     is an acceptable destination vertex.
    //
    // Summary of the errors in getDirection:
    //
    // yDir has no error.
    //
    // xDir has a relative error of DBL_ERR in two components, a relative error of 2 * DBL_ERR in
    // the other component, plus an overall relative length error of 4 * DBL_ERR (compared to yDir)
    // because "src" is assumed to be normalized only to within the tolerances of
    // S2Point::Normalize().
    //
    // p.DotProd(yDir) has a relative error of 1.5 * DBL_ERR and an absolute error of 1.5 * DBL_ERR
    // * yDir.Norm().
    //
    // p.DotProd(xDir) has a relative error of 5.5 * DBL_ERR and an absolute error of 3.5 * DBL_ERR
    // * yDir.Norm() (noting that xDir and yDir have the same length to within a relative error of 4
    // * DBL_ERR).
    //
    // It's possible to show by taking derivatives that these errors can affect the angle atan2(y,
    // x) by up 7.093 * DBL_ERR radians. Rounding up and including the call to atan2 gives a final
    // error bound of 10 * DBL_ERR.
    //
    // Summary of the errors in getSemiwidth:
    //
    // The distance a2 has a relative error of 5 * DBL_ERR plus an absolute error of 64 * DBL_ERR^2
    // because the points "src" and "p" may differ from unit length (by up to 4 * DBL_ERR). We have
    // already accounted for the absolute error above, leaving only the relative error.
    //
    // sin2R has a relative error of 2 * DBL_ERR.
    //
    // sin2A has a relative error of 12 * DBL_ERR assuming that a2 <= 2, i.e. distance(src, p) <=
    // 90 degrees. (The relative error gets arbitrarily larger as this distance approaches 180
    // degrees.)
    //
    // semiwidth has a relative error of 17 * DBL_ERR.
    //
    // Finally, (center +/- semiwidth) has a rounding error of up to 4 * DBL_ERR because in theory,
    // the result magnitude may be as large as 1.5 * PI which is larger than 4.0. This gives a total
    // error of:
    double error = (2 * 10 + 4) * DBL_ERR + 17 * DBL_ERR * semiwidth;
    return semiwidth + roundDirection * error;
  }

  private void avoidRange(S1Interval avoidInterval, boolean discOnLeft) {
    // If "avoidInterval" is a proper subset of "window", then in theory the result should be two
    // intervals. One interval points towards the given disc and passes on the correct side of it,
    // while the other interval points away from the disc. However the latter interval never
    // contains an acceptable output edge direction (as long as this class is being used correctly)
    // and can be safely ignored. This is true because (1) "window" is not full, which means that it
    // contains at least one vertex of the input polyline and is at most 180 degrees in length, and
    // (2) "discOnLeft" is computed with respect to the next edge of the input polyline, which means
    // that the next input vertex is either inside "avoidInterval" or somewhere in the 180 degrees
    // to its right/left according to "discOnLeft", which means that it cannot be contained by the
    // subinterval that we ignore.
    assert !window.isFull();
    if (window.contains(avoidInterval)) {
      if (discOnLeft) {
        window = new S1Interval(window.lo(), avoidInterval.lo());
      } else {
        window = new S1Interval(avoidInterval.hi(), window.hi());
      }
    } else {
      window = window.intersection(avoidInterval.complement());
    }
  }

  /**
   * Unfortunately, the discs to avoid cannot be processed until the direction of the output edge is
   * constrained to lie within an S1Interval of at most 180 degrees. This happens only when the
   * first target disc is added that does not contain the source vertex. Until that time we simply
   * store all the discs as ranges of directions to avoid.
   */
  private static class RangeToAvoid {
    /** Range of directions to avoid. */
    final S1Interval interval;

    /** Is this disc to the left of the output edge? */
    final boolean onLeft;

    /** Constructs a new RangeToAvoid containing the given values. */
    public RangeToAvoid(S1Interval avoidInterval, boolean discOnLeft) {
      interval = avoidInterval;
      onLeft = discOnLeft;
    }
  }
}
