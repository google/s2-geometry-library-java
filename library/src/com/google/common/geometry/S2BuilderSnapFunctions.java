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

import static com.google.common.geometry.S2.DBL_EPSILON;
import static com.google.common.geometry.S2.M_SQRT1_2;
import static com.google.common.geometry.S2.M_SQRT2;
import static java.lang.Math.ceil;
import static java.lang.Math.log10;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.round;

import com.google.common.base.Preconditions;

/** A library of useful implementations of {@link S2Builder.SnapFunction}. */
public final class S2BuilderSnapFunctions {
  private S2BuilderSnapFunctions() {}

  /*** The maximum supported snap radius (equivalent to about 7800km). */
  public static S1Angle maxSnapRadius() {
    // This value can't be larger than 85.7 degrees without changing the code related to
    // minEdgeLengthToSplitCA, and increasing it to 90 degrees or more would most likely
    // require significant changes to the algorithm.
    return S1Angle.degrees(70);
  }

  /**
   * A SnapFunction that snaps every vertex to itself. It should be used when vertices do not need
   * to be snapped to a discrete set of locations (such as E7 lat/lngs), or when maximum accuracy is
   * desired.
   *
   * <p>If the given "snapRadius" is zero, then all input vertices are preserved exactly. Otherwise,
   * S2Builder merges nearby vertices to ensure that no vertex pair is closer than "snapRadius".
   * Furthermore, vertices are separated from non-incident edges by at least
   * "minEdgeVertexSeparation", equal to (0.5 * snapRadius). For example, if the snapRadius is 1km,
   * then vertices will be separated from non-incident edges by at least 500m.
   */
  public static class IdentitySnapFunction implements S2Builder.SnapFunction {
    private final S1Angle snapRadius;

    /**
     * Returns a default {@link IdentitySnapFunction} that snaps every vertex to itself. Sets a
     * snapRadius of zero (i.e., no snapping).
     */
    public IdentitySnapFunction() {
      snapRadius = S1Angle.ZERO;
    }

    /**
     * Returns an {@link IdentitySnapFunction} that snaps every vertex to itself, aside from merging
     * vertices and ensuring vertex/edge separation as required for the given snapRadius.
     *
     * <p>REQUIRES: {@code snapRadius <= maxSnapRadius()}.
     */
    public IdentitySnapFunction(S1Angle snapRadius) {
      Preconditions.checkArgument(
          snapRadius.lessOrEquals(maxSnapRadius()),
          "snapRadius of %s exceeds maximum of %s", snapRadius, maxSnapRadius());
      this.snapRadius = snapRadius;
    }

    @Override
    public S1Angle snapRadius() {
      return snapRadius;
    }

    /** For the identity snap function, all vertex pairs are separated by at least snapRadius(). */
    @Override
    public S1Angle minVertexSeparation() {
      // Since SnapFunction does not move input points, output vertices are separated by the full
      // snapRadius().
      return snapRadius;
    }

    /**
     * For the identity snap function, edges are separated from all non-incident vertices by at
     * least 0.5 * snapRadius().
     */
    @Override
    public S1Angle minEdgeVertexSeparation() {
      // In the worst case configuration, the edge-vertex separation is half of the vertex
      // separation.
      return snapRadius.mul(0.5f);
    }

    @Override
    public S2Point snapPoint(S2Point point) {
      return point;
    }

    @Override
    public String toString() {
      return "IdentitySnapFunction with snapRadius " + snapRadius;
    }
  }

  /**
   * A SnapFunction that snaps vertices to S2CellId centers. This can be useful if you want to
   * encode your geometry compactly using {@link S2Polygon#encode(OutputStream)}, for example. You
   * can snap to the centers of cells at any level.
   *
   * <p>Every snap level has a corresponding minimum snap radius, which is simply the maximum
   * distance that a vertex can move when snapped. It is approximately equal to half of the maximum
   * diagonal length for cells at the chosen level. You can also set the snap radius to a larger
   * value; for example, you could snap to the centers of leaf cells (1cm resolution) but set the
   * snapRadius() to 10m. This would result in significant extra simplification, without moving
   * vertices unnecessarily (i.e., vertices that are at least 10m away from all other vertices will
   * move by less than 1cm).
   */
  public static class S2CellIdSnapFunction implements S2Builder.SnapFunction {
    private final int level;
    private final S1Angle snapRadius;

    /**
     * Constructor that defaults to snapping to S2CellId.MAX_LEVEL (i.e., the centers of leaf
     * cells), and uses the minimum allowable snap radius at that level.
     */
    public S2CellIdSnapFunction() {
      this.level = S2CellId.MAX_LEVEL;
      this.snapRadius = minSnapRadiusForLevel(this.level);
    }

    /**
     * Creates an S2CellIdSnapFunction that snaps vertices to S2Cell centers at the given level.
     * This method also sets "snapRadius" to the minimum value allowed at this level,
     * {@code minSnapRadiusForLevel(level) }. If you want to use a larger snap radius than the
     * minimum, (e.g., to simplify the geometry), use the constructor below that also accepts a
     * specified snapRadius.
     */
    public S2CellIdSnapFunction(int level) {
      Preconditions.checkArgument(level >= 0 && level <= S2CellId.MAX_LEVEL);
      this.level = level;
      this.snapRadius = minSnapRadiusForLevel(level);
    }

    /**
     * Creates an S2CellIdSnapFunction that snaps vertices to S2Cell centers at the given level,
     * and sets "snapRadius" to the given value. See {@link S2Builder} for more on the snap radius.
     * The snap radius must be at least the minimum value for the given level, and at most the
     * maximum snap radius.
     */
    public S2CellIdSnapFunction(int level, S1Angle snapRadius) {
      Preconditions.checkArgument(level >= 0 && level <= S2CellId.MAX_LEVEL);
      Preconditions.checkArgument(snapRadius.lessOrEquals(minSnapRadiusForLevel(level)));
      Preconditions.checkArgument(snapRadius.lessOrEquals(maxSnapRadius()));
      this.level = level;
      this.snapRadius = snapRadius;
    }

    /** Returns the S2 cell level that vertices will be snapped to. */
    public int level() {
      return level;
    }

    /** Returns the snap radius. */
    @Override
    public S1Angle snapRadius() {
      return snapRadius;
    }

    /**
     * Returns the minimum allowable snap radius for the given S2Cell level (approximately equal to
     * half of the maximum cell diagonal length).
     */
    public static S1Angle minSnapRadiusForLevel(int level) {
      // snapRadius() needs to be an upper bound on the true distance that a point can move when
      // snapped, taking into account numerical errors.
      //
      // The maximum error when converting from an S2Point to an S2CellId is
      // {@code S2Projections.PROJ.maxDiag.deriv() * DBL_EPSILON}. The maximum error when converting
      // an S2CellId center back to an S2Point is 1.5 * DBL_EPSILON. These add up to just slightly
      // less than 4 * DBL_EPSILON.
      return S1Angle.radians(0.5 * S2Projections.MAX_DIAG.getValue(level) + 4 * DBL_EPSILON);
    }

    /**
     * Returns the minimum S2Cell level (i.e., largest S2Cells) such that vertices will not move by
     * more than "snapRadius". This can be useful when choosing an appropriate level to snap to. The
     * return value is always a valid level (out of range values are silently clamped).
     *
     * <p>If you want to choose the snap level based on a distance, and then use the minimum
     * possible snap radius for the chosen level, do this:
     *
     * {@snippet :
     * S2CellIdSnapFunction f = new S2CellIdSnapFunction(
     *     S2CellIdSnapFunction.levelForMaxSnapRadius(distance));
     * }
     */
    public static int levelForMaxSnapRadius(S1Angle snapRadius) {
      // When choosing a level, we need to account for the error bound of 4 * DBL_EPSILON that is
      // added by minSnapRadiusForLevel().
      return S2Projections.MAX_DIAG.getMinLevel(2 * (snapRadius.radians() - 4 * DBL_EPSILON));
    }

    /**
     * For S2CellId snapping, the minimum separation between vertices depends on level() and
     * snapRadius(). It can vary between 0.5 * snapRadius() and snapRadius().
     */
    @Override
    public S1Angle minVertexSeparation() {
      // We have three different bounds for the minimum vertex separation: one is a constant bound,
      // one is proportional to snapRadius, and one is equal to snapRadius minus a constant. These
      // bounds give the best results for small, medium, and large snap radii respectively. We
      // return the maximum of the three bounds.
      //
      // 1. Constant bound: Vertices are always separated by at least minEdge(level), the minimum
      //    edge length for the chosen snap level.
      //
      // 2. Proportional bound: It can be shown that in the plane, the worst-case configuration has
      //    a vertex separation of 2 / sqrt(13) * snapRadius. This is verified in the unit test,
      //    except that on the sphere the ratio is slightly smaller at cell level 2 (0.54849 vs.
      //    0.55470). We reduce that value a bit more below to be conservative.
      //
      // 3. Best asymptotic bound: This bound bound is derived by observing we only select a new
      //    site when it is at least snapRadius() away from all existing sites, and the site can
      //    move by at most 0.5 * maxDiag(level) when snapped.
      S1Angle minEdge = S1Angle.radians(S2Projections.MIN_EDGE.getValue(level));
      S1Angle maxDiag = S1Angle.radians(S2Projections.MAX_DIAG.getValue(level));
      return S1Angle.max(
          minEdge,
          S1Angle.max(
              snapRadius.mul(0.548), // 2 / sqrt(13) in the plane
              snapRadius.sub(maxDiag.mul(0.5))));
    }

    /**
     * For S2CellId snapping, the minimum separation between edges and non-incident vertices depends
     * on level() and snapRadius(). It can be as low as 0.219 * snapRadius(), but is typically 0.5 *
     * snapRadius() or more.
     */
    @Override
    public S1Angle minEdgeVertexSeparation() {
      // Similar to minVertexSeparation(), in this case we have four bounds: a constant bound that
      // holds only at the minimum snap radius, a constant bound that holds for any snap radius, a
      // bound that is proportional to snapRadius, and a bound that approaches 0.5 * snapRadius
      // asymptotically.
      //
      // 1. Constant bounds:
      //
      //    (a) At the minimum snap radius for a given level, it can be shown that vertices are
      //    separated from edges by at least 0.5 * minDiag(level) in the plane. The unit test
      //    verifies this, except that on the sphere the worst case is slightly better:
      //    0.5652980068 * minDiag(level).
      //
      //    (b) Otherwise, for arbitrary snap radii the worst-case configuration in the plane has an
      //    edge-vertex separation of sqrt(3/19) * kMinDiag(level), where sqrt(3/19) is about
      //    0.3973597071. The unit test verifies that the bound is slightly better on the sphere:
      //    0.3973595687 * kMinDiag(level).
      //
      // 2. Proportional bound: In the plane, the worst-case configuration has an edge-vertex
      //    separation of 2 * sqrt(3/247) * snapRadius, which is about 0.2204155075. The unit test
      //    verifies this, except that on the sphere the bound is slightly worse for certain large
      //    S2Cells: the minimum ratio occurs at cell level 6, and is about 0.2196666953.
      //
      // 3. Best asymptotic bound: If snapRadius() is large compared to the minimum snap radius,
      //    then the best bound is achieved by 3 sites on a circular arc of radius "snapRadius",
      //    spaced "minVertexSeparation" apart. An input edge passing just to one side of the center
      //    of the circle intersects the Voronoi regions of the two end sites but not the Voronoi
      //    region of the center site, and gives an edge separation of (minVertexSeparation ** 2) /
      //    (2 * snapRadius). This bound approaches 0.5 * snapRadius for large snap radii, i.e. the
      //    minimum edge-vertex separation approaches half of the minimum vertex separation as the
      //    snap radius becomes large compared to the cell size.
      S1Angle minDiag = S1Angle.radians(S2Projections.MIN_DIAG.getValue(level));
      if (snapRadius().equals(minSnapRadiusForLevel(level))) {
        // This bound only holds when the minimum snap radius is being used.
        return minDiag.mul(0.565); // 0.500 in the plane
      }
      // Otherwise, these bounds hold for any snapRadius().
      S1Angle vertexSep = minVertexSeparation();
      return S1Angle.max(
          minDiag.mul(0.397), // sqrt(3 / 19) in the plane
          S1Angle.max(
              snapRadius.mul(0.219), // 2 * sqrt(3 / 247) in the plane
              vertexSep.mul(0.5 * vertexSep.div(snapRadius))));
    }

    @Override
    public S2Point snapPoint(S2Point point) {
      return S2CellId.fromPoint(point).parent(level).toPoint();
    }

    @Override
    public String toString() {
      return "S2CellIdSnapFunction with level " + level + " and snapRadius " + snapRadius;
    }
  }

  /**
   * A SnapFunction that snaps vertices to S2LatLng E5, E6, or E7 coordinates. These coordinates are
   * expressed in degrees multiplied by a power of 10 and then rounded to the nearest integer. For
   * example, in E6 coordinates the point (23.12345651, -45.65432149) would become (23123457,
   * -45654321).
   *
   * <p>The main argument of the SnapFunction is the exponent for the power of 10 that coordinates
   * should be multiplied by before rounding. For example, IntLatLngSnapFunction(7) is a function
   * that snaps to E7 coordinates. The exponent can range from 0 to 10.
   *
   * <p>Each exponent has a corresponding minimum snap radius, which is simply the maximum distance
   * that a vertex can move when snapped. It is approximately equal to 1/sqrt(2) times the nominal
   * point spacing; for example, for snapping to E7 the minimum snap radius is (1e-7 / sqrt(2))
   * degrees.
   *
   * <p>You can also set the snap radius to any value larger than this; this can result in
   * significant extra simplification (similar to using a larger exponent) but does not move
   * vertices unnecessarily.
   */
  public static class IntLatLngSnapFunction implements S2Builder.SnapFunction {

    /** The minimum exponent supported for snapping. */
    public static final int MIN_EXPONENT = 0;

    /** The maximum exponent supported for snapping. */
    public static final int MAX_EXPONENT = 10;

    private final int exponent;
    private final S1Angle snapRadius;
    private final double fromDegrees;
    private final double toDegrees;

    /**
     * Constructs an IntLatLngSnapFunction that snaps vertices to points whose (lat, lng)
     * coordinates are integers after converting to degrees and multiplying by 10 raised to the
     * given exponent. For example, (exponent == 7) yields E7 coordinates. Sets "snapRadius" to the
     * minimum value allowed for this exponent.
     *
     * <p>You may use a larger snap radius than the minimum (e.g., to simplify the geometry) by
     * using the constructor below which accepts a specific snapRadius.
     *
     * <p>REQUIRES: {@code MIN_EXPONENT <= exponent <= MAX_EXPONENT}.
     */
    public IntLatLngSnapFunction(int exponent) {
      Preconditions.checkArgument(exponent >= MIN_EXPONENT);
      Preconditions.checkArgument(exponent <= MAX_EXPONENT);
      this.exponent = exponent;
      this.snapRadius = minSnapRadiusForExponent(exponent);

      // Precompute the scale factors needed for snapping. Note that these calculations need to
      // exactly match the ones in S1Angle to ensure that the same S2Points are generated.
      double power = Math.pow(10, exponent);
      this.fromDegrees = power;
      this.toDegrees = 1 / power;
    }

    /**
     * Constructs an IntLatLngSnapFunction that snaps vertices to points whose (lat, lng)
     * coordinates are integers after converting to degrees and multiplying by 10 raised to the
     * given exponent. For example, (exponent == 7) yields E7 coordinates. Sets "snapRadius" to the
     * provided value, which must be at least the minimum valid value for the exponent, and at most
     * the maximum snap radius.
     *
     * <p>REQUIRES: {@code MIN_EXPONENT <= exponent <= MAX_EXPONENT}.
     */
    public IntLatLngSnapFunction(int exponent, S1Angle snapRadius) {
      Preconditions.checkArgument(exponent >= MIN_EXPONENT);
      Preconditions.checkArgument(exponent <= MAX_EXPONENT);
      Preconditions.checkState(snapRadius.greaterOrEquals(minSnapRadiusForExponent(exponent)));
      Preconditions.checkState(snapRadius.lessOrEquals(maxSnapRadius()));
      this.exponent = exponent;
      this.snapRadius = snapRadius;

      // Precompute the scale factors needed for snapping. Note that these calculations need to
      // exactly match the ones in S1Angle to ensure that the same S2Points are generated.
      double power = Math.pow(10, exponent);
      this.fromDegrees = power;
      this.toDegrees = 1 / power;
    }

    /** Returns the exponent that vertices will be snapped to. */
    public int exponent() {
      return exponent;
    }

    @Override
    public S1Angle snapRadius() {
      return snapRadius;
    }

    /**
     * Returns the minimum allowable snap radius for the given exponent (approximately equal to
     * (pow(10, -exponent) / sqrt(2)) degrees).
     */
    public static S1Angle minSnapRadiusForExponent(int exponent) {
      // snapRadius() needs to be an upper bound on the true distance that a point can move when
      // snapped, taking into account numerical errors.
      //
      // The maximum errors in latitude and longitude can be bounded as follows (as absolute errors
      // in terms of DBL_EPSILON):
      //
      //                                      Latitude      Longitude
      // Convert to S2LatLng:                    1.000          1.000
      // Convert to degrees:                     1.032          2.063
      // Scale by 10**exp:                       0.786          1.571
      // Round to integer: 0.5 * S1Angle.degrees(toDegrees)
      // Scale by 10**(-exp):                    1.375          2.749
      // Convert to radians:                     1.252          1.503
      // ------------------------------------------------------------
      // Total (except for rounding)             5.445          8.886
      //
      // The maximum error when converting the S2LatLng back to an S2Point is
      //
      //   sqrt(2) * (maximum error in latitude or longitude) + 1.5 * DBL_EPSILON
      //
      // which works out to (9 * sqrt(2) + 1.5) * DBL_EPSILON radians. Finally we need to consider
      // the effect of rounding to integer coordinates (much larger than the errors above), which
      // can change the position by up to (sqrt(2) * 0.5 * toDegrees) radians.
      double power = 1;
      for (int i = 0; i < exponent; ++i) {
        power *= 10;
      }
      return S1Angle.degrees(M_SQRT1_2 / power)
          .add(S1Angle.radians((9 * M_SQRT2 + 1.5) * DBL_EPSILON));
    }

    /**
     * Returns the minimum exponent such that vertices will not move by more than "snapRadius". This
     * can be useful when choosing an appropriate exponent for snapping. The return value is always
     * a valid exponent (out of range values are silently clamped).
     *
     * <p>If you want to choose the exponent based on a distance, and then use the minimum possible
     * snap radius for that exponent, do this:
     *
     * {@snippet :
     * IntLatLngSnapFunction f = new IntLatLngSnapFunction(
     *       IntLatLngSnapFunction.exponentForMaxSnapRadius(distance));
     * }
     */
    public static int exponentForMaxSnapRadius(S1Angle snapRadius) {
      // When choosing an exponent, we need to account for the error bound of (9 * sqrt(2) + 1.5) *
      // DBL_EPSILON added by minSnapRadiusForExponent().
      snapRadius = snapRadius.sub(S1Angle.radians((9 * M_SQRT2 + 1.5) * DBL_EPSILON));
      snapRadius = S1Angle.max(snapRadius, S1Angle.radians(1e-30));
      double exponent = log10(M_SQRT1_2 / snapRadius.degrees());

      // There can be small errors in the calculation above, so to ensure that this function is the
      // inverse of MinSnapRadiusForExponent() we subtract a small error tolerance.
      return max(MIN_EXPONENT, min(MAX_EXPONENT, (int) ceil(exponent - 2 * DBL_EPSILON)));
    }

    /**
     * For IntLatLng snapping, the minimum separation between vertices depends on exponent() and
     * snapRadius(). It can vary between snapRadius() and snapRadius(). WHAT?
     */
    @Override
    public S1Angle minVertexSeparation() {
      // We have two bounds for the minimum vertex separation: one is proportional to snapRadius,
      // and one is equal to snapRadius minus a constant. These bounds give the best results for
      // small and large snap radii respectively. We return the maximum of the two bounds.
      //
      // 1. Proportional bound: It can be shown that in the plane, the worst-case configuration has
      //    a vertex separation of (sqrt(2) / 3) * snapRadius. This is verified in the unit test,
      //    except that on the sphere the ratio is slightly smaller (0.471337 vs. 0.471404). We
      //    reduce that value a bit more below to be conservative.
      //
      // 2. Best asymptotic bound: This bound bound is derived by observing we only select a new
      //    site when it is at least snapRadius() away from all existing sites, and snapping a
      //    vertex can move it by up to ((1 / sqrt(2)) * toDegrees) degrees.
      return S1Angle.max(
          snapRadius.mul(0.471), // sqrt(2) / 3 in the plane
          snapRadius.sub(S1Angle.degrees(M_SQRT1_2 * toDegrees)));
    }

    /**
     * For IntLatLng snapping, the minimum separation between edges and non-incident vertices
     * depends on level() and snapRadius(). It can be as low as 0.222 * snapRadius(), but is
     * typically 0.39 * snapRadius() or more.
     */
    @Override
    public S1Angle minEdgeVertexSeparation() {
      // Similar to minVertexSeparation(), in this case we have three bounds: one is a constant
      // bound, one is proportional to snapRadius, and one approaches 0.5 * snapRadius
      // asymptotically.
      //
      // 1. Constant bound: In the plane, the worst-case configuration has an edge-vertex separation
      //    of ((1 / sqrt(13)) * toDegrees) degrees. The unit test verifies this, except that on the
      //    sphere the ratio is slightly lower when small exponents such as E1 are used (0.2772589
      //    vs. 0.2773501).
      //
      // 2. Proportional bound: In the plane, the worst-case configuration has an edge-vertex
      //    separation of (2 / 9) * snapRadius (0.222222222222). The unit test verifies this, except
      //    that on the sphere the bound can be slightly worse with large exponents (e.g., E9) due
      //    to small numerical errors (0.222222126756717).
      //
      // 3. Best asymptotic bound: If snapRadius() is large compared to the minimum snap radius,
      //    then the best bound is achieved by 3 sites on a circular arc of radius "snapRadius",
      //    spaced "minVertexSeparation" apart (see S2CellIdSnapFunction.minEdgeVertexSeparation).
      //    This bound approaches 0.5 * snapRadius as the snap radius becomes large relative to the
      // grid spacing.
      S1Angle vertexSep = minVertexSeparation();
      return S1Angle.max(
          S1Angle.degrees(toDegrees).mul(0.277), // 1/sqrt(13) in the plane
          S1Angle.max(
              snapRadius.mul(0.222), // 2/9 in the plane
              vertexSep.mul(0.5 * vertexSep.div(snapRadius))));
    }

    @Override
    public S2Point snapPoint(S2Point point) {
      Preconditions.checkState(exponent >= 0); // Make sure the snap function was initialized.
      S2LatLng input = new S2LatLng(point);
      long lat = round(input.lat().degrees() * fromDegrees);
      long lng = round(input.lng().degrees() * fromDegrees);
      return S2LatLng.fromDegrees(lat * toDegrees, lng * toDegrees).toPoint();
    }

    @Override
    public String toString() {
      return "IntLatLngSnapFunction with exponent " + exponent + " and snapRadius " + snapRadius;
    }
  }
}
