/*
 * Copyright 2023 Google Inc.
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

import static com.google.common.geometry.S2.M_SQRT1_2;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.min;
import static java.lang.Math.round;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;

import com.google.common.base.Preconditions;
import com.google.common.geometry.S2BuilderSnapFunctions.IntLatLngSnapFunction;
import com.google.common.geometry.S2BuilderSnapFunctions.S2CellIdSnapFunction;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for S2CellIdSnapFunction and IntLatLngSnapFunction in S2BuilderSnapFunctions. */
@RunWith(JUnit4.class)
public final class S2BuilderSnapFunctionsTest extends GeometryTestCase {
  /**
   * For tests that check ratios like (minVertexSep / snapRadius), this is the tolerance for changes
   * between the ratios computed during the current test run, and the expected ratios, which were
   * saved from a previous run.
   */
  private static final double RATIO_TOLERANCE = 1e-7;

  private static final S2CellId SEARCH_ROOT_ID = S2CellId.fromFace(0);
  private static final S2CellId SEARCH_FOCUS_ID = S2CellId.fromFace(0).child(3);

  /**
   * Tests S2CellIdSnapFunction.minSnapRadiusForLevel and levelForMaxSnapRadius for consistency over
   * all possible levels.
   */
  @Test
  public void testS2CellIdSnapFunctionLevelToFromSnapRadius() {
    for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
      S1Angle radius = S2CellIdSnapFunction.minSnapRadiusForLevel(level);
      assertEquals(level, S2CellIdSnapFunction.levelForMaxSnapRadius(radius));
      assertEquals(
          min(level + 1, S2CellId.MAX_LEVEL),
          S2CellIdSnapFunction.levelForMaxSnapRadius(radius.mul(0.999)));
    }
    assertEquals(0, S2CellIdSnapFunction.levelForMaxSnapRadius(S1Angle.radians(5)));
    assertEquals(
        S2CellId.MAX_LEVEL, S2CellIdSnapFunction.levelForMaxSnapRadius(S1Angle.radians(1e-30)));
  }

  /** Checks that S2CellIdSnapFunction snaps to the correct level. */
  @Test
  public void testS2CellIdSnapFunctionSnapPoint() {
    for (int iter = 0; iter < 1000; ++iter) {
      for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
        // This checks that points are snapped to the correct level, since S2CellId centers at
        // different levels are always different.
        S2CellIdSnapFunction f = new S2CellIdSnapFunction(level);
        S2Point p = data.getRandomCellId(level).toPoint();
        assertEquals(
            Platform.formatString("Iteration %d, snapping at level %d", iter, level),
            p,
            f.snapPoint(p));
      }
    }
  }

  /**
   * Tests IntLatLngSnapFunction.minSnapRadiusForExponent and exponentForMaxSnapRadius for
   * consistency over all possible exponents.
   */
  @Test
  public void testIntLatLngSnapFunctionExponentToFromSnapRadius() {
    for (int exponent = IntLatLngSnapFunction.MIN_EXPONENT;
        exponent <= IntLatLngSnapFunction.MAX_EXPONENT;
        ++exponent) {
      S1Angle radius = IntLatLngSnapFunction.minSnapRadiusForExponent(exponent);
      assertEquals(exponent, IntLatLngSnapFunction.exponentForMaxSnapRadius(radius));
      assertEquals(
          min(exponent + 1, IntLatLngSnapFunction.MAX_EXPONENT),
          IntLatLngSnapFunction.exponentForMaxSnapRadius(radius.mul(0.999)));
    }
    assertEquals(
        IntLatLngSnapFunction.MIN_EXPONENT,
        IntLatLngSnapFunction.exponentForMaxSnapRadius(S1Angle.radians(5)));
    assertEquals(
        IntLatLngSnapFunction.MAX_EXPONENT,
        IntLatLngSnapFunction.exponentForMaxSnapRadius(S1Angle.radians(1e-30)));
  }

  /**
   * Test that IntLatLngSnapFunction does not modify points that were generated using the
   * S2LatLng.from{E5,E6,E7} methods. This ensures that both functions are using bitwise-compatible
   * conversion methods.
   */
  @Test
  public void testIntLatLngSnapFunctionSnapPoint() {
    for (int iter = 0; iter < 1000; ++iter) {
      S2Point p = data.getRandomPoint();
      S2LatLng ll = new S2LatLng(p);
      S2Point p5 = S2LatLng.fromE5(ll.lat().e5(), ll.lng().e5()).toPoint();
      assertEquals(p5, new IntLatLngSnapFunction(5).snapPoint(p5));
      S2Point p6 = S2LatLng.fromE6(ll.lat().e6(), ll.lng().e6()).toPoint();
      assertEquals(p6, new IntLatLngSnapFunction(6).snapPoint(p6));
      S2Point p7 = S2LatLng.fromE7(ll.lat().e7(), ll.lng().e7()).toPoint();
      assertEquals(p7, new IntLatLngSnapFunction(7).snapPoint(p7));

      // Make sure that we're not snapping using some lower exponent.
      S2Point p7not6 = S2LatLng.fromE7(10 * ll.lat().e6() + 1, 10 * ll.lng().e6() + 1).toPoint();
      assertFalse(p7not6.equalsPoint(new IntLatLngSnapFunction(6).snapPoint(p7not6)));
    }
  }

  /**
   * The purpose of this "test" is to compute a lower bound to the fraction (minVertexSeparation() /
   * snapRadius()). Essentially this involves searching for two adjacent cells A and B such when one
   * of the corner vertices of B is snapped to the center of B, the distance to the center of A
   * decreases as much as possible. In other words, we want the ratio
   *
   * <pre>
   *   distance(center(A), center(B)) / distance(center(A), vertex(B))
   * </pre>
   *
   * to be as small as possible. We do this by considering one cell level at a time, and remembering
   * the cells that had the lowest ratios. When we proceed from one level to the next, we consider
   * all the children of those cells and keep the best ones.
   *
   * <p>The reason we can restrict the search to children of cells at the previous level is that the
   * ratio above is essentially a function of the local distortions created by projecting the S2
   * cube space onto the sphere. These distortions change smoothly over the sphere, so by keeping a
   * fairly large number of candidates ("numToKeep"), we are essentially keeping all the neighbors
   * of the optimal cell as well.
   */
  @Test
  public void testS2CellIdSnapFunctionMinVertexSeparationSnapRadiusRatio() {
    Platform.printf(System.out, "\n\ntestS2CellIdSnapFunctionMinVertexSeparationSnapRadiusRatio\n");
    double bestScore = 1e10;
    HashSet<S2CellId> bestCells = new HashSet<>();
    for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
      double score = getS2CellIdMinVertexSeparation(level, bestCells);
      bestScore = min(bestScore, score);
    }
    Platform.printf(
        System.out, "S2CellIdSnapFunction minVertexSep / snapRadius ratio: %.15f\n", bestScore);
    // Detect unexpected changes.
    assertDoubleNear(0.548490277027825, bestScore, RATIO_TOLERANCE);
  }

  /**
   * Computes the minimum edge separation (as a fraction of kMinDiag) for any snap radius at each
   * level.
   */
  @Test
  public void testS2CellIdSnapFunctionMinEdgeVertexSeparationForLevel() {
    Platform.printf(System.out, "\n\ntestS2CellIdSnapFunctionMinEdgeVertexSeparationForLevel\n");
    double score =
        getS2CellIdMinEdgeSeparation(
            "minSepForLevel",
            (level, edgeSep, minSnapRadius, maxSnapRadius) ->
                edgeSep.radians() / S2Projections.MIN_DIAG.getValue(level));
    Platform.printf(
        System.out, "S2CellIdSnapFunction minEdgeVertexSep / kMinDiag ratio: %.15f\n", score);
    // Detect unexpected changes.
    assertDoubleNear(0.397359568667803, score, RATIO_TOLERANCE);
  }

  /**
   * Computes the minimum edge separation (as a fraction of kMinDiag) for the special case where the
   * minimum snap radius is being used.
   */
  @Test
  public void testS2CellIdSnapFunctionMinEdgeVertexSeparationAtMinSnapRadius() {
    Platform.printf(
        System.out, "\n\ntestS2CellIdSnapFunctionMinEdgeVertexSeparationAtMinSnapRadius\n");
    double score =
        getS2CellIdMinEdgeSeparation(
            "minSepAtMinRadius",
            (level, edgeSep, minSnapRadius, maxSnapRadius) -> {
              double minRadiusAtLevel = S2Projections.MAX_DIAG.getValue(level) / 2;
              return (minSnapRadius.radians() <= (1 + 1e-10) * minRadiusAtLevel)
                  ? (edgeSep.radians() / S2Projections.MIN_DIAG.getValue(level))
                  : 100.0;
            });
    Platform.printf(
        System.out,
        "S2CellIdSnapFunction minEdgeVertexSep / kMinDiag at MinSnapRadiusForLevel: %.15f\n",
        score);
    // Detect unexpected changes.
    assertDoubleNear(0.565298006776224, score, RATIO_TOLERANCE);
  }

  /**
   * Computes the minimum edge separation expressed as a fraction of the maximum snap radius that
   * could yield that edge separation.
   */
  @Test
  public void testS2CellIdSnapFunctionMinEdgeVertexSeparationSnapRadiusRatio() {
    Platform.printf(
        System.out, "\n\testS2CellIdSnapFunctionMinEdgeVertexSeparationSnapRadiusRatio\n");
    double score =
        getS2CellIdMinEdgeSeparation(
            "minSepSnapRadiusRatio",
            (level, edgeSep, minSnapRadius, maxSnapRadius) ->
                edgeSep.radians() / maxSnapRadius.radians());
    Platform.printf(
        System.out, "S2CellIdSnapFunction minEdgeVertexSep / snapRadius ratio: %.15f\n", score);
    // Detect unexpected changes.
    assertDoubleNear(0.219666695288891, score, RATIO_TOLERANCE);
  }

  @Test
  public void testIntLatLngSnapFunctionMinVertexSeparationSnapRadiusRatio() {
    Platform.printf(
        System.out, "\n\ntestIntLatLngSnapFunctionMinVertexSeparationSnapRadiusRatio\n");
    double bestScore = 1e10;
    HashSet<IntLatLng> bestConfigs = new HashSet<>();
    long scale = 18;
    for (int lat0 = 0; lat0 <= 9; ++lat0) {
      bestConfigs.add(new IntLatLng(lat0, 0));
    }
    for (int exp = 0; exp <= 10; ++exp, scale *= 10) {
      double score = getLatLngMinVertexSeparation(scale, 10 * scale, bestConfigs);
      bestScore = min(bestScore, score);
    }
    Platform.printf(
        System.out, "IntLatLngSnapFunction minVertexSep / snapRadius ratio: %.15f\n", bestScore);
    // Detect unexpected changes.
    assertDoubleNear(0.471337477576603, bestScore, RATIO_TOLERANCE);
  }

  /**
   * Computes the minimum edge separation (as a fraction of kMinDiag) for any snap radius at each
   * level.
   */
  @Test
  public void testIntLatLngSnapFunctionMinEdgeVertexSeparationForLevel() {
    Platform.printf(System.out, "\n\ntestIntLatLngSnapFunctionMinEdgeVertexSeparationForLevel\n");
    double score =
        getLatLngMinEdgeSeparation(
            "minSepForLevel", (scale, edgeSep, maxSnapRadius) -> edgeSep.radians() / (PI / scale));
    Platform.printf(
        System.out, "IntLatLngSnapFunction minEdgeVertexSep / eUnit ratio: %.15f\n", score);
    // Detect unexpected changes.
    assertDoubleNear(0.277258917722462, score, RATIO_TOLERANCE);
  }

  /**
   * Computes the minimum edge separation expressed as a fraction of the maximum snap radius that
   * could yield that edge separation.
   */
  @Test
  public void testIntLatLngSnapFunctionMinEdgeVertexSeparationSnapRadiusRatio() {
    Platform.printf(
        System.out, "\n\ntestIntLatLngSnapFunctionMinEdgeVertexSeparationSnapRadiusRatio\n");
    double score =
        getLatLngMinEdgeSeparation(
            "minSepSnapRadiusRatio",
            (scale, edgeSep, maxSnapRadius) -> edgeSep.radians() / maxSnapRadius.radians());
    Platform.printf(
        System.out, "IntLatLngSnapFunction minEdgeVertexSep / snapRadius ratio: %.15f\n", score);
    // Detect unexpected changes.
    assertDoubleNear(0.222222126756717, score, RATIO_TOLERANCE);
  }

  /** Like a {@code Pair<Double, S2CellId>}. */
  private static class ScoredCellId implements Comparable<ScoredCellId> {
    final double score;
    final S2CellId id;

    public ScoredCellId(double score, S2CellId id) {
      this.score = score;
      this.id = id;
    }

    @Override
    public int compareTo(ScoredCellId that) {
      if (this.score < that.score) {
        return -1;
      } else if (this.score > that.score) {
        return 1;
      }
      return this.id.compareTo(that.id);
    }
  }

  /** Like a {@code Pair<Double, IntLatLng>}. */
  private static class ScoredIntLatLng implements Comparable<ScoredIntLatLng> {
    final double score;
    final IntLatLng latLng;

    public ScoredIntLatLng(double score, IntLatLng latLng) {
      this.score = score;
      this.latLng = latLng;
    }

    @Override
    public int compareTo(ScoredIntLatLng that) {
      if (this.score < that.score) {
        return -1;
      } else if (this.score > that.score) {
        return 1;
      }
      return this.latLng.compareTo(that.latLng);
    }
  }

  /** Like a {@code Pair<Double, LatLngConfig>}. */
  private static class ScoredLatLngConfig implements Comparable<ScoredLatLngConfig> {
    final double score;
    final LatLngConfig latLngConfig;

    public ScoredLatLngConfig(double score, LatLngConfig latLngConfig) {
      this.score = score;
      this.latLngConfig = latLngConfig;
    }

    @Override
    public int compareTo(ScoredLatLngConfig that) {
      if (this.score < that.score) {
        return -1;
      } else if (this.score > that.score) {
        return 1;
      }
      return this.latLngConfig.compareTo(that.latLngConfig);
    }
  }

  /** Helper function that computes the vertex separation between "id0" and its neighbors. */
  private static void updateS2CellIdMinVertexSeparation(
      S2CellId id0, List<ScoredCellId> scores) {
    S2Point site0 = id0.toPoint();
    List<S2CellId> neighbors = new ArrayList<>();
    id0.getAllNeighbors(id0.level(), neighbors);
    for (S2CellId id1 : neighbors) {
      S2Point site1 = id1.toPoint();
      S1Angle vertexSep = new S1Angle(site0, site1);
      S1Angle maxSnapRadius = getMaxVertexDistance(site0, id1);
      assertGreaterOrEqual(maxSnapRadius, S2CellIdSnapFunction.minSnapRadiusForLevel(id0.level()));
      double r = vertexSep.radians() / maxSnapRadius.radians();
      scores.add(new ScoredCellId(r, id0));
    }
  }

  /** Helper for testS2CellIdSnapFunctionMinVertexSeparationSnapRadiusRatio. */
  private static double getS2CellIdMinVertexSeparation(int level, HashSet<S2CellId> bestCells) {
    // The worst-case separation ratios always occur when the snapRadius is not much larger than the
    // minimum, since this allows the site spacing to be reduced by as large a fraction as possible.
    //
    // For the minimum vertex separation ratio, we choose a site and one of its 8-way neighbors,
    // then look at the ratio of the distance to the center of that neighbor to the distance to the
    // furthest corner of that neighbor (which is the largest possible snap radius for this
    // configuration).
    List<ScoredCellId> scores = new ArrayList<>();
    if (level == 0) {
      updateS2CellIdMinVertexSeparation(SEARCH_ROOT_ID, scores);
    } else {
      for (S2CellId parent : bestCells) {
        for (S2CellId id0 = parent.childBegin(); !id0.equals(parent.childEnd()); id0 = id0.next()) {
          updateS2CellIdMinVertexSeparation(id0, scores);
        }
      }
    }

    // Now sort and deduplicate the entries, print out the "numToPrint" best ones, and keep the best
    // "numToKeep" of them to seed the next round.
    scores.sort(cellDistanceComparator);
    S2BuilderUtil.deduplicateSortedList(scores);
    bestCells.clear();
    int numToKeep = 300;
    int numToPrint = 1;
    for (ScoredCellId entry : scores) {
      S2CellId id = entry.id;
      if (--numToPrint >= 0) {
        R2Vector uv = id.getCenterUV();
        Platform.printf(
            System.out,
            "S2CellId Level %2d: minVertexSepRatio = %.15f u=%.6f v=%.6f %s\n",
            level, entry.score, uv.get(0), uv.get(1), id.toToken());
      }
      if (SEARCH_FOCUS_ID.contains(id) || id.contains(SEARCH_FOCUS_ID)) {
        if (bestCells.add(id) && --numToKeep <= 0) {
          break;
        }
      }
    }

    return scores.get(0).score;
  }

  /** Finds the largest distance from the point 'p' to any of the four vertices of 'id'. */
  private static S1Angle getMaxVertexDistance(S2Point p, S2CellId id) {
    S2Cell cell = new S2Cell(id);
    return S1Angle.max(
        S1Angle.max(new S1Angle(p, cell.getVertex(0)), new S1Angle(p, cell.getVertex(1))),
        S1Angle.max(new S1Angle(p, cell.getVertex(2)), new S1Angle(p, cell.getVertex(3))));
  }

  // TODO(torrey): getCircumradius seems to be a generally useful method. Perhaps define "too far"
  // etc. and move it to a production library for both Java and C++.

  /**
   * Returns the radius of a circle that goes through the three points, or PI if they are too far
   * apart.
   */
  private static S1Angle getCircumradius(S2Point a, S2Point b, S2Point c) {
    // We return this value if the circumradius is very large.
    S1Angle kTooBig = S1Angle.radians(PI);
    double turnAngle = S2.turnAngle(a, b, c);
    if (abs(Platform.IEEEremainder(turnAngle, PI)) < 1e-2) {
      return kTooBig;
    }

    double a2 = b.sub(c).norm2();
    double b2 = c.sub(a).norm2();
    double c2 = a.sub(b).norm2();
    if (a2 > 2 || b2 > 2 || c2 > 2) {
      return kTooBig;
    }
    double ma = a2 * (b2 + c2 - a2);
    double mb = b2 * (c2 + a2 - b2);
    double mc = c2 * (a2 + b2 - c2);
    S2Point p = S2Point.sum(a.mul(ma), b.mul(mb), c.mul(mc)).div(ma + mb + mc);
    return new S1Angle(p, a);
  }

  /**
   * Gets the neighbors of the given S2CellId 'id', and the neighbors of those neighbors (excluding
   * 'id' itself.
   */
  private static Set<S2CellId> getNeighbors(S2CellId id) {
    int kNumLayers = 2;
    Set<S2CellId> neighbors = new HashSet<>();
    neighbors.add(id);
    for (int layer = 0; layer < kNumLayers; ++layer) {
      List<S2CellId> newNeighbors = new ArrayList<>();
      for (S2CellId nbr : neighbors) {
        nbr.getAllNeighbors(id.level(), newNeighbors);
      }
      neighbors.addAll(newNeighbors);
      neighbors.remove(id);
    }
    return neighbors;
  }

  /**
   * S2CellIdMinEdgeSeparationFunction defines an objective function that will be optimized by
   * {@link #getS2CellIdMinEdgeSeparation(String, S2CellIdMinEdgeSeparationFunction, int, HashSet)}
   * by finding worst-case configurations of S2CellIds. We use this to find the worst cases under
   * various conditions (e.g., when the minimum snap radius at a given level is being used). The
   * objective function is called for a specific configuration of vertices that are snapped at the
   * given S2CellId level. "edgeSep" is the edge-vertex distance that is achieved by this
   * configuration, and "minSnapRadius" and "maxSnapRadius" are the minimum and maximum snap radii
   * for which this configuration is valid (i.e., where the desired snapping will take place).
   */
  private interface S2CellIdMinEdgeSeparationFunction {
    double apply(int level, S1Angle edgeSep, S1Angle minSnapRadius, S1Angle maxSnapRadius);
  }

  /**
   * Returns the minimum value of the given objective function over sets of nearby vertices that are
   * designed to minimize the edge-vertex separation when an edge is snapped.
   */
  private static double getS2CellIdMinEdgeSeparation(
      String label,
      S2CellIdMinEdgeSeparationFunction objective,
      int level,
      HashSet<S2CellId> bestCells) {
    // To find minimum edge separations, we choose a cell ("id0") and two nearby cells ("id1" and
    // "id2"), where "nearby" is defined by getNeighbors(). Let "site0", "site1", and "site2" be the
    // centers of these cells. The idea is to consider an input edge E that intersects the Voronoi
    // regions of "site1" and "site2" (and therefore snaps to an edge E' between these sites) but
    // does not intersect the Voronoi region of "site0" (and therefore can't be snapped to site0).
    // The goal is to search for snapped edges E' that approach site0 as closely as possible.
    //
    // To do this, we first compute the circumradius of the three cell centers ("site0", "site1",
    // and "site2"); this is the minimum snap radius in order for it to be possible to construct an
    // edge E that snaps to "site1" and "site2" but not to "site0". We also compute the distance
    // from "site0" to the snapped edge. Next we find the corner vertex of "id1" and "id2" that is
    // furthest from "site0"; the smaller of these two distances is the maximum snap radius such
    // that "site1" and "site2" can be chosen as sites after choosing "site0". If the maximum is
    // less than the minimum, then this configuration is rejected; otherwise we evaluate the given
    // objective function and keep the configurations that result in the smallest values.
    //
    // The optimization process works by keeping track of the set of S2CellIds that yielded the
    // best results at the previous level, and exploring all the nearby neighbor combinations of the
    // children of those cells at the next level. In order to get better coverage, we keep track of
    // the best score and configuration (i.e. the two neighboring cells "id1" and "id2") for each
    // initial cell "id0".
    HashMap<S2CellId, Double> bestScores = new HashMap<>();
    HashMap<S2CellId, Pair<S2CellId, S2CellId>> bestConfigs = new HashMap<>();
    for (S2CellId parent : bestCells) {
      for (S2CellId id0 = parent.childBegin(level);
          !id0.equals(parent.childEnd(level));
          id0 = id0.next()) {
        S2Point site0 = id0.toPoint();
        Set<S2CellId> neighbors = getNeighbors(id0);
        for (S2CellId id1 : neighbors) {
          S2Point site1 = id1.toPoint();
          S1Angle maxV1 = getMaxVertexDistance(site0, id1);
          for (S2CellId id2 : neighbors) {
            if (id2.lessOrEquals(id1)) {
              continue;
            }
            S2Point site2 = id2.toPoint();
            S1Angle minSnapRadius = getCircumradius(site0, site1, site2);
            if (minSnapRadius.greaterThan(S2BuilderSnapFunctions.maxSnapRadius())) {
              continue;
            }
            // Note that it is only the original points *before* snapping that need to be at least
            // "snapRadius" away from "site0". The points after snapping ("site1" and "site2") may
            // be closer.
            S1Angle maxV2 = getMaxVertexDistance(site0, id2);
            S1Angle maxSnapRadius = S1Angle.min(maxV1, maxV2);
            if (minSnapRadius.greaterThan(maxSnapRadius)) {
              continue;
            }
            assertGreaterOrEqual(maxSnapRadius, S2CellIdSnapFunction.minSnapRadiusForLevel(level));

            // This is a valid configuration, so evaluate it.
            S1Angle edgeSep = S2EdgeUtil.getDistance(site0, site1, site2);
            double score = objective.apply(level, edgeSep, minSnapRadius, maxSnapRadius);
            double bestScore = bestScores.getOrDefault(id0, 0.0);
            if (bestScore == 0 || bestScore > score) {
              bestScores.put(id0, score);
              bestConfigs.put(id0, Pair.of(id1, id2));
            }
          }
        }
      }
    }

    // Now sort the entries, print out the "numToPrint" best ones, and generate a set of candidates
    // for the next round by generating all the 8-way neighbors of the best candidates, and keeping
    // up to "numToKeep" of them. The results vary slightly according to how many candidates we
    // keep, but the variations are much smaller than the conservative assumptions made by the
    // S2CellIdSnapFunction implementation.
    int numToKeep = 100;
    int numToPrint = 3;
    List<ScoredCellId> sorted = new ArrayList<>();
    for (Map.Entry<S2CellId, Double> entry : bestScores.entrySet()) {
      sorted.add(new ScoredCellId(entry.getValue(), entry.getKey()));
    }
    sorted.sort(cellDistanceComparator);
    bestCells.clear();
    Platform.printf(System.out, "S2CellIdMinEdgeSeparation Level %d:\n", level);
    for (ScoredCellId entry : sorted) {
      S2CellId id = entry.id;
      if (--numToPrint >= 0) {
        R2Vector uv = id.getCenterUV();
        Pair<S2CellId, S2CellId> neighbors = bestConfigs.get(id);
        Platform.printf(
            System.out,
            "  %s = %.15f u=%7.4f v=%7.4f %s %s %s\n",
            label,
            entry.score,
            uv.get(0),
            uv.get(1),
            id.toToken(),
            neighbors.first.toToken(),
            neighbors.second.toToken());
      }
      List<S2CellId> neighbors = new ArrayList<>();
      neighbors.add(id);
      id.getAllNeighbors(id.level(), neighbors);
      for (S2CellId nbr : neighbors) {
        // The S2Cell hierarchy has many regions that are symmetrical. We can eliminate most of the
        // "duplicates" by restricting the search to cells in kS2CellIdFocus.
        if (SEARCH_FOCUS_ID.contains(nbr) || nbr.contains(SEARCH_FOCUS_ID)) {
          if (bestCells.add(nbr) && --numToKeep <= 0) {
            return sorted.get(0).score;
          }
        }
      }
    }
    return sorted.get(0).score;
  }

  /** Helper for testS2CellIdSnapFunctionMinEdgeVertexSeparationForLevel. */
  private static double getS2CellIdMinEdgeSeparation(
      String label, S2CellIdMinEdgeSeparationFunction objective) {
    double bestScore = 1e10;
    HashSet<S2CellId> bestCells = new HashSet<>();
    bestCells.add(SEARCH_ROOT_ID);
    for (int level = 0; level <= S2CellId.MAX_LEVEL; ++level) {
      double score = getS2CellIdMinEdgeSeparation(label, objective, level, bestCells);
      bestScore = min(bestScore, score);
    }
    return bestScore;
  }

  private static double getLatLngMinVertexSeparation(
      long oldScale, long scale, HashSet<IntLatLng> bestConfigs) {
    // The worst-case separation ratios always occur when the snapRadius is not much larger than
    // the minimum, since this allows the site spacing to be reduced by as large a fraction as
    // possible.
    //
    // For the minimum vertex separation ratio, we choose a site and one of its 8-way neighbors,
    // then look at the ratio of the distance to the center of that neighbor to the distance to the
    // furthest corner of that neighbor (which is the largest possible snap radius for this
    // configuration).
    S1Angle minSnapRadiusAtScale = S1Angle.radians(M_SQRT1_2 * PI / scale);
    List<ScoredIntLatLng> scores = new ArrayList<>();
    double scaleFactor = (double) scale / oldScale;
    for (IntLatLng parent : bestConfigs) {
      IntLatLng newParent = parent.rescale(scaleFactor);
      for (int dlat0 = -7; dlat0 <= 7; ++dlat0) {
        IntLatLng ll0 = newParent.add(new IntLatLng(dlat0, 0));
        if (!ll0.isValid(scale) || ll0.lat < 0) {
          continue;
        }
        S2Point site0 = ll0.toPoint(scale);
        for (int dlat1 = 0; dlat1 <= 2; ++dlat1) {
          for (int dlng1 = 0; dlng1 <= 5; ++dlng1) {
            IntLatLng ll1 = IntLatLng.add(ll0, new IntLatLng(dlat1, dlng1));
            if (ll1.isEqualTo(ll0) || !ll1.hasValidVertices(scale)) {
              continue;
            }
            S1Angle maxSnapRadius = ll1.getMaxVertexDistance(site0, scale);
            if (maxSnapRadius.lessThan(minSnapRadiusAtScale)) {
              continue;
            }
            S2Point site1 = ll1.toPoint(scale);
            S1Angle vertexSep = new S1Angle(site0, site1);
            double r = vertexSep.radians() / maxSnapRadius.radians();
            scores.add(new ScoredIntLatLng(r, ll0));
          }
        }
      }
    }

    // Now sort and deduplicate the entries, print out the "numToPrint" best ones, and keep the best
    // "numToKeep" of them to seed the next round.
    scores.sort(intLatLngDistanceComparator);
    S2BuilderUtil.deduplicateSortedList(scores);

    bestConfigs.clear();
    int numToKeep = 100;
    int numToPrint = 1;
    for (ScoredIntLatLng entry : scores) {
      if (--numToPrint >= 0) {
        Platform.printf(
            System.out,
            "LatLngMinVertexSeparation Scale %14d: minVertexSepRatio = %.15f, %s\n",
            scale, entry.score, S2TextFormat.toString(entry.latLng.toPoint(scale)));
      }
      if (bestConfigs.add(entry.latLng) && --numToKeep <= 0) {
        break;
      }
    }
    return scores.get(0).score;
  }

  private interface LatLngMinEdgeSeparationFunction {
    double apply(long scale, S1Angle edgeSep, S1Angle maxSnapRadius);
  }

  private static double getLatLngMinEdgeSeparation(
      String label,
      LatLngMinEdgeSeparationFunction objective,
      long scale,
      List<LatLngConfig> bestConfigs) {
    S1Angle minSnapRadiusAtScale = S1Angle.radians(M_SQRT1_2 * PI / scale);
    List<ScoredLatLngConfig> scores = new ArrayList<>();
    for (LatLngConfig oldScaleParent : bestConfigs) {
      // To reduce duplicates, we require that site0 always has longitude 0.
      Preconditions.checkState(oldScaleParent.ll0.lng == 0);
      // Get oldScaleParent scaled to the specified scale.
      LatLngConfig parent = oldScaleParent.rescale(scale);

      for (int dlat0 = -1; dlat0 <= 1; ++dlat0) {
        IntLatLng ll0 = IntLatLng.add(parent.ll0, new IntLatLng(dlat0, 0));
        // To reduce duplicates, we require that site0.latitude >= 0.
        if (!ll0.isValid(scale) || ll0.lat < 0) {
          continue;
        }
        S2Point site0 = ll0.toPoint(scale);
        for (int dlat1 = -1; dlat1 <= 1; ++dlat1) {
          for (int dlng1 = -2; dlng1 <= 2; ++dlng1) {
            IntLatLng ll1 = IntLatLng.add(parent.ll1, new IntLatLng(dlat0 + dlat1, dlng1));
            if (ll1.isEqualTo(ll0) || !ll1.hasValidVertices(scale)) {
              continue;
            }
            // Only consider neighbors within 2 latitude units of site0.
            if (abs(ll1.lat - ll0.lat) > 2) {
              continue;
            }

            S2Point site1 = ll1.toPoint(scale);
            S1Angle maxV1 = ll1.getMaxVertexDistance(site0, scale);
            for (int dlat2 = -1; dlat2 <= 1; ++dlat2) {
              for (int dlng2 = -2; dlng2 <= 2; ++dlng2) {
                IntLatLng ll2 = IntLatLng.add(parent.ll2, new IntLatLng(dlat0 + dlat2, dlng2));
                if (!ll2.hasValidVertices(scale)) {
                  continue;
                }
                // Only consider neighbors within 2 latitude units of site0.
                if (abs(ll2.lat - ll0.lat) > 2) {
                  continue;
                }
                // To reduce duplicates, we require ll1 < ll2 lexicographically and
                // site2.longitude >= 0. (It's *not* okay to require site1.longitude >= 0, because
                // then some configurations with site1.latitude == site2.latitude would be missed.)
                if (ll2.isLessThanOrEqual(ll1) || ll2.lng < 0) {
                  continue;
                }

                S2Point site2 = ll2.toPoint(scale);
                S1Angle minSnapRadius = getCircumradius(site0, site1, site2);
                if (minSnapRadius.greaterThan(S2BuilderSnapFunctions.maxSnapRadius())) {
                  continue;
                }
                // Only the original points *before* snapping that need to be at least "snapRadius"
                // away from "site0". The points after snapping ("site1" and "site2") may be closer.
                S1Angle maxV2 = ll2.getMaxVertexDistance(site0, scale);
                S1Angle maxSnapRadius = S1Angle.min(maxV1, maxV2);
                if (minSnapRadius.greaterThan(maxSnapRadius)) {
                  continue;
                }
                if (maxSnapRadius.lessThan(minSnapRadiusAtScale)) {
                  continue;
                }

                // This is a valid configuration, so evaluate it.
                S1Angle edgeSep = S2EdgeUtil.getDistance(site0, site1, site2);
                double score = objective.apply(scale, edgeSep, maxSnapRadius);
                LatLngConfig config = new LatLngConfig(scale, ll0, ll1, ll2);
                scores.add(new ScoredLatLngConfig(score, config));
              }
            }
          }
        }
      }
    }
    Preconditions.checkState(!scores.isEmpty());

    // Now sort and deduplicate the entries, print out the "numToPrint" best ones, and keep the best
    // "numToKeep" of them to seed the next round.
    scores.sort(latLngConfigDistanceComparator);
    S2BuilderUtil.deduplicateSortedList(scores);

    bestConfigs.clear();
    int numToKeep = 200;
    int numToPrint = 3;
    Platform.printf(System.out, "LatLngMinEdgeSeparation Scale %d:\n", scale);
    for (ScoredLatLngConfig entry : scores) {
      LatLngConfig config = entry.latLngConfig;
      if (--numToPrint >= 0) {
        Platform.printf(
            System.out,
            "  %s = %.15f %s %s %s\n",
            label,
            entry.score,
            S2TextFormat.toString(config.ll0.toPoint(config.scale)),
            S2TextFormat.toString(config.ll1.toPoint(config.scale)),
            S2TextFormat.toString(config.ll2.toPoint(config.scale)));
      }
      // Optional: filter the candidates to concentrate on a specific region (e.g., the north pole).
      bestConfigs.add(config);
      if (--numToKeep <= 0) {
        break;
      }
    }
    return scores.get(0).score;
  }

  private static double getLatLngMinEdgeSeparation(
      String label, LatLngMinEdgeSeparationFunction objective) {
    double bestScore = 1e10;
    List<LatLngConfig> bestConfigs = new ArrayList<>();
    long scale = 6; // Initially points are 30 degrees apart.
    long maxLng = scale;
    long maxLat = scale / 2;
    for (int lat0 = 0; lat0 <= maxLat; ++lat0) {
      for (int lat1 = lat0 - 2; lat1 <= min(maxLat, lat0 + 2); ++lat1) {
        for (int lng1 = 0; lng1 <= maxLng; ++lng1) {
          for (int lat2 = lat1; lat2 <= min(maxLat, lat0 + 2); ++lat2) {
            for (int lng2 = 0; lng2 <= maxLng; ++lng2) {
              IntLatLng ll0 = new IntLatLng(lat0, 0);
              IntLatLng ll1 = new IntLatLng(lat1, lng1);
              IntLatLng ll2 = new IntLatLng(lat2, lng2);
              if (ll2.isLessThanOrEqual(ll1)) {
                continue;
              }
              bestConfigs.add(new LatLngConfig(scale, ll0, ll1, ll2));
            }
          }
        }
      }
    }
    long targetScale = 180;
    for (int exp = 0; exp <= 10; ++exp, targetScale *= 10) {
      while (scale < targetScale) {
        scale = min((long) (1.8 * scale), targetScale);
        double score = getLatLngMinEdgeSeparation(label, objective, scale, bestConfigs);
        if (scale == targetScale) {
          bestScore = min(bestScore, score);
        }
      }
    }
    return bestScore;
  }

  /**
   * Comparator for ScoredCellIds, comparing by score (distance), breaking ties by comparing cell
   * ids.
   */
  private static final Comparator<ScoredCellId> cellDistanceComparator =
      (left, right) -> {
        if (left.score == right.score) {
          // Deterministic tie breaking.
          return left.id.compareTo(right.id);
        }
        return Double.compare(left.score, right.score);
      };

  /**
   * Comparator for ScoredIntLatLng, comparing by score (distance), breaking ties by comparing
   * IntLatLngs. Only makes sense if the IntLatLngs are using the same scale.
   */
  private static final Comparator<ScoredIntLatLng> intLatLngDistanceComparator =
      (left, right) -> {
        if (left.score == right.score) {
          // Deterministic tie breaking.
          return left.latLng.compareTo(right.latLng);
        }
        return Double.compare(left.score, right.score);
      };

  /**
   * Comparator for ScoredLatLngConfig. Compares by score (distance), breaking ties by comparing
   * LatLngConfigs. Only makes sense if they have the same scale.
   */
  private static final Comparator<ScoredLatLngConfig> latLngConfigDistanceComparator =
      (left, right) -> {
        if (left.score == right.score) {
          // Deterministic tie breaking.
          return left.latLngConfig.compareTo(right.latLngConfig);
        }
        return Double.compare(left.score, right.score);
      };

  /**
   * A scaled S2LatLng with integer coordinates, similar to E7 coordinates, except that the scale is
   * variable. The coordinates are multiplied by (PI / scale) to convert them to radians. For
   * example, IntLatLngs using a scale = 180,000 have coordinates in millidegrees.
   *
   * <p>See also {@link LatLngConfig} below.
   */
  @SuppressWarnings("AmbiguousMethodReference")
  private static class IntLatLng implements Comparable<IntLatLng> {
    public final long lat;
    public final long lng;

    /** Constructs an IntLatLng with the given values. */
    public IntLatLng(long lat, long lng) {
      this.lat = lat;
      this.lng = lng;
    }

    /** Returns true if this IntLatLng is equal to 'other'. */
    public boolean isEqualTo(IntLatLng other) {
      return this.lat == other.lat && this.lng == other.lng;
    }

    @Override
    public int compareTo(IntLatLng other) {
      if (this.lat == other.lat) {
        return Long.compare(this.lng, other.lng);
      }
      return Long.compare(this.lat, other.lat);
    }

    /** Returns true if this IntLatLng compares as less than or equal to 'other'. */
    public boolean isLessThanOrEqual(IntLatLng other) {
      return compareTo(other) <= 0;
    }

    /** Returns a new IntLatLng with 'lat' and 'lng' equal the sum of 'a' and 'b'. */
    public static IntLatLng add(IntLatLng a, IntLatLng b) {
      return new IntLatLng(a.lat + b.lat, a.lng + b.lng);
    }

    /** Returns a new IntLatLng with 'lat' and 'lng' equal the sum of 'other' and 'this'. */
    public IntLatLng add(IntLatLng other) {
      return add(this, other);
    }

    /** Returns a new IntLatLng with 'lat' and 'lng' equal to those in 'll' multiplied by 'm'. */
    public static IntLatLng mul(IntLatLng ll, int m) {
      return new IntLatLng(ll.lat * m, ll.lng * m);
    }

    /** Returns a new IntLatLng with 'lat' and 'lng' equal to this, multiplied by 'm'. */
    public IntLatLng mul(int m) {
      return mul(this, m);
    }

    /**
     * Returns true if the lat and lng values are within the bounds of the given 'scale'. A
     * coordinate value of "scale" corresponds to 180 degrees.
     */
    public boolean isValid(long scale) {
      return (abs(lat) <= scale / 2 && abs(lng) <= scale);
    }

    /**
     * Like isValid, but excludes latitudes of 90 and longitudes of 180. A coordinate value of
     * "scale" corresponds to 180 degrees.
     */
    public boolean hasValidVertices(long scale) {
      return (abs(lat) < scale / 2 && abs(lng) < scale);
    }

    /** Returns a new IntLatLng scaled to the given scaleFactor. */
    public IntLatLng rescale(double scaleFactor) {
      return new IntLatLng(round(scaleFactor * lat), round(scaleFactor * lng));
    }

    /** Converts this IntLatLng to an S2Point using the provided 'scale'. */
    public S2Point toPoint(long scale) {
      return S2LatLng.fromRadians(lat * (PI / scale), lng * (PI / scale)).toPoint();
    }

    /**
     * Treats this IntLatLng like a square with edge length (PI / scale) radians, centered on the
     * (lat,lng) coordinates. For 'i' in [0..3], returns S2Points at the four corners of the square,
     * in CCW order starting from the lower left. Consecutive corners are separated by one scale
     * unit, with this IntLatLng at the center, half a scale unit from each edge of the square.
     * Adjacent IntLatLngs share two corners.
     */
    public S2Point getVertex(long scale, int i) {
      int dlat = (i == 0 || i == 3) ? -1 : 1;
      int dlng = (i == 0 || i == 1) ? -1 : 1;
      return IntLatLng.add(mul(2), new IntLatLng(dlat, dlng)).toPoint(2 * scale);
    }

    /**
     * Computes the largest distance on the surface of the sphere between this IntLatLng and the
     * four 'corners' around it defined by {@link #getVertex(long, int)}. Due to distortions of
     * the S2 projection, those spherical distances are not equal.
     */
    public S1Angle getMaxVertexDistance(S2Point p, long scale) {
      return S1Angle.max(
          S1Angle.max(
              new S1Angle(p, getVertex(scale, 0)), new S1Angle(p, getVertex(scale, 1))),
          S1Angle.max(
              new S1Angle(p, getVertex(scale, 2)), new S1Angle(p, getVertex(scale, 3))));
    }
  }

  /**
   * A triple of scaled S2LatLng coordinates. The coordinates are multiplied by (PI / scale) to
   * convert them to radians.
   */
  private static class LatLngConfig implements Comparable<LatLngConfig> {
    final long scale;
    final IntLatLng ll0;
    final IntLatLng ll1;
    final IntLatLng ll2;

    public LatLngConfig(long scale, IntLatLng ll0, IntLatLng ll1, IntLatLng ll2) {
      this.scale = scale;
      this.ll0 = ll0;
      this.ll1 = ll1;
      this.ll2 = ll2;
    }

    /** Returns a new LatLngConfig based on the current one, but scaled to the given 'scale'. */
    public LatLngConfig rescale(long newScale) {
      double scaleFactor = (double) newScale / scale;
      return new LatLngConfig(
          scale,
          ll0.rescale(scaleFactor),
          ll1.rescale(scaleFactor),
          ll2.rescale(scaleFactor));
    }

    @Override
    public int compareTo(LatLngConfig other) {
      Preconditions.checkState(scale == other.scale);
      int c = this.ll0.compareTo(other.ll0);
      if (c != 0) {
        return c;
      }
      c = this.ll1.compareTo(other.ll1);
      if (c != 0) {
        return c;
      }
      return this.ll2.compareTo(other.ll2);
    }
  }
}
