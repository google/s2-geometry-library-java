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

import static com.google.common.geometry.S2Projections.AVG_ANGLE_SPAN;
import static com.google.common.geometry.S2Projections.AVG_AREA;
import static com.google.common.geometry.S2Projections.AVG_DIAG;
import static com.google.common.geometry.S2Projections.AVG_EDGE;
import static com.google.common.geometry.S2Projections.AVG_WIDTH;
import static com.google.common.geometry.S2Projections.MAX_ANGLE_SPAN;
import static com.google.common.geometry.S2Projections.MAX_AREA;
import static com.google.common.geometry.S2Projections.MAX_DIAG;
import static com.google.common.geometry.S2Projections.MAX_DIAG_ASPECT;
import static com.google.common.geometry.S2Projections.MAX_EDGE;
import static com.google.common.geometry.S2Projections.MAX_EDGE_ASPECT;
import static com.google.common.geometry.S2Projections.MAX_WIDTH;
import static com.google.common.geometry.S2Projections.MIN_ANGLE_SPAN;
import static com.google.common.geometry.S2Projections.MIN_AREA;
import static com.google.common.geometry.S2Projections.MIN_DIAG;
import static com.google.common.geometry.S2Projections.MIN_EDGE;
import static com.google.common.geometry.S2Projections.MIN_WIDTH;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.log;
import static java.lang.Math.log1p;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.annotations.GwtIncompatible;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for S2Cell. */
@RunWith(JUnit4.class)
public class S2CellTest extends GeometryTestCase {

  public static final boolean DEBUG_MODE = true;

  @Test
  public void testFaces() {
    Map<S2Point, Integer> edgeCounts = new HashMap<>();
    Map<S2Point, Integer> vertexCounts = new HashMap<>();
    for (int face = 0; face < 6; ++face) {
      S2CellId id = S2CellId.fromFacePosLevel(face, 0, 0);
      S2Cell cell = new S2Cell(id);
      assertEquals(cell.id(), id);
      assertEquals(cell.face(), face);
      assertEquals(0, cell.level());
      // Top-level faces have alternating orientations to get RHS coordinates.
      assertEquals(cell.orientation(), face & S2.SWAP_MASK);
      assertTrue(!cell.isLeaf());
      for (int k = 0; k < 4; ++k) {
        if (edgeCounts.containsKey(cell.getEdgeRaw(k))) {
          edgeCounts.put(cell.getEdgeRaw(k), edgeCounts.get(cell.getEdgeRaw(k)) + 1);
        } else {
          edgeCounts.put(cell.getEdgeRaw(k), 1);
        }

        if (vertexCounts.containsKey(cell.getVertexRaw(k))) {
          vertexCounts.put(cell.getVertexRaw(k), vertexCounts.get(cell.getVertexRaw(k)) + 1);
        } else {
          vertexCounts.put(cell.getVertexRaw(k), 1);
        }
        assertAlmostEquals(cell.getVertexRaw(k).dotProd(cell.getEdgeRaw(k)), 0);
        assertAlmostEquals(cell.getVertexRaw(k + 1).dotProd(cell.getEdgeRaw(k)), 0);
        assertAlmostEquals(
            cell.getVertexRaw(k).crossProd(cell.getVertexRaw(k + 1)).normalize()
                .dotProd(cell.getEdge(k)),
            1.0);
      }
    }
    // Check that edges have multiplicity 2 and vertices have multiplicity 3.
    for (Integer i : edgeCounts.values()) {
      assertEquals(2, i.intValue());
    }
    for (Integer i : vertexCounts.values()) {
      assertEquals(3, i.intValue());
    }
  }

  static class LevelStats {
    double count;
    double minArea;
    double maxArea;
    double avgArea;
    double minWidth;
    double maxWidth;
    double avgWidth;
    double minEdge;
    double maxEdge;
    double avgEdge;
    double maxEdgeAspect;
    double minDiag;
    double maxDiag;
    double avgDiag;
    double maxDiagAspect;
    double minAngleSpan;
    double maxAngleSpan;
    double avgAngleSpan;
    double minApproxRatio;
    double maxApproxRatio;

    LevelStats() {
      count = 0;
      minArea = 100;
      maxArea = 0;
      avgArea = 0;
      minWidth = 100;
      maxWidth = 0;
      avgWidth = 0;
      minEdge = 100;
      maxEdge = 0;
      avgEdge = 0;
      maxEdgeAspect = 0;
      minDiag = 100;
      maxDiag = 0;
      avgDiag = 0;
      maxDiagAspect = 0;
      minAngleSpan = 100;
      maxAngleSpan = 0;
      avgAngleSpan = 0;
      minApproxRatio = 100;
      maxApproxRatio = 0;
    }
  }

  static List<LevelStats> levelStats = new ArrayList<>(S2CellId.MAX_LEVEL + 1);

  static {
    for (int i = 0; i < S2CellId.MAX_LEVEL + 1; ++i) {
      levelStats.add(new LevelStats());
    }
  }

  static void gatherStats(S2Cell cell) {
    LevelStats s = levelStats.get(cell.level());
    double exactArea = cell.exactArea();
    double approxArea = cell.approxArea();
    double minEdge = 100;
    double maxEdge = 0;
    double avgEdge = 0;
    double minDiag = 100;
    double maxDiag = 0;
    double minWidth = 100;
    double maxWidth = 0;
    double minAngleSpan = 100;
    double maxAngleSpan = 0;
    for (int i = 0; i < 4; ++i) {
      double edge = cell.getVertexRaw(i).angle(cell.getVertexRaw(i + 1));
      minEdge = min(edge, minEdge);
      maxEdge = max(edge, maxEdge);
      avgEdge += 0.25 * edge;
      S2Point mid = cell.getVertexRaw(i).add(cell.getVertexRaw(i + 1));
      double width = S2.M_PI_2 - mid.angle(cell.getEdgeRaw(i ^ 2));
      minWidth = min(width, minWidth);
      maxWidth = max(width, maxWidth);
      if (i < 2) {
        double diag = cell.getVertexRaw(i).angle(cell.getVertexRaw(i ^ 2));
        minDiag = min(diag, minDiag);
        maxDiag = max(diag, maxDiag);
        double angleSpan = cell.getEdgeRaw(i).angle(cell.getEdgeRaw(i ^ 2).neg());
        minAngleSpan = min(angleSpan, minAngleSpan);
        maxAngleSpan = max(angleSpan, maxAngleSpan);
      }
    }
    s.count += 1;
    s.minArea = min(exactArea, s.minArea);
    s.maxArea = max(exactArea, s.maxArea);
    s.avgArea += exactArea;
    s.minWidth = min(minWidth, s.minWidth);
    s.maxWidth = max(maxWidth, s.maxWidth);
    s.avgWidth += 0.5 * (minWidth + maxWidth);
    s.minEdge = min(minEdge, s.minEdge);
    s.maxEdge = max(maxEdge, s.maxEdge);
    s.avgEdge += avgEdge;
    s.maxEdgeAspect = max(maxEdge / minEdge, s.maxEdgeAspect);
    s.minDiag = min(minDiag, s.minDiag);
    s.maxDiag = max(maxDiag, s.maxDiag);
    s.avgDiag += 0.5 * (minDiag + maxDiag);
    s.maxDiagAspect = max(maxDiag / minDiag, s.maxDiagAspect);
    s.minAngleSpan = min(minAngleSpan, s.minAngleSpan);
    s.maxAngleSpan = max(maxAngleSpan, s.maxAngleSpan);
    s.avgAngleSpan += 0.5 * (minAngleSpan + maxAngleSpan);
    double approxRatio = approxArea / exactArea;
    s.minApproxRatio = min(approxRatio, s.minApproxRatio);
    s.maxApproxRatio = max(approxRatio, s.maxApproxRatio);
  }

  public void checkSubdivide(S2Cell cell) {
    gatherStats(cell);
    if (cell.isLeaf()) {
      return;
    }

    S2Cell[] children = new S2Cell[4];
    for (int i = 0; i < children.length; ++i) {
      children[i] = new S2Cell();
    }
    assertTrue(cell.subdivide(children));
    S2CellId childId = cell.id().childBegin();
    double exactArea = 0;
    double approxArea = 0;
    double averageArea = 0;
    for (int i = 0; i < 4; ++i, childId = childId.next()) {
      exactArea += children[i].exactArea();
      approxArea += children[i].approxArea();
      averageArea += children[i].averageArea();

      // Check that the child geometry is consistent with its cell id.
      assertEquals(children[i].id(), childId);
      assertTrue(children[i].getCenter().aequal(childId.toPoint(), 1e-15));
      S2Cell direct = new S2Cell(childId);
      assertEquals(children[i].face(), direct.face());
      assertEquals(children[i].level(), direct.level());
      assertEquals(children[i].orientation(), direct.orientation());
      assertEquals(children[i].getCenterRaw(), direct.getCenterRaw());
      for (int k = 0; k < 4; ++k) {
        assertEquals(children[i].getVertexRaw(k), direct.getVertexRaw(k));
        assertEquals(children[i].getEdgeRaw(k), direct.getEdgeRaw(k));
      }

      // Test contains() and mayIntersect().
      assertTrue(cell.contains(children[i]));
      assertTrue(cell.mayIntersect(children[i]));
      assertFalse(children[i].contains(cell));
      assertTrue(cell.contains(children[i].getCenterRaw()));
      for (int j = 0; j < 4; ++j) {
        assertTrue(cell.contains(children[i].getVertexRaw(j)));
        if (j != i) {
          assertFalse(children[i].contains(children[j].getCenterRaw()));
          assertFalse(children[i].mayIntersect(children[j]));
        }
      }

      // Test getCapBound and getRectBound.
      S2Cap parentCap = cell.getCapBound();
      S2LatLngRect parentRect = cell.getRectBound();
      if (cell.contains(new S2Point(0, 0, 1)) || cell.contains(new S2Point(0, 0, -1))) {
        assertTrue(parentRect.lng().isFull());
      }
      S2Cap childCap = children[i].getCapBound();
      S2LatLngRect childRect = children[i].getRectBound();
      assertTrue(childCap.contains(children[i].getCenter()));
      assertTrue(childRect.contains(children[i].getCenterRaw()));
      assertTrue(parentCap.contains(children[i].getCenter()));
      assertTrue(parentRect.contains(children[i].getCenterRaw()));
      for (int j = 0; j < 4; ++j) {
        assertTrue(childCap.contains(children[i].getVertex(j)));
        assertTrue(childRect.contains(children[i].getVertex(j)));
        assertTrue(childRect.contains(children[i].getVertexRaw(j)));
        assertTrue(parentCap.contains(children[i].getVertex(j)));
        if (!parentRect.contains(children[i].getVertex(j))) {
          System.out.println("cell: " + cell + " i: " + i + " j: " + j);
          System.out.println("Children " + i + ": " + children[i]);
          System.out.println("Parent rect: " + parentRect);
          System.out.println("Vertex raw(j) " + children[i].getVertex(j));
          System.out.println("Latlng of vertex: " + new S2LatLng(children[i].getVertex(j)));
        }
        assertTrue(parentRect.contains(children[i].getVertex(j)));
        if (!parentRect.contains(children[i].getVertexRaw(j))) {
          System.out.println("cell: " + cell + " i: " + i + " j: " + j);
          System.out.println("Children " + i + ": " + children[i]);
          System.out.println("Parent rect: " + parentRect);
          System.out.println("Vertex raw(j) " + children[i].getVertexRaw(j));
          System.out.println("Latlng of vertex: " + new S2LatLng(children[i].getVertexRaw(j)));
        }
        assertTrue(parentRect.contains(children[i].getVertexRaw(j)));
        if (j != i) {
          // The bounding caps and rectangles should be tight enough so that they exclude at least
          // two vertices of each adjacent cell.
          int capCount = 0;
          int rectCount = 0;
          for (int k = 0; k < 4; ++k) {
            if (childCap.contains(children[j].getVertex(k))) {
              ++capCount;
            }
            if (childRect.contains(children[j].getVertexRaw(k))) {
              ++rectCount;
            }
          }
          assertTrue(capCount <= 2);
          if (childRect.latLo().radians() > -S2.M_PI_2 && childRect.latHi().radians() < S2.M_PI_2) {
            // Bounding rectangles may be too large at the poles because the pole itself has an
            // arbitrary fixed longitude.
            assertTrue(rectCount <= 2);
          }
        }
      }

      // Check all children for the first few levels, and then sample randomly. Also subdivide one
      // corner cell, one edge cell, and one center cell so that we have a better chance of sample
      // the minimum metric values. We also always subdivide the cells containing a few chosen
      // points so that we have a better chance of sampling the minimum and maximum metric values.
      // kMaxSizeUV is the absolute value of the u- and v-coordinate where the cell size at a given
      // level is maximal.
      double kMaxSizeUV = 0.3964182625366691;
      R2Vector[] specialUv = {
        new R2Vector(S2.DBL_EPSILON, S2.DBL_EPSILON), // Face center
        new R2Vector(S2.DBL_EPSILON, 1), // Edge midpoint
        new R2Vector(1, 1), // Face corner
        new R2Vector(kMaxSizeUV, kMaxSizeUV), // Largest cell area
        new R2Vector(S2.DBL_EPSILON, kMaxSizeUV)
      }; // Longest edge/diagonal
      boolean forceSubdivide = false;
      for (int k = 0; k < specialUv.length; k++) {
        if (children[i].getBoundUV().contains(specialUv[k])) {
          forceSubdivide = true;
        }
      }
      if (forceSubdivide || cell.level() < (DEBUG_MODE ? 5 : 6) || data.oneIn(DEBUG_MODE ? 5 : 4)) {
        checkSubdivide(children[i]);
      }
    }

    // Check sum of child areas equals parent area.
    //
    // For exactArea(), the best relative error we can expect is about 1e-6 because the precision of
    // the unit vector coordinates is only about 1e-15 and the edge length of a leaf cell is about
    // 1e-9.
    //
    // For approxArea(), the areas are accurate to within a few percent.
    //
    // For averageArea(), the areas themselves are not very accurate, but the average area of a
    // parent is exactly 4 times the area of a child.

    assertTrue(abs(log(exactArea / cell.exactArea())) <= abs(log1p(1e-6)));
    assertTrue(abs(log(approxArea / cell.approxArea())) <= abs(log(1.03)));
    assertTrue(abs(log(averageArea / cell.averageArea())) <= abs(log1p(1e-15)));
  }

  public void checkMinMaxAvg(
      String label,
      int level,
      double count,
      double absError,
      double minValue,
      double maxValue,
      double avgValue,
      S2.Metric minMetric,
      S2.Metric maxMetric,
      S2.Metric avgMetric) {

    // All metrics are minimums, maximums, or averages of differential quantities, and therefore
    // will not be exact for cells at any finite level. The differential minimum is always a lower
    // bound, and the maximum is always an upper bound, but these minimums and maximums may not be
    // achieved for two different reasons. First, the cells at each level are sampled and we may
    // miss the most extreme examples. Second, the actual metric for a cell is obtained by
    // integrating the differential quantity, which is not constant across the cell. Therefore cells
    // at low levels (bigger cells) have smaller variations.
    //
    // The "tolerance" below is an attempt to model both of these effects. At low levels, error is
    // dominated by the variation of differential quantities across the cells, while at high levels
    // error is dominated by the effects of random sampling.
    double tolerance =
        (maxMetric.getValue(level) - minMetric.getValue(level))
            / sqrt(min(count, 0.5 * (double) (1L << level)));
    if (tolerance == 0) {
      tolerance = absError;
    }

    double minError = minValue - minMetric.getValue(level);
    double maxError = maxMetric.getValue(level) - maxValue;
    double avgError = abs(avgMetric.getValue(level) - avgValue);
    Platform.printf(
        System.out,
        "%-10s (%6.0f samples, tolerance %8.3g) - min (%9.3g : %9.3g) "
            + "max (%9.3g : %9.3g), avg (%9.3g : %9.3g)\n",
        label,
        count,
        tolerance,
        minError / minValue,
        minError / tolerance,
        maxError / maxValue,
        maxError / tolerance,
        avgError / avgValue,
        avgError / tolerance);

    assertTrue(minMetric.getValue(level) <= minValue + absError);
    assertTrue(minMetric.getValue(level) >= minValue - tolerance);
    System.out.println("Level: " + maxMetric.getValue(level) + " max " + (maxValue + tolerance));
    assertTrue(maxMetric.getValue(level) <= maxValue + tolerance);
    assertTrue(maxMetric.getValue(level) >= maxValue - absError);
    assertDoubleNear(avgMetric.getValue(level), avgValue, 10 * tolerance);
  }

  @Test
  public void testSubdivide() {
    for (int face = 0; face < 6; ++face) {
      checkSubdivide(S2Cell.fromFace(face));
    }

    // The maximum edge *ratio* is the ratio of the longest edge of any cell to the shortest edge of
    // any cell at the same level (and similarly for the maximum diagonal ratio).
    //
    // The maximum edge *aspect* is the maximum ratio of the longest edge of a cell to the shortest
    // edge of that same cell (and similarly for the maximum diagonal aspect).

    Platform.printf(
        System.out, "Level    Area      Edge          Diag          Approx       Average\n");
    Platform.printf(
        System.out, "        Ratio  Ratio Aspect  Ratio Aspect    Min    Max    Min    Max\n");
    for (int i = 0; i <= S2CellId.MAX_LEVEL; ++i) {
      LevelStats s = levelStats.get(i);
      if (s.count > 0) {
        s.avgArea /= s.count;
        s.avgWidth /= s.count;
        s.avgEdge /= s.count;
        s.avgDiag /= s.count;
        s.avgAngleSpan /= s.count;
      }
      Platform.printf(
          System.out,
          "%5d  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
          i,
          s.maxArea / s.minArea,
          s.maxEdge / s.minEdge,
          s.maxEdgeAspect,
          s.maxDiag / s.minDiag,
          s.maxDiagAspect,
          s.minApproxRatio,
          s.maxApproxRatio,
          S2Cell.averageArea(i) / s.maxArea,
          S2Cell.averageArea(i) / s.minArea);
    }

    // Now check the validity of the S2 length and area metrics.
    for (int i = 0; i <= S2CellId.MAX_LEVEL; ++i) {
      LevelStats s = levelStats.get(i);
      if (s.count == 0) {
        continue;
      }

      Platform.printf(System.out, "Level %2d - metric (error/actual : error/tolerance)\n", i);

      // The various length calculations are only accurate to 1e-15 or so, so we need to allow for
      // this amount of discrepancy with the theoretical minimums and maximums. The area calculation
      // is accurate to about 1e-15 times the cell width.
      checkMinMaxAvg(
          "area",
          i,
          s.count,
          1e-15 * s.minWidth,
          s.minArea,
          s.maxArea,
          s.avgArea,
          MIN_AREA,
          MAX_AREA,
          AVG_AREA);
      checkMinMaxAvg(
          "width",
          i,
          s.count,
          1e-15,
          s.minWidth,
          s.maxWidth,
          s.avgWidth,
          MIN_WIDTH,
          MAX_WIDTH,
          AVG_WIDTH);
      checkMinMaxAvg(
          "edge",
          i,
          s.count,
          1e-15,
          s.minEdge,
          s.maxEdge,
          s.avgEdge,
          MIN_EDGE,
          MAX_EDGE,
          AVG_EDGE);
      checkMinMaxAvg(
          "diagonal",
          i,
          s.count,
          1e-15,
          s.minDiag,
          s.maxDiag,
          s.avgDiag,
          MIN_DIAG,
          MAX_DIAG,
          AVG_DIAG);
      checkMinMaxAvg(
          "angle span",
          i,
          s.count,
          1e-15,
          s.minAngleSpan,
          s.maxAngleSpan,
          s.avgAngleSpan,
          MIN_ANGLE_SPAN,
          MAX_ANGLE_SPAN,
          AVG_ANGLE_SPAN);

      // The aspect ratio calculations are ratios of lengths and are therefore less accurate at
      // higher subdivision levels.
      assertTrue(s.maxEdgeAspect <= MAX_EDGE_ASPECT + 1e-15 * (1 << i));
      assertTrue(s.maxDiagAspect <= MAX_DIAG_ASPECT + 1e-15 * (1 << i));
    }
  }

  @Test
  public void testCellVsLoopRectBound() {
    // This test verifies that the S2Cell and S2Loop bounds contain each other to within their
    // maximum errors.
    //
    // The S2Cell and S2Loop calculations for the latitude of a vertex can differ by up to 2 *
    // DBL_EPSILON, therefore the S2Cell bound should never exceed the S2Loop bound by more than
    // this (the reverse is not true, because the S2Loop code sometimes thinks that the maximum
    // occurs along an edge). Similarly, the longitude bounds can differ by up to 4 * DBL_EPSILON
    // since the S2Cell bound has an error of 2 * DBL_EPSILON and then expands by this amount, while
    // the S2Loop bound does no expansion at all.

    // Possible additional S2Cell error compared to S2Loop error:
    final S2LatLng kCellError = S2LatLng.fromRadians(2 * S2.DBL_EPSILON, 4 * S2.DBL_EPSILON);
    // Possible additional S2Loop error compared to S2Cell error:
    final S2LatLng kLoopError = S2EdgeUtil.RectBounder.maxErrorForTests();

    for (int iter = 0; iter < 1000; ++iter) {
      S2Cell cell = new S2Cell(data.getRandomCellId());
      S2Loop loop = new S2Loop(cell);
      S2LatLngRect cellBound = cell.getRectBound();
      S2LatLngRect loopBound = loop.getRectBound();
      assertTrue(loopBound.expanded(kCellError).contains(cellBound));
      assertTrue(cellBound.expanded(kLoopError).contains(loopBound));
    }
  }

  @Test
  public void testRectBoundIsLargeEnough() {
    // Construct many points that are nearly on an S2Cell edge, and verify that whenever the cell
    // contains a point P then its bound contains S2LatLng(P).
    for (int iter = 0; iter < 1000; /* advanced in loop below */ ) {
      S2Cell cell = new S2Cell(data.getRandomCellId());
      int i1 = data.nextInt(4);
      int i2 = i1 + 1;
      S2Point v1 = cell.getVertex(i1);
      S2Point v2 =
          data.samplePoint(S2Cap.fromAxisAngle(cell.getVertex(i2), S1Angle.radians(1e-15)));
      S2Point p = S2EdgeUtil.interpolate(data.nextDouble(), v1, v2);
      if (new S2Loop(cell).contains(p)) {
        assertTrue(cell.getRectBound().contains(new S2LatLng(p)));
        ++iter;
      }
    }
  }

  @Test
  public void testConsistentWithS2CellIdFromPoint() {
    // Construct many points that are nearly on an S2Cell edge, and verify that
    // S2Cell(S2CellId(p)).Contains(p) is always true.

    for (int iter = 0; iter < 1000; ++iter) {
      S2Cell cell = new S2Cell(data.getRandomCellId());
      int i = data.nextInt(4);
      S2Point v1 = cell.getVertex(i);
      S2Point v2 =
          data.samplePoint(S2Cap.fromAxisAngle(cell.getVertex(i + 1), S1Angle.radians(1e-15)));
      S2Point p = S2EdgeUtil.interpolate(data.nextDouble(), v1, v2);
      S2CellId cellId = S2CellId.fromPoint(p);
      S2Cell c = new S2Cell(cellId);
      assertTrue("iter " + iter + ": S2Cell.fromPoint" + p.toString() + " = " + cellId.toToken(),
          c.contains(p));
    }
  }

  @Test
  public void testAmbiguousContainsPoint() {
    // This tests a case where S2CellId returns the "wrong" cell for a point that is very close to
    // the cell edge. (testConsistentWithS2CellIdFromPoint generates more examples like this.)
    //
    // The S2Point below should have x = 0, but conversion from latlng to (x,y,z) gives x = 6.1e-17.
    // When xyz is converted to uv, this gives u = -6.1e-17.  However when converting to st, which
    // is centered at 0.5 rather than 0, the low precision bits of u are lost and we wind up with
    // s = 0.5. S2CellId.fromPoint() then chooses an arbitrary neighboring cell.
    //
    // This tests that S2Cell.contains() expands the cell bounds sufficiently so that the returned
    // cell is still considered to contain "p".
    S2Point p = S2LatLng.fromDegrees(-2, 90).toPoint();
    S2CellId cellId = S2CellId.fromPoint(p).parent(1);
    S2Cell cell = new S2Cell(cellId);
    assertTrue(cell.contains(p));
  }

  @Test
  public void testGetDistanceToPoint() {
    for (int iter = 0; iter < 1000; ++iter) {
      S2Cell cell = new S2Cell(data.getRandomCellId());
      S2Point target = data.getRandomPoint();
      S1Angle expectedToBoundary = getDistanceToCellBruteForce(cell, target).toAngle();
      S1Angle expectedToInterior = cell.contains(target) ? S1Angle.ZERO : expectedToBoundary;
      S1Angle expectedMax = getMaxDistanceToPointBruteForce(cell, target).toAngle();
      S1Angle actualToBoundary = cell.getBoundaryDistance(target).toAngle();
      S1Angle actualToInterior = cell.getDistance(target).toAngle();
      S1Angle actualMax = cell.getMaxDistance(target).toAngle();
      // The error peaks near Pi/2 for edge distance, and near Pi for vertex distance.
      assertEquals(expectedToBoundary.radians(), actualToBoundary.radians(), 1e-12);
      assertEquals(expectedToInterior.radians(), actualToInterior.radians(), 1e-12);
      assertEquals(expectedMax.radians(), actualMax.radians(), 1e-12);
      if (expectedToBoundary.radians() <= S2.M_PI / 3) {
        assertEquals(expectedToBoundary.radians(), actualToBoundary.radians(), 1e-15);
        assertEquals(expectedToInterior.radians(), actualToInterior.radians(), 1e-15);
      }
      if (expectedMax.radians() <= S2.M_PI / 3) {
        assertEquals(expectedMax.radians(), actualMax.radians(), 1e-15);
      }
    }
  }

  @Test
  public void testGetDistanceMax() {
    for (int iter = 0; iter < 1000; ++iter) {
      S2Cell cell = new S2Cell(data.getRandomCellId());
      S2Point target = data.getRandomPoint();
      S1Angle expected = getMaxDistanceToPointBruteForce(cell, target).toAngle();
      S1Angle actual = cell.getMaxDistance(target).toAngle();
      // The error has a peak near Pi/2 for edge distance, and another peak near Pi for vertex
      // distance.
      assertEquals(expected.radians(), actual.radians(), 1e-12);
      if (expected.radians() <= PI / 3) {
        assertEquals(expected.radians(), actual.radians(), 1e-15);
      }
    }
  }

  @Test
  public void testGetDistanceToCell() {
    // Adjacent S2Cells on the same face have distance exactly zero.
    S2Cell a = new S2Cell(new S2CellId(0x8ac0000000000000L)); // Level 3, corner of face
    S2Cell b = new S2Cell(new S2CellId(0x8b40000000000000L)); // Level 3, one below 'a'

    // Cells on the same face sharing an edge.
    assertEquals(S1ChordAngle.ZERO, a.getDistance(b));

    // Adjacent S2Cells on different faces should also have distance exactly zero.
    S2Cell c = new S2Cell(new S2CellId(0x0ac0000000000000L)); // Level 3, corner of adjacent face

    // Sharing an edge on the face boundary.
    assertEquals(S1ChordAngle.ZERO, a.getDistance(c));
    // Sharing a point on the face boundary.
    assertEquals(S1ChordAngle.ZERO, b.getDistance(c));
  }

  private static S1ChordAngle getDistanceToCellBruteForce(S2Cell cell, S2Point target) {
    S1ChordAngle minDistance = S1ChordAngle.INFINITY;
    for (int i = 0; i < 4; ++i) {
      S2Point a = cell.getVertex(i);
      S2Point b = cell.getVertex(i + 1);
      minDistance = S2EdgeUtil.updateMinDistance(target, a, b, minDistance);
    }
    return minDistance;
  }

  private static S1ChordAngle getMaxDistanceToPointBruteForce(S2Cell cell, S2Point target) {
    if (cell.contains(target.neg())) {
      return S1ChordAngle.STRAIGHT;
    }
    S1ChordAngle maxDistance = S1ChordAngle.NEGATIVE;
    for (int i = 0; i < 4; ++i) {
      maxDistance =
          S2EdgeUtil.updateMaxDistance(
              target, cell.getVertex(i), cell.getVertex(i + 1), maxDistance);
    }
    return maxDistance;
  }

  @Test
  public void testGetDistanceToEdge() {
    for (int iter = 0; iter < 1000; ++iter) {
      S2Cell cell = new S2Cell(data.getRandomCellId());
      MutableEdge edge = new MutableEdge();
      chooseEdgeNearCell(cell, edge);
      S1Angle expectedMin = getDistanceToEdgeBruteForce(cell, edge.a, edge.b).toAngle();
      S1Angle expectedMax = getMaxDistanceToEdgeBruteForce(cell, edge.a, edge.b).toAngle();
      S1Angle actualMin = cell.getDistanceToEdge(edge.a, edge.b).toAngle();
      S1Angle actualMax = cell.getMaxDistance(edge.a, edge.b).toAngle();
      // The error peaks near Pi/2 for edge distance, and near Pi for vertex distance.
      if (expectedMin.radians() > S2.M_PI / 2) {
        // Max error for S1ChordAngle as it approaches Pi is just less than 3e-8. Specifically,
        // expected vs. actual values *at* Pi differ by 2.98e-8.
        assertEquals(expectedMin.radians(), actualMin.radians(), 3e-8);
      } else if (expectedMin.radians() <= S2.M_PI / 3) {
        assertEquals(expectedMin.radians(), actualMin.radians(), 1e-15);
      } else {
        assertEquals(expectedMin.radians(), actualMin.radians(), 1e-12);
      }

      assertEquals(expectedMax.radians(), actualMax.radians(), 1e-12);
      if (expectedMax.radians() <= S2.M_PI / 3) {
        assertEquals(expectedMax.radians(), actualMax.radians(), 1e-15);
      }
    }
  }

  private static S1ChordAngle getDistanceToEdgeBruteForce(S2Cell cell, S2Point a, S2Point b) {
    if (cell.contains(a) || cell.contains(b)) {
      return S1ChordAngle.ZERO;
    }
    S1ChordAngle minDist = S1ChordAngle.INFINITY;
    for (int i = 0; i < 4; i++) {
      S2Point v0 = cell.getVertex(i);
      S2Point v1 = cell.getVertex(i + 1);
      // If the edge crosses through the cell, max distance is 0.
      if (S2EdgeUtil.robustCrossing(a, b, v0, v1) >= 0) {
        return S1ChordAngle.ZERO;
      }
      minDist = S2EdgeUtil.updateMinDistance(a, v0, v1, minDist);
      minDist = S2EdgeUtil.updateMinDistance(b, v0, v1, minDist);
      minDist = S2EdgeUtil.updateMinDistance(v0, a, b, minDist);
    }
    return minDist;
  }

  @Test
  public void testGetMaxDistanceToEdge() {
    // Test an edge for which its antipode crosses the cell. Validates both the standard and brute
    // force implementations for this case.
    S2Cell cell = S2Cell.fromFacePosLevel(0, 0, 20);
    S2Point a = S2EdgeUtil.interpolate(2.0, cell.getCenter(), cell.getVertex(0)).neg();
    S2Point b = S2EdgeUtil.interpolate(2.0, cell.getCenter(), cell.getVertex(2)).neg();

    S1ChordAngle actual = cell.getMaxDistance(a, b);
    S1ChordAngle expected = getMaxDistanceToEdgeBruteForce(cell, a, b);

    // Verify both the expected and actual values ~= STRAIGHT.
    double expectedRadians = S1ChordAngle.STRAIGHT.toAngle().radians();
    assertEquals(expectedRadians, expected.toAngle().radians(), 1e-15);
    assertEquals(expectedRadians, actual.toAngle().radians(), 1e-15);
  }

  private static S1ChordAngle getMaxDistanceToEdgeBruteForce(S2Cell cell, S2Point a, S2Point b) {
    // If any antipodal endpoint is within the cell, the max distance is Pi.
    if (cell.contains(a.neg()) || cell.contains(b.neg())) {
      return S1ChordAngle.STRAIGHT;
    }

    S1ChordAngle maxDist = S1ChordAngle.NEGATIVE;
    for (int i = 0; i < 4; ++i) {
      S2Point v0 = cell.getVertex(i);
      S2Point v1 = cell.getVertex(i + 1);
      // If the antipodal edge crosses through the cell, max distance is Pi.
      if (S2EdgeUtil.robustCrossing(a.neg(), b.neg(), v0, v1) >= 0) {
        return S1ChordAngle.STRAIGHT;
      }
      maxDist = S2EdgeUtil.updateMaxDistance(a, v0, v1, maxDist);
      maxDist = S2EdgeUtil.updateMaxDistance(b, v0, v1, maxDist);
      maxDist = S2EdgeUtil.updateMaxDistance(v0, a, b, maxDist);
    }
    return maxDist;
  }

  private void chooseEdgeNearCell(S2Cell cell, MutableEdge edge) {
    S2Cap cap = cell.getCapBound();
    if (data.oneIn(5)) {
      // Choose a point anywhere on the sphere.
      edge.a = data.getRandomPoint();
    } else {
      // Choose a point inside or somewhere near the cell.
      S1Angle angle = S1Angle.radians(1.5 * cap.angle().radians());
      edge.a = data.samplePoint(S2Cap.fromAxisAngle(cap.axis(), angle));
    }
    // Now choose a maximum edge length ranging from very short to very long relative to the cell
    // size, and choose the other endpoint.
    double maxLength =
        min(S2.M_PI_2, 100 * pow(1e-4, data.nextDouble()) * cap.angle().radians());
    edge.b = data.samplePoint(S2Cap.fromAxisAngle(edge.a, S1Angle.radians(maxLength)));

    if (data.oneIn(20)) {
      // Occasionally replace edge with antipodal edge.
      edge.a = edge.a.neg();
      edge.b = edge.b.neg();
    }
  }

  @Test
  public void testGetMaxDistanceToCellAntipodal() {
    S2Point p = S2LatLng.fromDegrees(0, 0).toPoint();
    S2Cell cell = new S2Cell(p);
    S2Cell antipodalCell = new S2Cell(p.neg());
    S1ChordAngle dist = cell.getMaxDistance(antipodalCell);
    assertEquals(S1ChordAngle.STRAIGHT, dist);
  }

  @Test
  public void testGetMaxDistanceToCell() {
    double tolerance = 1e-8;
    for (int i = 0; i < 1000; i++) {
      S2CellId id = data.getRandomCellId();
      S2Cell cell = new S2Cell(id);
      S2CellId testId = data.getRandomCellId();
      S2Cell testCell = new S2Cell(testId);
      S2CellId antipodalLeafId = S2CellId.fromPoint(testCell.getCenter().neg());
      S2Cell antipodalTestCell = new S2Cell(antipodalLeafId.parent(testCell.level()));

      S1ChordAngle distFromMin = S1ChordAngle.sub(
          S1ChordAngle.STRAIGHT, cell.getDistance(antipodalTestCell));
      S1ChordAngle distFromMax = cell.getMaxDistance(testCell);
      assertEquals(distFromMin.radians(), distFromMax.radians(), tolerance);
    }
  }

  @GwtIncompatible("Javascript doesn't support Java serialization.")
  @Test
  public void testS2CellSerialization() {
    doSerializationTest(new S2Cell(S2LatLng.fromDegrees(0.1, 0.2)));
  }
}
