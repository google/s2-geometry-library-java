/*
 * Copyright 2026 Google Inc.
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

import static com.google.common.geometry.S2TextFormat.makePolygonOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2BuilderSnapFunctions.IdentitySnapFunction;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for {@link S2MultiOperation}. */
@RunWith(JUnit4.class)
public final class S2MultiOperationTest extends GeometryTestCase {

  private static final S1Angle TOLERANCE = S1Angle.degrees(1e-3);

  @Test
  public void testDisjointPolygons() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:1, 1:1, 1:0");
    S2Polygon p2 = makePolygonOrDie("2:2, 2:3, 3:3, 3:2");

    S2Polygon result = S2MultiOperation.union(p1, p2);

    S2Polygon expected = makePolygonOrDie("0:0, 0:1, 1:1, 1:0; 2:2, 2:3, 3:3, 3:2");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testOverlappingPolygons() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:2, 2:2, 2:0");
    S2Polygon p2 = makePolygonOrDie("1:1, 1:3, 3:3, 3:1");

    S2Polygon result = S2MultiOperation.union(p1, p2);

    S2Polygon expected = makePolygonOrDie("0:0, 0:2, 1:2, 1:3, 3:3, 3:1, 2:1, 2:0");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testPolygonsWithHoles() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:3, 3:3, 3:0; 1:1, 1:2, 2:2, 2:1");
    S2Polygon p2 = makePolygonOrDie("1.2:1.2, 1.2:1.8, 1.8:1.8, 1.8:1.2");

    S2Polygon result = S2MultiOperation.union(p1, p2);

    S2Polygon expected =
        makePolygonOrDie(
            "0:0, 0:3, 3:3, 3:0; 1:1, 1:2, 2:2, 2:1; 1.2:1.2, 1.2:1.8, 1.8:1.8, 1.8:1.2");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testHoleFilledByOverlap() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:3, 3:3, 3:0; 1:1, 1:2, 2:2, 2:1");
    S2Polygon p2 = makePolygonOrDie("1.5:1.5, 1.5:2.5, 2.5:2.5, 2.5:1.5");

    S2Polygon result = S2MultiOperation.union(p1, p2);

    S2Polygon expected =
        makePolygonOrDie("0:0, 0:3, 3:3, 3:0; 1:1, 1:2, 1.5:2, 1.5:1.5, 2:1.5, 2:1");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testBulkGridUnion() {
    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(S2MultiOperation.OpType.UNION, layer);
    for (int i = 10; i < 15; ++i) {
      for (int j = 0; j < 5; ++j) {
        op.add(makePolygonOrDie(unitSquareAtCornerString(i, j)));
      }
    }
    op.build();

    String expectedStr =
        "10:0, 10:1, 10:2, 10:3, 10:4, 10:5, " + "15:5, " + "15:4, 15:3, 15:2, 15:1, 15:0";
    S2Polygon expected = makePolygonOrDie(expectedStr);
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testGridUnionWithHole() {
    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(S2MultiOperation.OpType.UNION, layer);
    for (int i = 10; i < 13; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (i == 11 && j == 1) {
          continue;
        }
        op.add(makePolygonOrDie(unitSquareAtCornerString(i, j)));
      }
    }
    op.build();

    String expectedStr =
        "10:0, 10:1, 10:2, 10:3, " + "13:3, " + "13:2, 13:1, 13:0; " + "11:1, 11:2, 12:2, 12:1";
    S2Polygon expected = makePolygonOrDie(expectedStr);
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testIntersectionOverlapping() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:2, 2:2, 2:0");
    S2Polygon p2 = makePolygonOrDie("1:1, 1:3, 3:3, 3:1");

    S2Polygon result = S2MultiOperation.intersection(p1, p2);

    S2Polygon expected = makePolygonOrDie("1:1, 1:2, 2:2, 2:1");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testIntersectionContained() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:3, 3:3, 3:0");
    S2Polygon p2 = makePolygonOrDie("1:1, 1:2, 2:2, 2:1");

    S2Polygon result = S2MultiOperation.intersection(p1, p2);

    S2Polygon expected = makePolygonOrDie("1:1, 1:2, 2:2, 2:1");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testIntersectionThreePolygons() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:3, 3:3, 3:0");
    S2Polygon p2 = makePolygonOrDie("1:1, 1:4, 4:4, 4:1");
    S2Polygon p3 = makePolygonOrDie("2:2, 2:5, 5:5, 5:2");

    S2Polygon result = S2MultiOperation.intersection(p1, p2, p3);

    S2Polygon expected = makePolygonOrDie("2:2, 2:3, 3:3, 3:2");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testOddContained() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:3, 3:3, 3:0");
    S2Polygon p2 = makePolygonOrDie("1:1, 1:2, 2:2, 2:1");

    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(S2MultiOperation.OpType.SYMMETRIC_DIFFERENCE, layer);
    op.add(p1);
    op.add(p2);
    op.build();

    S2Polygon expected = makePolygonOrDie("0:0, 0:3, 3:3, 3:0; 1:1, 1:2, 2:2, 2:1");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testDifference() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:3, 3:3, 3:0");
    S2Polygon p2 = makePolygonOrDie("1:1, 1:4, 4:4, 4:1");

    S2Polygon result = S2MultiOperation.difference(p1, Arrays.asList(p2));

    S2Polygon expected = makePolygonOrDie("0:0, 0:3, 1:3, 1:1, 3:1, 3:0");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testDifferenceMultipleSubtrahends() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:5, 5:5, 5:0");
    S2Polygon p2 = makePolygonOrDie("1:1, 1:2, 2:2, 2:1");
    S2Polygon p3 = makePolygonOrDie("3:3, 3:4, 4:4, 4:3");

    S2Polygon result = S2MultiOperation.difference(p1, Arrays.asList(p2, p3));

    S2Polygon expected =
        makePolygonOrDie("0:0, 0:5, 5:5, 5:0; 1:1, 1:2, 2:2, 2:1; 3:3, 3:4, 4:4, 4:3");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testDifferenceIntersectingSubtrahends() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:5, 5:5, 5:0");
    S2Polygon p2 = makePolygonOrDie("1:1, 1:3, 3:3, 3:1");
    S2Polygon p3 = makePolygonOrDie("2:2, 2:4, 4:4, 4:2");

    S2Polygon result = S2MultiOperation.difference(p1, Arrays.asList(p2, p3));

    S2Polygon expected =
        makePolygonOrDie("0:0, 0:5, 5:5, 5:0; 1:1, 3:1, 3:2, 4:2, 4:4, 2:4, 2:3, 1:3");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testClearAndReuse() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:1, 1:1, 1:0");
    S2Polygon p2 = makePolygonOrDie("2:2, 2:3, 3:3, 3:2");
    S2Polygon p3 = makePolygonOrDie("4:4, 4:5, 5:5, 5:4");

    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(S2MultiOperation.OpType.UNION, layer);
    op.add(p1);
    op.build();

    assertTrue(result.approxEquals(p1, TOLERANCE));

    op.clear();

    op.add(p2);
    op.add(p3);
    op.build();

    S2Polygon expected = S2MultiOperation.union(p2, p3);
    assertTrue(result.approxEquals(expected, TOLERANCE));
    assertFalse(result.contains(S2TextFormat.makePointOrDie("0.5:0.5")));
  }

  @Test
  public void testUnionPolygonsEnclosingOrigin() {
    S2Polygon p1 = makePolygonOrDie("80:0, 80:120, 80:-120");
    S2Polygon p2 = makePolygonOrDie("70:0, 70:120, 70:-120");

    assertTrue(p1.contains(S2.origin()));
    assertTrue(p2.contains(S2.origin()));

    S2Polygon result = S2MultiOperation.union(p1, p2);

    assertTrue(result.approxEquals(p2, TOLERANCE));
  }

  @Test
  public void testUnionWithComplement() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:1, 1:1, 1:0");
    S2Polygon p2 = makePolygonOrDie("2:2, 2:3, 3:3, 3:2");

    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(S2MultiOperation.OpType.UNION, layer);
    op.add(p1);
    op.addComplement(p2);
    op.build();

    S2Polygon expected = new S2Polygon();
    S2PolygonLayer expectedLayer = new S2PolygonLayer();
    expectedLayer.setPolygon(expected);
    S2MultiOperation expectedOp =
        new S2MultiOperation(S2MultiOperation.OpType.UNION, expectedLayer);
    expectedOp.addComplement(p2);
    expectedOp.build();

    assertTrue(result.approxEquals(expected, TOLERANCE));
    assertTrue(result.contains(S2TextFormat.makePointOrDie("0.5:0.5")));
    assertTrue(result.contains(S2TextFormat.makePointOrDie("5:5")));
    assertFalse(result.contains(S2TextFormat.makePointOrDie("2.5:2.5")));
  }

  @Test
  public void testIntersectionWithComplement() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:3, 3:3, 3:0");
    S2Polygon p2 = makePolygonOrDie("1:1, 1:2, 2:2, 2:1");

    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(S2MultiOperation.OpType.INTERSECTION, layer);
    op.add(p1);
    op.addComplement(p2);
    op.build();

    S2Polygon expected = makePolygonOrDie("0:0, 0:3, 3:3, 3:0; 1:1, 1:2, 2:2, 2:1");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testSetOptionsAfterAddingThrows() {
    S2Polygon p = makePolygonOrDie("0:0, 0:1, 1:1, 1:0");
    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(S2MultiOperation.OpType.UNION, layer);
    op.add(p);
    S2WindingOperation.Options options = new S2WindingOperation.Options();
    assertThrows(IllegalStateException.class, () -> op.setOptions(options));
  }

  @Test
  public void testSetOptionsBeforeAddingWorks() {
    S2Polygon p = makePolygonOrDie("0:0, 0:1, 1:1, 1:0");
    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(S2MultiOperation.OpType.UNION, layer);
    op.setOptions(new S2WindingOperation.Options());
    op.add(p);
    op.build();
    assertTrue(result.approxEquals(p, TOLERANCE));
  }

  @Test
  public void testAddNon2DShapeThrows() {
    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(S2MultiOperation.OpType.UNION, layer);

    S2Polyline polyline = S2TextFormat.makePolylineOrDie("0:0, 1:1");
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(polyline);

    assertThrows(IllegalArgumentException.class, () -> op.add(index));
  }

  @Test
  public void testAddPolygonDoesNotBuildColdIndexAndHandlesOriginContainment() {
    // S2.origin() is at LatLng (0, ~90E).
    S2Polygon nonOriginPoly = makePolygonOrDie("0:0, 0:2, 2:2, 2:0");
    S2Polygon originContainingShell = makePolygonOrDie("-10:80, -10:100, 10:100, 10:80");
    S2Polygon originHoleSubtrahend = makePolygonOrDie("-5:85, -5:95, 5:95, 5:85");
    S2Polygon shellWithOriginInHole =
        makePolygonOrDie("-20:70, -20:110, 20:110, 20:70; -5:85, -5:95, 5:95, 5:85");

    assertFalse(nonOriginPoly.index().isFresh());
    assertFalse(originContainingShell.index().isFresh());
    assertFalse(originHoleSubtrahend.index().isFresh());
    assertFalse(shellWithOriginInHole.index().isFresh());

    S2Polygon emptyPoly = new S2Polygon();
    assertTrue(emptyPoly.isEmpty());

    // 1. Difference where both minuend and subtrahend contain S2.origin(), plus an empty polygon
    // complement (whose complement is the full sphere).
    S2Polygon diffResult = new S2Polygon();
    S2PolygonLayer diffLayer = new S2PolygonLayer();
    diffLayer.setPolygon(diffResult);
    S2MultiOperation diffOp = new S2MultiOperation(S2MultiOperation.OpType.INTERSECTION, diffLayer);
    diffOp.add(originContainingShell);
    diffOp.addComplement(originHoleSubtrahend);
    diffOp.addComplement(emptyPoly);
    diffOp.build();

    // Verify cold polygon indexes were never built.
    assertFalse(originContainingShell.index().isFresh());
    assertFalse(originHoleSubtrahend.index().isFresh());

    S2Polygon expectedDiff =
        makePolygonOrDie("-10:80, -10:100, 10:100, 10:80; -5:85, -5:95, 5:95, 5:85");
    assertTrue(diffResult.approxEquals(expectedDiff, TOLERANCE));

    // 2. Union of a multi-loop polygon (whose hole contains S2.origin()), a disjoint polygon, and
    // an empty polygon.
    S2Polygon unionResult = new S2Polygon();
    S2PolygonLayer unionLayer = new S2PolygonLayer();
    unionLayer.setPolygon(unionResult);
    S2MultiOperation unionOp = new S2MultiOperation(S2MultiOperation.OpType.UNION, unionLayer);
    unionOp.add(emptyPoly);
    unionOp.add(shellWithOriginInHole);
    unionOp.add(nonOriginPoly);
    unionOp.build();

    assertFalse(shellWithOriginInHole.index().isFresh());
    assertFalse(nonOriginPoly.index().isFresh());

    S2Polygon expectedUnion = new S2Polygon();
    expectedUnion.initToUnion(shellWithOriginInHole, nonOriginPoly);
    assertTrue(unionResult.approxEquals(expectedUnion, TOLERANCE));

    // 3. Intersecting a non-empty polygon with an empty polygon must return an empty polygon.
    S2Polygon emptyIntersection = S2MultiOperation.intersection(nonOriginPoly, emptyPoly);
    assertTrue(emptyIntersection.isEmpty());

    // 4. Subtracting a non-empty subtrahend from an empty minuend must return an empty polygon.
    S2Polygon emptyMinuendDiff = new S2Polygon();
    S2PolygonLayer emptyMinuendLayer = new S2PolygonLayer();
    emptyMinuendLayer.setPolygon(emptyMinuendDiff);
    S2MultiOperation emptyMinuendOp =
        new S2MultiOperation(S2MultiOperation.OpType.INTERSECTION, emptyMinuendLayer);
    emptyMinuendOp.add(emptyPoly);
    emptyMinuendOp.addComplement(nonOriginPoly);
    emptyMinuendOp.build();
    assertTrue(emptyMinuendDiff.isEmpty());

    // 5. Union with the complement of an empty polygon (the full sphere) must return a full
    // polygon.
    S2Polygon fullUnionResult = new S2Polygon();
    S2PolygonLayer fullUnionLayer = new S2PolygonLayer();
    fullUnionLayer.setPolygon(fullUnionResult);
    S2MultiOperation fullUnionOp =
        new S2MultiOperation(S2MultiOperation.OpType.UNION, fullUnionLayer);
    fullUnionOp.add(nonOriginPoly);
    fullUnionOp.addComplement(emptyPoly);
    fullUnionOp.build();
    assertTrue(fullUnionResult.isFull());

    // 6. Symmetric difference with the complement of an empty polygon must return the complement
    // of the other polygon.
    S2Polygon symDiffComplementResult = new S2Polygon();
    S2PolygonLayer symDiffComplementLayer = new S2PolygonLayer();
    symDiffComplementLayer.setPolygon(symDiffComplementResult);
    S2MultiOperation symDiffComplementOp =
        new S2MultiOperation(S2MultiOperation.OpType.SYMMETRIC_DIFFERENCE, symDiffComplementLayer);
    symDiffComplementOp.add(nonOriginPoly);
    symDiffComplementOp.addComplement(emptyPoly);
    symDiffComplementOp.build();
    assertTrue(symDiffComplementResult.contains(S2.origin()));
    assertFalse(symDiffComplementResult.contains(S2TextFormat.makePointOrDie("1:1")));
  }

  @Test
  public void testAddShapeIndexColdVsWarmLaxPolygonPaths() {
    S2Polygon poly1 =
        makePolygonOrDie("70:0, 70:120, 70:-120"); // Encloses S2.origin() (~89.4N, 165.5E)
    S2Polygon poly2 = makePolygonOrDie("0:0, 0:2, 2:2, 2:0"); // Does not contain S2.origin()
    assertTrue(poly1.contains(S2.origin()));
    assertFalse(poly2.contains(S2.origin()));

    S2LaxPolygonShape lax1 = S2LaxPolygonShape.create(poly1);
    S2LaxPolygonShape lax2 = S2LaxPolygonShape.create(poly2);

    // Path A: Cold S2ShapeIndex (!index.isFresh()) uses brute-force without building index
    S2ShapeIndex coldIndex1 = new S2ShapeIndex();
    coldIndex1.add(lax1);
    S2ShapeIndex coldIndex2 = new S2ShapeIndex();
    coldIndex2.add(lax2);
    assertFalse(coldIndex1.isFresh());
    assertFalse(coldIndex2.isFresh());

    S2Polygon coldResult = new S2Polygon();
    S2PolygonLayer coldLayer = new S2PolygonLayer();
    coldLayer.setPolygon(coldResult);
    S2MultiOperation coldOp = new S2MultiOperation(S2MultiOperation.OpType.UNION, coldLayer);
    coldOp.add(coldIndex1);
    coldOp.add(coldIndex2);
    coldOp.build();

    assertFalse(
        "Cold S2ShapeIndex should not be built by S2MultiOperation.add", coldIndex1.isFresh());
    assertFalse(
        "Cold S2ShapeIndex should not be built by S2MultiOperation.add", coldIndex2.isFresh());
    assertTrue(coldResult.contains(S2.origin()));

    // Path B: Warm S2ShapeIndex (index.isFresh()) uses S2ContainsPointQuery
    S2ShapeIndex warmIndex1 = new S2ShapeIndex();
    warmIndex1.add(lax1);
    warmIndex1.applyUpdates();
    S2ShapeIndex warmIndex2 = new S2ShapeIndex();
    warmIndex2.add(lax2);
    warmIndex2.applyUpdates();
    assertTrue(warmIndex1.isFresh());
    assertTrue(warmIndex2.isFresh());

    S2Polygon warmResult = new S2Polygon();
    S2PolygonLayer warmLayer = new S2PolygonLayer();
    warmLayer.setPolygon(warmResult);
    S2MultiOperation warmOp = new S2MultiOperation(S2MultiOperation.OpType.UNION, warmLayer);
    warmOp.add(warmIndex1);
    warmOp.add(warmIndex2);
    warmOp.build();
    assertTrue(warmResult.contains(S2.origin()));

    S2Polygon expected = new S2Polygon();
    expected.initToUnion(poly1, poly2);
    assertTrue(coldResult.approxEquals(expected, TOLERANCE));
    assertTrue(warmResult.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testAddComplementShapeIndex() {
    S2Polygon outer = makePolygonOrDie("0:0, 0:3, 3:3, 3:0");
    S2Polygon hole = makePolygonOrDie("1:1, 1:2, 2:2, 2:1");
    S2ShapeIndex holeIndex = new S2ShapeIndex();
    holeIndex.add(S2LaxPolygonShape.create(hole));

    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(S2MultiOperation.OpType.INTERSECTION, layer);
    op.add(outer);
    op.addComplement(holeIndex);
    op.build();

    S2Polygon expected = makePolygonOrDie("0:0, 0:3, 3:3, 3:0; 1:1, 2:1, 2:2, 1:2");
    assertTrue(result.approxEquals(expected, TOLERANCE));
    assertTrue(result.contains(S2TextFormat.makePointOrDie("0.5:0.5")));
    assertFalse(result.contains(S2TextFormat.makePointOrDie("1.5:1.5")));
  }

  @Test
  public void testAddComplementShapeIndexEnclosingOrigin() {
    S2Polygon shell = makePolygonOrDie("60:0, 60:120, 60:-120");
    S2Polygon originHole = makePolygonOrDie("80:0, 80:120, 80:-120");
    assertTrue(shell.contains(S2.origin()));
    assertTrue(originHole.contains(S2.origin()));

    S2ShapeIndex originHoleIndex = new S2ShapeIndex();
    originHoleIndex.add(S2LaxPolygonShape.create(originHole));

    S2Polygon diffResult = new S2Polygon();
    S2PolygonLayer diffLayer = new S2PolygonLayer();
    diffLayer.setPolygon(diffResult);
    S2MultiOperation diffOp = new S2MultiOperation(S2MultiOperation.OpType.INTERSECTION, diffLayer);
    diffOp.add(shell);
    diffOp.addComplement(originHoleIndex);
    diffOp.build();

    assertFalse(diffResult.contains(S2.origin()));
    assertTrue(diffResult.contains(S2TextFormat.makePointOrDie("70:0")));
  }

  private static String unitSquareAtCornerString(int i, int j) {
    return i + ":" + j + ", " + i + ":" + (j + 1) + ", " + (i + 1) + ":" + (j + 1) + ", " + (i + 1)
        + ":" + j;
  }

  @Test
  public void testDefaultOptionsAndOpType() {
    S2WindingOperation.Options defaultOptions = S2MultiOperation.defaultOptions();
    assertEquals(S2EdgeUtil.INTERSECTION_MERGE_RADIUS, defaultOptions.snapFunction().snapRadius());

    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(S2MultiOperation.OpType.UNION, layer);
    assertEquals(S2MultiOperation.OpType.UNION, op.opType());
  }

  @Test
  public void testDefaultSnappingBehavior() {
    S1Angle mergeRadius = S2EdgeUtil.INTERSECTION_MERGE_RADIUS;
    double delta = 0.5 * mergeRadius.radians(); // Strictly within INTERSECTION_MERGE_RADIUS

    // Diamond 1 centered around longitude 0.000, with its eastern corner at longitude 0.005.
    S2Point p1North = S2LatLng.fromRadians(0.01, 0.0).toPoint();
    S2Point p1West = S2LatLng.fromRadians(0.0, -0.01).toPoint();
    S2Point p1South = S2LatLng.fromRadians(-0.01, 0.0).toPoint();
    S2Point p1East = S2LatLng.fromRadians(0.0, 0.005).toPoint();

    // Diamond 2 centered around longitude 0.010, with its western corner at longitude 0.005 +
    // delta.
    S2Point p2West = S2LatLng.fromRadians(0.0, 0.005 + delta).toPoint();
    S2Point p2North = S2LatLng.fromRadians(0.01, 0.01 + delta).toPoint();
    S2Point p2East = S2LatLng.fromRadians(0.0, 0.015 + delta).toPoint();
    S2Point p2South = S2LatLng.fromRadians(-0.01, 0.01 + delta).toPoint();

    S2Polygon d1 = new S2Polygon(new S2Loop(Arrays.asList(p1North, p1West, p1South, p1East)));
    S2Polygon d2 = new S2Polygon(new S2Loop(Arrays.asList(p2West, p2South, p2East, p2North)));

    // Default union merges nearby corners within INTERSECTION_MERGE_RADIUS into a shared vertex.
    S2Polygon multiOpDefault = S2MultiOperation.union(d1, d2);
    assertEquals(2, multiOpDefault.numLoops());
    assertTrue(
        "Default S2MultiOperation.union should merge nearby corners into a shared vertex",
        hasSharedVertexBetweenLoops(multiOpDefault));

    S2Polygon multiOpDefaultIterable = S2MultiOperation.union(Arrays.asList(d1, d2));
    assertEquals(2, multiOpDefaultIterable.numLoops());
    assertTrue(
        "Default S2MultiOperation.union(Iterable) should merge nearby corners into a shared vertex",
        hasSharedVertexBetweenLoops(multiOpDefaultIterable));

    // Explicit zero snap radius: non-coincident vertices remain distinct.
    S2WindingOperation.Options zeroSnapOptions =
        new S2WindingOperation.Options(new IdentitySnapFunction(S1Angle.ZERO));
    S2Polygon multiOpZeroSnap = S2MultiOperation.union(Arrays.asList(d1, d2), zeroSnapOptions);
    assertEquals(2, multiOpZeroSnap.numLoops());
    assertFalse(
        "Explicit snapRadius=0 should keep non-coincident vertices distinct",
        hasSharedVertexBetweenLoops(multiOpZeroSnap));
  }

  private static boolean hasSharedVertexBetweenLoops(S2Polygon polygon) {
    for (int i = 0; i < polygon.numLoops(); i++) {
      S2Loop loopI = polygon.loop(i);
      for (int j = i + 1; j < polygon.numLoops(); j++) {
        S2Loop loopJ = polygon.loop(j);
        for (int vi = 0; vi < loopI.numVertices(); vi++) {
          for (int vj = 0; vj < loopJ.numVertices(); vj++) {
            if (loopI.vertex(vi).equalsPoint(loopJ.vertex(vj))) {
              return true;
            }
          }
        }
      }
    }
    return false;
  }

  @Test
  public void testZeroPolygons() {
    S2Polygon unionResult = S2MultiOperation.union();
    assertTrue(unionResult.isEmpty());

    S2Polygon unionIterableResult = S2MultiOperation.union(ImmutableList.of());
    assertTrue(unionIterableResult.isEmpty());

    S2Polygon intersectionResult = S2MultiOperation.intersection();
    assertTrue(intersectionResult.isFull());

    S2Polygon intersectionIterableResult = S2MultiOperation.intersection(ImmutableList.of());
    assertTrue(intersectionIterableResult.isFull());
  }

  @Test
  public void testSinglePolygon() {
    S2Polygon p = makePolygonOrDie("0:0, 0:2, 2:2, 2:0");

    S2Polygon unionResult = S2MultiOperation.union(p);
    assertTrue(unionResult.approxEquals(p, TOLERANCE));

    S2Polygon unionIterableResult = S2MultiOperation.union(ImmutableList.of(p));
    assertTrue(unionIterableResult.approxEquals(p, TOLERANCE));

    S2Polygon intersectionResult = S2MultiOperation.intersection(p);
    assertTrue(intersectionResult.approxEquals(p, TOLERANCE));

    S2Polygon intersectionIterableResult = S2MultiOperation.intersection(ImmutableList.of(p));
    assertTrue(intersectionIterableResult.approxEquals(p, TOLERANCE));
  }

  @Test
  public void testDifferenceZeroOrOneSubtrahends() {
    S2Polygon minuend = makePolygonOrDie("0:0, 0:3, 3:3, 3:0");
    S2Polygon subtrahend = makePolygonOrDie("1:1, 1:2, 2:2, 2:1");

    S2Polygon diffEmpty = S2MultiOperation.difference(minuend, ImmutableList.of());
    assertTrue(diffEmpty.approxEquals(minuend, TOLERANCE));

    S2Polygon diffSingle = S2MultiOperation.difference(minuend, ImmutableList.of(subtrahend));
    S2Polygon expected = makePolygonOrDie("0:0, 0:3, 3:3, 3:0; 1:1, 1:2, 2:2, 2:1");
    assertTrue(diffSingle.approxEquals(expected, TOLERANCE));
  }

  @Test
  public void testPermutationInvariance() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:2, 2:2, 2:0");
    S2Polygon p2 = makePolygonOrDie("1:1, 1:3, 3:3, 3:1");
    S2Polygon p3 = makePolygonOrDie("2:2, 2:4, 4:4, 4:2");
    S2Polygon p4 = makePolygonOrDie("5:5, 5:6, 6:6, 6:5");
    S2Polygon p5 = makePolygonOrDie("1.5:1.5, 1.5:2.5, 2.5:2.5, 2.5:1.5");

    List<S2Polygon> inputs = Arrays.asList(p1, p2, p3, p4, p5);
    S2Polygon baseline = S2MultiOperation.union(inputs);

    Random rand = new Random(42);
    List<S2Polygon> permutedInputs = new ArrayList<>(inputs);
    for (int i = 0; i < 20; i++) {
      Collections.shuffle(permutedInputs, rand);
      S2Polygon permutedResult = S2MultiOperation.union(permutedInputs);
      assertEquals(
          "S2MultiOperation.union must produce identical polygon results under input permutation",
          baseline,
          permutedResult);
    }
  }

  @Test
  public void testSingleUseIterableUnion() {
    S2Polygon p1 = makePolygonOrDie("0:0, 0:1, 1:1, 1:0");
    S2Polygon p2 = makePolygonOrDie("2:2, 2:3, 3:3, 3:2");

    Iterator<S2Polygon> singleUseIterator = Arrays.asList(p1, p2).iterator();
    Iterable<S2Polygon> singleUseIterable = () -> singleUseIterator;

    S2Polygon result = S2MultiOperation.union(singleUseIterable);
    S2Polygon expected = makePolygonOrDie("0:0, 0:1, 1:1, 1:0; 2:2, 2:3, 3:3, 3:2");
    assertTrue(result.approxEquals(expected, TOLERANCE));
  }
}
