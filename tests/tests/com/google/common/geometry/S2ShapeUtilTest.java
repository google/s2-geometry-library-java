/*
 * Copyright 2014 Google Inc.
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

import static com.google.common.collect.ImmutableSet.toImmutableSet;
import static com.google.common.geometry.S2TextFormat.makeLoop;
import static com.google.common.geometry.S2TextFormat.makePoint;
import static com.google.common.geometry.S2TextFormat.makePolygon;
import static com.google.common.geometry.S2TextFormat.makePolygonOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Lists;
import com.google.common.geometry.S2Polygon.Shape;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2Shape.ReferencePoint;
import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import com.google.common.geometry.S2ShapeUtil.RangeIterator;
import com.google.common.geometry.S2ShapeUtil.S2EdgeVectorShape;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Verifies S2ShapeUtil. */
@RunWith(JUnit4.class)
public class S2ShapeUtilTest extends GeometryTestCase {
  private final Shape poly1 = makePolygon("10:20, 90:0, 20:30").shape();
  private final Shape poly2 = makePolygon("10:20, 90:0, 20:31").shape();

  @Test
  public void testEdgeAccess() {
    S2EdgeVectorShape shape = new S2EdgeVectorShape();
    final int kNumEdges = 100;
    List<S2Edge> edges = Lists.newArrayList();
    for (int i = 0; i < kNumEdges; ++i) {
      S2Point a = data.getRandomPoint();
      S2Point b = data.getRandomPoint();
      shape.add(a, b);
      edges.add(new S2Edge(a, b));
    }
    assertEquals(kNumEdges, shape.numEdges());
    MutableEdge edge = new MutableEdge();
    for (int i = 0; i < kNumEdges; ++i) {
      shape.getEdge(i, edge);
      assertEquals(edges.get(i), new S2Edge(edge.getStart(), edge.getEnd()));
    }
  }

  @Test
  public void testSingletonConstructor() {
    S2Point a = S2Point.X_POS;
    S2Point b = S2Point.Y_POS;
    S2EdgeVectorShape shape = new S2EdgeVectorShape(a, b);
    assertEquals(1, shape.numEdges());
    MutableEdge edge = new MutableEdge();
    shape.getEdge(0, edge);
    assertEquals(a, edge.getStart());
    assertEquals(b, edge.getEnd());
  }

  @Test
  public void testEmptyShape() {
    S2Shape shape = new S2EdgeVectorShape();
    assertFalse(shape.hasInterior());
    assertEquals(0, shape.numEdges());
    assertEquals(0, shape.numChains());
    assertEquals(1, shape.dimension());
  }

  @Test
  public void testSingleEdgeShape() {
    S2Shape shape = new S2EdgeVectorShape(makePoint("0:0"), makePoint("1:1"));
    assertFalse(shape.hasInterior());
    assertEquals(1, shape.numEdges());
    checkFirstNEdges(shape, "0:0|1:1");
    assertEquals(1, shape.numChains());
    checkFirstNChainStarts(shape, 0);
    checkFirstNChainLengths(shape, 1);
    checkFirstNChainEdges(shape, 0, "0:0|1:1");
    assertEquals(1, shape.dimension());
  }

  @Test
  public void testDoubleEdgeShape() {
    S2EdgeVectorShape shape = new S2EdgeVectorShape(makePoint("0:0"), makePoint("1:1"));
    shape.add(makePoint("2:2"), makePoint("3:3"));
    assertFalse(shape.hasInterior());
    assertEquals(2, shape.numEdges());
    checkFirstNEdges(shape, "0:0|1:1, 2:2|3:3");
    assertEquals(2, shape.numChains());
    checkFirstNChainStarts(shape, 0, 1);
    checkFirstNChainLengths(shape, 1, 1);
    checkFirstNChainEdges(shape, 0, "0:0|1:1");
    checkFirstNChainEdges(shape, 1, "2:2|3:3");
    assertEquals(1, shape.dimension());
  }

  @Test
  public void testEqualsShapes() {
    assertTrue(S2ShapeUtil.equals(poly1, poly1));
    assertFalse(S2ShapeUtil.equals(poly1, poly2));
  }

  @Test
  public void testEqualsShapeLists() {
    assertFalse(S2ShapeUtil.equals(Arrays.asList(), Arrays.asList(poly1)));
    assertTrue(S2ShapeUtil.equals(Arrays.asList(poly1), Arrays.asList(poly1)));
    assertFalse(S2ShapeUtil.equals(Arrays.asList(poly1), Arrays.asList(poly2)));
    assertTrue(S2ShapeUtil.equals(Arrays.asList(poly1, poly2), Arrays.asList(poly1, poly2)));
  }

  @Test
  public void testEqualsShapeIndices() {
    assertTrue(S2ShapeUtil.equals(poly1.polygon().index, poly1.polygon().index));
    S2ShapeIndex i2 = new S2ShapeIndex();
    i2.add(poly1);
    assertTrue(S2ShapeUtil.equals(poly1.polygon().index, i2));
    assertFalse(S2ShapeUtil.equals(poly1.polygon().index, poly2.polygon().index));
    assertFalse(S2ShapeUtil.equals(poly1.polygon().index, new S2ShapeIndex()));
  }

  @Test
  public void testEqualsCells() {
    S2ShapeIndex.Cell c1 = poly1.polygon().index.iterator().entry();
    assertTrue(S2ShapeUtil.equals(c1, poly2.polygon().index.iterator().entry()));
    S2ShapeIndex i2 = new S2ShapeIndex();
    i2.add(poly1);
    i2.add(poly1);
    S2ShapeIndex.Cell c2 = i2.iterator().entry();
    assertFalse(S2ShapeUtil.equals(c1, c2));
  }

  @Test
  public void testEqualsClippedShapes() {
    S2ClippedShape c1 = poly1.polygon().index.iterator().entry().clipped(0);
    assertTrue(S2ShapeUtil.equals(c1, c1));
    S2ClippedShape c2 = new S2Polygon(S2Loop.full()).index.iterator().entry().clipped(0);
    assertFalse(S2ShapeUtil.equals(c1, c2));
  }

  @Test
  public void testReferenceEmptyPolygon() {
    assertFalse(new S2Polygon().shape().getReferencePoint().contained());
  }

  @Test
  public void testReferenceFullPolygon() {
    assertTrue(new S2Polygon(S2Loop.full()).shape().getReferencePoint().contained());
  }

  @Test
  public void testReferenceDegenerateLoops() {
    List<S2Loop> loops =
        new ArrayList<>(
            Arrays.asList(
                makeInvalidLoop("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1"),
                makeInvalidLoop("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0"),
                makeInvalidLoop("5:5, 6:6")));
    S2Polygon invalidPolygon = uncheckedCreate(() -> new S2Polygon(loops));
    assertFalse(invalidPolygon.shape().getReferencePoint().contained());
  }

  @Test
  public void testReferenceInvertedLoops() {
    S2Polygon poly = new S2Polygon();
    poly.initOriented(
        new ArrayList<>(Arrays.asList(makeLoop("1:2, 1:1, 2:2"), makeLoop("3:4, 3:3, 4:4"))));
    assertTrue(S2ShapeUtil.containsBruteForce(poly.shape(), S2.origin()));
  }

  @Test
  public void testReferencePartiallyDegenerateLoops() {
    List<S2Point> loop = new ArrayList<>();
    for (int iter = 0; iter < 100; iter++) {
      // First we construct a long convoluted edge chain that follows the S2CellId Hilbert curve.
      // At some random point along the curve, we insert a small triangular loop.
      int numVertices = 100;
      S2CellId start = data.getRandomCellId(S2CellId.MAX_LEVEL - 1);
      S2CellId end = start.advanceWrap(numVertices);
      S2CellId loopCellId = start.advanceWrap(data.random(numVertices - 2) + 1);
      List<S2Point> triangle = new ArrayList<>();
      loop.clear();
      for (S2CellId cellId = start; !cellId.equals(end); cellId = cellId.nextWrap()) {
        if (cellId.equals(loopCellId)) {
          // Insert a small triangular loop. We save the loop so that we can test whether it
          // contains the origin later.
          triangle.add(cellId.child(0).toPoint());
          triangle.add(cellId.child(1).toPoint());
          triangle.add(cellId.child(2).toPoint());
          loop.addAll(triangle);
          loop.add(cellId.child(0).toPoint());
        } else {
          loop.add(cellId.toPoint());
        }
      }
      // Now we retrace our steps, except that we skip the three edges that form the triangular
      // loop above.
      for (S2CellId cellId = end; !cellId.equals(start); cellId = cellId.prevWrap()) {
        if (cellId.equals(loopCellId)) {
          loop.add(cellId.child(0).toPoint());
        } else {
          loop.add(cellId.toPoint());
        }
      }

      // Since the S2Loop reference point and contains methods aren't tolerant of degeneracies,
      // just use it as a bag of edges, and delegate to S2ShapeUtil methods that are tolerant.
      S2Loop triangleLoop = new S2Loop(triangle);
      S2Loop degenLoop = uncheckedCreate(() -> new S2Loop(loop));
      ReferencePoint ref = S2ShapeUtil.getReferencePoint(degenLoop);
      assertEquals(S2ShapeUtil.containsBruteForce(triangleLoop, ref.point()), ref.contained());
    }
  }

  @Test
  public void testCountVertices() {
    // Test index built only out of three points.
    S2ShapeIndex index = S2TextFormat.makeIndexOrDie("1:1 | 2:2 | 3:3 # #");
    assertEquals(3, S2ShapeUtil.countVertices(index));

    // Test index built out of two points and a two edge polylines.
    index = S2TextFormat.makeIndexOrDie("1:1 | 2:2 # 3:3, 4:4, 5:5 #");
    assertEquals(5, S2ShapeUtil.countVertices(index));

    // Test index built out of two points, one two edge polylines, and four edge polygons.
    index = S2TextFormat.makeIndexOrDie("1:1 | 2:2 # 3:3, 4:4, 5:5 # 6:6, 7:7, 8:8, 9:9");
    assertEquals(9, S2ShapeUtil.countVertices(index));

    // Test that degenerate polylines count correctly.
    index = S2TextFormat.makeIndexOrDie("# 3:3, 3:3, 3:3 #");
    assertEquals(3, S2ShapeUtil.countVertices(index));

    // Test that degenerate polygons count correctly.
    index = S2TextFormat.makeIndexOrDie("# # 4:4, 4:4, 4:4, 4:4");
    assertEquals(4, S2ShapeUtil.countVertices(index));
  }

  @Test
  public void testLowerBound() {
    int[] data = {4, 4, 4, 5, 8, 8};
    assertEquals("before range of 4s", 0, S2ShapeUtil.lowerBound(0, 6, i -> data[i] < -1));
    assertEquals("starts range of 4s", 0, S2ShapeUtil.lowerBound(0, 6, i -> data[i] < 4));
    assertEquals("after range of 4s", 3, S2ShapeUtil.lowerBound(0, 6, i -> data[i] < 5));
    assertEquals("first in range of 8s", 4, S2ShapeUtil.lowerBound(0, 6, i -> data[i] < 7));
    assertEquals("after range of 8s", 6, S2ShapeUtil.lowerBound(0, 6, i -> data[i] < 10));
  }

  @Test
  public void testUpperBound() {
    int[] data = {1, 3, 5, 7, 7};
    assertEquals("after equal range", 5, S2ShapeUtil.upperBound(0, 5, i -> data[i] > 8));
    assertEquals("ends equal range", 5, S2ShapeUtil.upperBound(0, 5, i -> data[i] > 7));
    assertEquals("before equal range", 2, S2ShapeUtil.upperBound(0, 5, i -> data[i] > 3));
    assertEquals("before start", 0, S2ShapeUtil.upperBound(0, 5, i -> data[i] > 0));
  }

  @Test
  public void testRangeIteratorNext() {
    // Create an index with one point each on S2CellId faces 0, 1, and 2.
    S2ShapeIndex index = new S2ShapeIndex();
    // TODO(eengle): Should be S2TextFormat.makeIndexOrDie("0:0 | 0:90 | 90:0 # #"), but j2cl fails.
    index.add(S2Point.Shape.fromList(Arrays.asList(S2Point.X_POS, S2Point.Y_POS, S2Point.Z_POS)));
    RangeIterator<Cell> it = new RangeIterator<>(index.iterator());
    assertEquals(0, it.id().face());
    it.next();
    assertEquals(1, it.id().face());
    it.next();
    assertEquals(2, it.id().face());
    it.next();
    assertEquals(S2CellId.sentinel(), it.id());
    assertTrue(it.done());
  }

  @Test
  public void testRangeIteratorEmptyIndex() {
    // TODO(eengle): Should be S2TextFormat.makeIndexOrDie("# #"), but fails j2cl.
    S2ShapeIndex empty = new S2ShapeIndex();
    S2ShapeIndex nonEmpty = new S2ShapeIndex();
    // TODO(eengle): Should be S2TextFormat.makeIndexOrDie("0:0 # #"), but fails j2cl.
    nonEmpty.add(S2Point.Shape.singleton(S2Point.X_POS));
    RangeIterator<Cell> emptyIt = new RangeIterator<>(empty.iterator());
    RangeIterator<Cell> nonEmptyIt = new RangeIterator<>(nonEmpty.iterator());
    assertFalse(nonEmptyIt.done());
    assertTrue(emptyIt.done());

    emptyIt.seekTo(nonEmptyIt);
    assertTrue(emptyIt.done());

    emptyIt.seekBeyond(nonEmptyIt);
    assertTrue(emptyIt.done());

    emptyIt.seekTo(emptyIt);
    assertTrue(emptyIt.done());

    emptyIt.seekBeyond(emptyIt);
    assertTrue(emptyIt.done());

    nonEmptyIt.seekBeyond(nonEmptyIt);
    assertFalse(nonEmptyIt.id().isValid());
  }

  @Test
  public void testFindSelfIntersectionBasic() {
    // Coordinates are (lat,lng), which can be visualized as (y,x).
    checkHasCrossing("0:0, 0:1, 0:2, 1:2, 1:1, 1:0", false);
    checkHasCrossing("0:0, 0:1, 0:2, 1:2, 0:1, 1:0", true); // duplicate vertex
    checkHasCrossing("0:0, 0:1, 1:0, 1:1", true); // edge crossing
    checkHasCrossing("0:0, 1:1, 0:1; 0:0, 1:1, 1:0", true); // duplicate edge on different chains
    checkHasCrossing("0:0, 1:1, 0:1; 1:1, 0:0, 1:0", true); // reversed edge
    checkHasCrossing("0:0, 0:2, 2:2, 2:0; 1:1, 0:2, 3:1, 2:0", true); // vertex crossing
  }

  /**
   * Return true if any loop crosses any other loop (including vertex crossings and duplicate
   * edges), or any loop has a self-intersection (including duplicate vertices).
   */
  private static boolean hasSelfIntersection(S2ShapeIndex index) {
    S2Error error = new S2Error();
    if (S2ShapeUtil.findSelfIntersection(index, error)) {
      System.out.println("  Found self intersection: " + error);
      return true;
    }
    return false;
  }

  /**
   * This function recursively verifies that HasCrossing returns the given result for all possible
   * cyclic permutations of the loop vertices for the given set of loops.
   */
  private void checkHasCrossingPermutations(List<S2Loop> loops, int i, boolean hasCrossing) {
    if (i == loops.size()) {
      S2ShapeIndex index = new S2ShapeIndex();
      S2Polygon polygon = uncheckedCreate(() -> new S2Polygon(loops));
      index.add(polygon.shape());
      assertEquals(hasCrossing, hasSelfIntersection(index));
      polygon.release(loops);
    } else {
      S2Loop origLoop = loops.get(i);
      for (int j = 0; j < origLoop.numVertices(); ++j) {
        List<S2Point> vertices = new ArrayList<>();
        for (int k = 0; k < origLoop.numVertices(); ++k) {
          vertices.add(origLoop.vertex(j + k));
        }
        // Allow invalid loops here, as we expect edge crossings.
        loops.set(i, uncheckedCreate(() -> new S2Loop(vertices)));
        checkHasCrossingPermutations(loops, i + 1, hasCrossing);
      }
      loops.set(i, origLoop);
    }
  }

  /**
   * Given a string representing a polygon, and a boolean indicating whether this polygon has any
   * self-intersections or loop crossings, verify that hasSelfIntersection returns the expected
   * result for all possible cyclic permutations of the loop vertices.
   */
  private void checkHasCrossing(String polygonStr, boolean hasCrossing) {
    // Disable assertions so we can create invalid polygons.
    S2Polygon polygon = uncheckedCreate(() -> makePolygonOrDie(polygonStr));
    List<S2Loop> loops = new ArrayList<>();
    polygon.release(loops);
    checkHasCrossingPermutations(loops, 0, hasCrossing);
  }

  // ShapeToS2Polygon tests

  // Creates a (lax) polygon shape from the provided loops, and ensures that the S2Polygon produced
  // by shapeToS2Polygon represents the same polygon. Allows invalid S2Polygons to be created and
  // compared to the corresponding S2LaxPolygonShape.
  private void verifyShapeToS2Polygon(
      List<List<S2Point>> loops, int expectedNumLoops, int expectedNumVertices) {
    S2LaxPolygonShape laxPolygon = S2LaxPolygonShape.create(loops);
    S2Polygon polygon = uncheckedCreate(() -> S2ShapeUtil.shapeToS2Polygon(laxPolygon));

    assertEquals(expectedNumLoops, polygon.numLoops());
    assertEquals(expectedNumVertices, polygon.getNumVertices());

    // S2Polygon init methods may reorder loops.
    var actualLoops =
        polygon.getLoops().stream().map(S2Loop::orientedVertices).collect(toImmutableSet());
    assertEquals(ImmutableSet.copyOf(loops), actualLoops);
  }

  @Test
  public void polygonWithHoleToS2Polygon() {
    // a polygon with one shell and one hole
    var shell = S2TextFormat.parsePointsOrDie("0:0, 0:10, 10:10, 10:0");
    var hole = S2TextFormat.parsePointsOrDie("4:4, 6:4, 6:6, 4:6");
    verifyShapeToS2Polygon(ImmutableList.of(shell, hole), 2, 8);
  }

  @Test
  public void multiPolygonToS2Polygon() {
    // a polygon with multiple shells
    var shell1 = S2TextFormat.parsePointsOrDie("0:0, 0:2, 2:2, 2:0");
    var shell2 = S2TextFormat.parsePointsOrDie("0:4, 0:6, 3:6");
    verifyShapeToS2Polygon(ImmutableList.of(shell1, shell2), 2, 7);
  }

  @Test
  public void twoHolesToS2Polygon() {
    // a polygon shell with two holes
    var shell = S2TextFormat.parsePointsOrDie("0:0, 0:10, 10:10, 10:0");
    var hole1 = S2TextFormat.parsePointsOrDie("1:1, 3:3, 1:3");
    var hole2 = S2TextFormat.parsePointsOrDie("2:6, 4:7, 2:8");
    verifyShapeToS2Polygon(ImmutableList.of(shell, hole1, hole2), 3, 10);
  }

  @Test
  public void intersectingLoopsToS2Polygon() {
    // TODO(user): This will start failing when S2Polygon asserts validity, update
    // this test and the shapeToS2Polygon Javadoc when that happens.
    var shell1 = S2TextFormat.parsePointsOrDie("0:0, 0:2, 2:2, 2:0");
    var shell2 = S2TextFormat.parsePointsOrDie("1:1, 1:3, 3:3, 3:1");
    verifyShapeToS2Polygon(ImmutableList.of(shell1, shell2), 2, 8);
  }

  @Test
  public void fullPolygonToS2Polygon() {
    // verify that a full polygon is converted correctly
    var laxPolygon = S2TextFormat.makeLaxPolygonOrDie("full");
    var polygon = S2ShapeUtil.shapeToS2Polygon(laxPolygon);
    assertEquals(1, polygon.numLoops());
    assertEquals(1, polygon.getNumVertices());
    assertTrue(polygon.isFull());
    var loops = S2Loop.full().chains();
    for (int i = 0; i < polygon.numLoops(); i++) {
      S2Loop loop = polygon.loop(i);
      for (int j = 0; j < loop.numVertices(); j++) {
        assertEquals(loop.orientedVertex(j), loops.get(i).get(j));
      }
    }
  }
}
