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

import com.google.common.annotations.GwtCompatible;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.geometry.S2Polygon.Shape;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2Shape.ReferencePoint;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import com.google.common.geometry.S2ShapeUtil.S2EdgeVectorShape;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/** Verifies S2ShapeUtil. */
@GwtCompatible
public class S2ShapeUtilTest extends GeometryTestCase {
  private final Shape poly1 = makePolygon("10:20, 90:0, 20:30").shape();
  private final Shape poly2 = makePolygon("10:20, 90:0, 20:31").shape();

  public void testEdgeAccess() {
    S2EdgeVectorShape shape = new S2EdgeVectorShape();
    final int kNumEdges = 100;
    List<S2Edge> edges = Lists.newArrayList();
    for (int i = 0; i < kNumEdges; ++i) {
      S2Point a = randomPoint();
      S2Point b = randomPoint();
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

  public void testEmptyShape() {
    S2Shape shape = new S2EdgeVectorShape();
    assertFalse(shape.hasInterior());
    assertEquals(0, shape.numEdges());
    assertEquals(0, shape.numChains());
    assertEquals(1, shape.dimension());
  }

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

  public void testEqualsShapes() {
    assertTrue(S2ShapeUtil.equals(poly1, poly1));
    assertFalse(S2ShapeUtil.equals(poly1, poly2));
  }

  public void testEqualsShapeLists() {
    assertFalse(S2ShapeUtil.equals(Arrays.asList(), Arrays.asList(poly1)));
    assertTrue(S2ShapeUtil.equals(Arrays.asList(poly1), Arrays.asList(poly1)));
    assertFalse(S2ShapeUtil.equals(Arrays.asList(poly1), Arrays.asList(poly2)));
    assertTrue(S2ShapeUtil.equals(Arrays.asList(poly1, poly2), Arrays.asList(poly1, poly2)));
  }

  public void testEqualsShapeIndices() {
    assertTrue(S2ShapeUtil.equals(poly1.polygon().index, poly1.polygon().index));
    S2ShapeIndex i2 = new S2ShapeIndex();
    i2.add(poly1);
    assertTrue(S2ShapeUtil.equals(poly1.polygon().index, i2));
    assertFalse(S2ShapeUtil.equals(poly1.polygon().index, poly2.polygon().index));
    assertFalse(S2ShapeUtil.equals(poly1.polygon().index, new S2ShapeIndex()));
  }

  public void testEqualsCells() {
    S2ShapeIndex.Cell c1 = poly1.polygon().index.iterator().entry();
    assertTrue(S2ShapeUtil.equals(c1, poly2.polygon().index.iterator().entry()));
    S2ShapeIndex i2 = new S2ShapeIndex();
    i2.add(poly1);
    i2.add(poly1);
    S2ShapeIndex.Cell c2 = i2.iterator().entry();
    assertFalse(S2ShapeUtil.equals(c1, c2));
  }

  public void testEqualsClippedShapes() {
    S2ClippedShape c1 = poly1.polygon().index.iterator().entry().clipped(0);
    assertTrue(S2ShapeUtil.equals(c1, c1));
    S2ClippedShape c2 = new S2Polygon(S2Loop.full()).index.iterator().entry().clipped(0);
    assertFalse(S2ShapeUtil.equals(c1, c2));
  }

  public void testReferenceEmptyPolygon() {
    assertFalse(new S2Polygon().shape().getReferencePoint().contained());
  }

  public void testReferenceFullPolygon() {
    assertTrue(new S2Polygon(S2Loop.full()).shape().getReferencePoint().contained());
  }

  public void testReferenceDegenerateLoops() {
    List<S2Loop> loops = new ArrayList<>(Arrays.asList(
        makeLoop("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1"),
        makeLoop("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0"),
        makeLoop("5:5, 6:6")));
    assertFalse(new S2Polygon(loops).shape().getReferencePoint().contained());
  }

  public void testReferenceInvertedLoops() {
    S2Polygon poly = new S2Polygon();
    poly.initOriented(new ArrayList<>(Arrays.asList(
        makeLoop("1:2, 1:1, 2:2"),
        makeLoop("3:4, 3:3, 4:4"))));
    assertTrue(S2ShapeUtil.containsBruteForce(poly.shape(), S2.origin()));
  }

  public void testReferencePartiallyDegenerateLoops() {
    List<S2Point> loop = new ArrayList<>();
    for (int iter = 0; iter < 100; iter++) {
      // First we construct a long convoluted edge chain that follows the S2CellId Hilbert curve.
      // At some random point along the curve, we insert a small triangular loop.
      int numVertices = 100;
      S2CellId start = getRandomCellId(S2CellId.MAX_LEVEL - 1);
      S2CellId end = start.advanceWrap(numVertices);
      S2CellId loopCellId = start.advanceWrap(random(numVertices - 2) + 1);
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
      ReferencePoint ref = S2ShapeUtil.getReferencePoint(new S2Loop(loop));
      assertEquals(S2ShapeUtil.containsBruteForce(triangleLoop, ref), ref.contained());
    }
  }

  public void testShapeToShapeId() {
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(poly1);
    index.add(poly1);
    index.add(poly2);

    Multimap<S2Shape, Integer> actual = S2ShapeUtil.shapeToShapeId(index);

    assertEquals(3, actual.size());
    assertEquals(Lists.newArrayList(0, 1), actual.get(poly1));
    assertEquals(Lists.newArrayList(2), actual.get(poly2));
  }

  public void testLowerBound() {
    int[] data = {4, 4, 4, 5, 8};
    assertEquals("before equal range", 0, S2ShapeUtil.lowerBound(0, 5, i -> -1 > data[i]));
    assertEquals("starts equal range", 0, S2ShapeUtil.lowerBound(0, 5, i -> 4 > data[i]));
    assertEquals("after equal range", 3, S2ShapeUtil.lowerBound(0, 5, i -> 5 > data[i]));
    assertEquals("after end", 5, S2ShapeUtil.lowerBound(0, 5, i -> 10 > data[i]));
  }

  public void testUpperBound() {
    int[] data = {1, 3, 5, 7, 7};
    assertEquals("after equal range", 5, S2ShapeUtil.upperBound(0, 5, i -> 8 < data[i]));
    assertEquals("ends equal range", 5, S2ShapeUtil.upperBound(0, 5, i -> 7 < data[i]));
    assertEquals("before equal range", 2, S2ShapeUtil.upperBound(0, 5, i -> 3 < data[i]));
    assertEquals("before start", 0, S2ShapeUtil.upperBound(0, 5, i -> 0 < data[i]));
  }
}
