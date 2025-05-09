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

import static com.google.common.geometry.S2CellId.FACE_CELLS;
import static com.google.common.geometry.S2TextFormat.makeLoop;
import static com.google.common.geometry.S2TextFormat.makePolygon;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.base.Splitter;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2LaxPolygonShape}. */
@RunWith(JUnit4.class)
public class S2LaxPolygonShapeTest extends GeometryTestCase {
  private static final S2Polygon POLYGON = makePolygon("0:0, 0:1, 1:0");
  private static final S2Polygon SNAPPED_POLYGON = new S2Polygon(new S2Loop(ImmutableList.of(
      FACE_CELLS[0].toPoint(),
      FACE_CELLS[1].toPoint(),
      FACE_CELLS[2].toPoint())));

  @Test
  public void testEmptyPolygon() {
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(new S2Polygon());
    assertEquals(0, shape.numChains());
    assertEquals(0, shape.numEdges());
    assertEquals(2, shape.dimension());
    assertTrue(shape.isEmpty());
    assertFalse(shape.isFull());
    assertFalse(shape.getReferencePoint().contained());
  }

  @Test
  public void testFullPolygon() {
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(new S2Polygon(S2Loop.full()));
    assertEquals(1, shape.numChains());
    assertEquals(0, shape.numEdges());
    assertEquals(2, shape.dimension());
    assertFalse(shape.isEmpty());
    assertTrue(shape.isFull());
    assertTrue(shape.getReferencePoint().contained());
  }

  @Test
  public void testSingleVertexPolygon() {
    // S2Polygon doesn't support single-vertex loops, so construct an S2LaxPolygonShape directly.
    List<List<S2Point>> loops = Arrays.asList(parsePoints("0:0"));
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(loops);
    assertEquals(1, shape.numChains());
    assertEquals(1, shape.numEdges());
    assertEquals(0, shape.getChainStart(0));
    assertEquals(1, shape.getChainLength(0));
    S2Point v = loops.get(0).get(0);
    assertEquals(v, edge(shape, 0).getStart());
    assertEquals(v, edge(shape, 0).getEnd());
    assertEquals(v, chainEdge(shape, 0, 0).getStart());
    assertEquals(v, chainEdge(shape, 0, 0).getEnd());
    assertEquals(2, shape.dimension());
    assertFalse(shape.isEmpty());
    assertFalse(shape.isFull());
    assertFalse(shape.getReferencePoint().contained());
  }

  @Test
  public void testSingleLoopPolygon() {
    // Test S2Polygon constructor.
    List<S2Point> vertices = parsePoints("0:0, 0:1, 1:1, 1:0");
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(new S2Polygon(new S2Loop(vertices)));
    assertEquals(1, shape.numChains());
    assertEquals(vertices.size(), shape.numEdges());
    assertEquals(1, shape.numChains());
    assertEquals(0, shape.getChainStart(0));
    assertEquals(vertices.size(), shape.getChainLength(0));
    for (int i = 0; i < vertices.size(); i++) {
      S2Point a = vertices.get(i);
      S2Point b = vertices.get((i + 1) % vertices.size());
      assertEquals(a, shape.getChainVertex(0, i));
      assertEquals(a, edge(shape, i).getStart());
      assertEquals(a, chainEdge(shape, 0, i).getStart());
      assertEquals(b, edge(shape, i).getEnd());
      assertEquals(b, chainEdge(shape, 0, i).getEnd());
    }
    assertEquals(2, shape.dimension());
    assertFalse(shape.isEmpty());
    assertFalse(shape.isFull());
    assertFalse(S2ShapeUtil.containsBruteForce(shape, S2.origin()));
  }

  @Test
  public void testMultiLoopPolygon() {
    // Test List<List<S2Point>> constructor.  Make sure that the loops are oriented so that the
    // interior of the polygon is always on the left.
    List<List<S2Point>> loops = Arrays.asList(
        parsePoints("0:0, 0:3, 3:3"), // CCW
        parsePoints("1:1, 2:2, 1:2")); // CW
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(loops);

    assertEquals(loops.size(), shape.numChains());
    int numVertices = 0;
    for (int i = 0; i < loops.size(); i++) {
      List<S2Point> loop = loops.get(i);
      assertEquals(loop.size(), shape.getChainLength(i));
      assertEquals(numVertices, shape.getChainStart(i));
      assertEquals(loop.size(), shape.getChainLength(i));
      for (int j = 0; j < loop.size(); j++) {
        assertEquals(loop.get(j), shape.getChainVertex(i, j));
        int k = numVertices + j;
        assertEquals(loop.get(j), edge(shape, k).getStart());
        assertEquals(loop.get((j + 1) % loop.size()), edge(shape, k).getEnd());
      }
      numVertices += loop.size();
    }
    assertEquals(numVertices, shape.numEdges());
    assertEquals(2, shape.dimension());
    assertFalse(shape.isEmpty());
    assertFalse(shape.isFull());
    assertFalse(S2ShapeUtil.containsBruteForce(shape, S2.origin()));
  }

  @Test
  public void testMultiLoopS2Polygon() {
    // Verify that the orientation of loops representing holes is reversed when converting from an
    // S2Polygon to an S2LaxPolygonShape.
    S2Polygon polygon = makePolygon("0:0, 0:3, 3:3; 1:1, 1:2, 2:2");
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(polygon);
    for (int i = 0; i < polygon.numLoops(); i++) {
      S2Loop loop = polygon.loop(i);
      for (int j = 0; j < loop.numVertices(); j++) {
        assertEquals(loop.orientedVertex(j), shape.getChainVertex(i, j));
      }
    }
  }

  @Test
  public void testManyLoopPolygon() {
    // Test a polygon with enough loops so that cumulative_vertices_ is used.
    List<List<S2Point>> loops = new ArrayList<>();
    for (int i = 0; i < 100; i++) {
      loops.add(S2Loop.makeRegularVertices(
          S2LatLng.fromDegrees(0, i).toPoint(), S1Angle.degrees(0.1), data.random(3)));
    }
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(loops);

    assertEquals(loops.size(), shape.numChains());
    int numVertices = 0;
    for (int i = 0; i < loops.size(); i++) {
      List<S2Point> loop = loops.get(i);
      assertEquals(loop.size(), shape.getChainLength(i));
      assertEquals(numVertices, shape.getChainStart(i));
      assertEquals(loop.size(), shape.getChainLength(i));
      for (int j = 0; j < loop.size(); j++) {
        assertEquals(loop.get(j), shape.getChainVertex(i, j));
        int k = numVertices + j;
        assertEquals(loop.get(j), edge(shape, k).getStart());
        assertEquals(loop.get((j + 1) % loop.size()), edge(shape, k).getEnd());
      }
      numVertices += loop.size();
    }
    assertEquals(numVertices, shape.numEdges());
  }

  @Test
  public void testDegenerateLoops() {
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(Arrays.asList(
        parsePoints("1:1, 1:2, 2:2, 1:2, 1:3, 1:2, 1:1"),
        parsePoints("0:0, 0:3, 0:6, 0:9, 0:6, 0:3, 0:0"),
        parsePoints("5:5, 6:6")));
    assertFalse(shape.getReferencePoint().contained());
  }

  @Test
  public void testInvertedLoops() {
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(Arrays.asList(
        parsePoints("1:2, 1:1, 2:2"),
        parsePoints("3:4, 3:3, 4:4")));
    assertTrue(S2ShapeUtil.containsBruteForce(shape, S2.origin()));
  }

  @Test
  public void testCompareToS2Loop() {
    for (int iter = 0; iter < 100; iter++) {
      S2FractalBuilder fractal = new S2FractalBuilder(data.rand);
      fractal.setMaxLevel(data.random(5));
      fractal.setFractalDimension(data.uniform(1, 2));
      S2Point center = data.getRandomPoint();
      S2Loop loop = fractal.makeLoop(data.getRandomFrameAt(center), S1Angle.degrees(5));

      // Compare S2Loop to S2LaxPolygonShape.
      compareS2LoopToShape(loop, S2LaxPolygonShape.create(Arrays.asList(loop.vertices())));
    }
  }

  @Test
  public void testEncodeDecode() throws IOException {
    checkEncodeDecode(S2LaxPolygonShape.EMPTY);
    checkEncodeDecode(S2LaxPolygonShape.FULL);

    // 1 loop.
    S2Loop shell = makeLoop("0:0, 0:4, 4:4, 4:0");
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(new S2Polygon(shell));
    checkEncodeDecode(shape);

    // 2 loops.
    S2Loop hole = makeLoop("1:1, 1:3, 3:3, 3:1");
    shape = S2LaxPolygonShape.create(new S2Polygon(Lists.newArrayList(shell, hole)));
    checkEncodeDecode(shape);
  }

  @Test
  public void testCoderFast() {
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(POLYGON);
    assertShapesEqual(shape, roundtrip(S2LaxPolygonShape.FAST_CODER, shape));
  }

  @Test
  public void testCoderCompact() {
    S2LaxPolygonShape shape = S2LaxPolygonShape.create(SNAPPED_POLYGON);
    assertShapesEqual(shape, roundtrip(S2LaxPolygonShape.COMPACT_CODER, shape));
  }

  private void compareS2LoopToShape(S2Loop loop, S2Shape shape) {
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(shape);
    S2Cap cap = loop.getCapBound();
    S2ContainsPointQuery query = new S2ContainsPointQuery(index);
    for (int iter = 0; iter < 100; iter++) {
      S2Point point = data.samplePoint(cap);
      assertEquals(loop.contains(point), query.shapeContains(0, point));
    }
  }

  private static List<S2Point> parsePoints(String str) {
    return Lists.transform(
        Splitter.on(",").trimResults().splitToList(str), S2TextFormat::makePoint);
  }

  private static MutableEdge edge(S2Shape shape, int edgeId) {
    MutableEdge edge = new MutableEdge();
    shape.getEdge(edgeId, edge);
    return edge;
  }

  private static MutableEdge chainEdge(S2Shape shape, int chainId, int edgeId) {
    MutableEdge edge = new MutableEdge();
    shape.getChainEdge(chainId, edgeId, edge);
    return edge;
  }

  private void checkEncodeDecode(S2LaxPolygonShape expected) throws IOException {
    for (S2Coder<S2LaxPolygonShape> coder :
        Lists.newArrayList(S2LaxPolygonShape.FAST_CODER, S2LaxPolygonShape.COMPACT_CODER)) {
      ByteArrayOutputStream out = new ByteArrayOutputStream();
      coder.encode(expected, out);
      Bytes data = Bytes.fromByteArray(out.toByteArray());
      Cursor cursor = data.cursor();
      S2LaxPolygonShape actual = coder.decode(data, cursor);

      // Test the edges in random order.
      int numEdges = expected.numEdges();
      List<Integer> edgeIds = new ArrayList<>(numEdges);
      for (int i = 0; i < numEdges; i++) {
        edgeIds.add(i);
      }
      Collections.shuffle(edgeIds, rand());

      S2Shape.ChainPosition expectedChainPosition = new S2Shape.ChainPosition();
      S2Shape.ChainPosition actualChainPosition = new S2Shape.ChainPosition();
      edgeIds.forEach(edgeId -> {
        expected.getChainPosition(edgeId, expectedChainPosition);
        actual.getChainPosition(edgeId, actualChainPosition);
        assertTrue(expectedChainPosition.isEqualTo(actualChainPosition));
      });

      assertTrue(S2ShapeUtil.equals(expected, actual));
      assertEquals(out.size(), cursor.position);
    }
  }
}
