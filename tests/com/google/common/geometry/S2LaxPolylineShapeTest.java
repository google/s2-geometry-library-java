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

import com.google.common.annotations.GwtIncompatible;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.Lists;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

@GwtIncompatible("S2ShapeAspect incompatible")
public class S2LaxPolylineShapeTest extends GeometryTestCase {
  private final S2Point a = makePoint("0:0");
  private final S2Point b = makePoint("0:1");
  private final S2Point c = makePoint("1:1");
  private final ImmutableMap<S2Point, String> labels = ImmutableMap.of(a, "a", b, "b", c, "c");
  private List<S2Point> vertices = Arrays.asList(a, b, c);
  
  public void testNoVertices() {
    S2LaxPolylineShape shape = S2LaxPolylineShape.create(ImmutableList.of());
    assertTrue(shape.isEmpty());
    assertFalse(shape.isFull());
    assertEquals("[[]]", edges(shape));
  }
  
  public void testOneVertex() {
    S2LaxPolylineShape line = S2LaxPolylineShape.create(ImmutableList.of(S2Point.X_POS));
    // Note that the C++ lax polyline stores a vertex, but there is no way to access it.
    // Whereas we discard it during construction. This would be an invisible difference, except we
    // also have a 'numVertices' method. It's package-private, but we should still assert this
    // behavior.
    assertEquals(0, line.numVertices());
    assertEquals(0, line.numEdges());
    assertEquals(1, line.numChains());
    assertTrue(line.isEmpty());
    assertFalse(line.isFull());
  }
  
  public void testTwoEdges() {
    S2LaxPolylineShape shape = S2LaxPolylineShape.create(vertices);
    assertEquals(1, shape.dimension());
    assertFalse(shape.isEmpty());
    assertFalse(shape.isFull());
    assertEquals("[[[a, b], [b, c]]]", edges(shape));
  }
  
  public void testDegenerateSingle() {
    assertEquals("[[[a, b], [b, a]]]", edges(S2LaxPolylineShape.createMulti(Arrays.asList(
        Arrays.asList(c),
        Arrays.asList(a, b, a)))));
  }
  
  public void testDegenerateMulti() {
    // Note that the C++ lax polyline only supports 1 chain today, so cannot encode this case yet.
    assertEquals("[[[b, a]], [[a, b], [b, c]], [[c, b]]]",
        edges(S2LaxPolylineShape.createMulti(Arrays.asList(
            Arrays.asList(b, a),
            Arrays.asList(a, b, c),
            Arrays.asList(b),
            Arrays.asList(c, b)))));
  }

  public void testEncodeDecode() throws IOException {
    checkEncodeDecode(S2LaxPolylineShape.create(makePolyline("0:0, 0:1, 1:1").vertices()));
    checkEncodeDecode(S2LaxPolylineShape.EMPTY);
  }

  private String edges(S2ShapeAspect.Mixed shape) {
    assertEquals(1, shape.dimension());
    MutableEdge edge = new MutableEdge();
    List<String[][]> chains = new ArrayList<>();
    int lastEdge = 0;
    for (int chainId = 0; chainId < shape.numChains(); chainId++) {
      List<String[]> chain = new ArrayList<>();
      for (int edgeOffset = 0; edgeOffset < shape.getChainLength(chainId); edgeOffset++) {
        shape.getChainEdge(chainId, edgeOffset, edge);
        S2Point a = edge.a;
        S2Point b = edge.b;
        assertEquals(a, shape.getChainVertex(chainId, edgeOffset));
        assertEquals(a, shape.chains().get(chainId).get(edgeOffset));
        chain.add(new String[] {labels.get(a), labels.get(b)});
        int edgeIndex = shape.getChainStart(chainId) + edgeOffset;
        assertEquals(lastEdge++, edgeIndex);
        shape.getEdge(edgeIndex, edge);
        assertEquals(a, edge.a);
        assertEquals(b, edge.b);
        assertEquals(a, shape.vertices().get(shape.vertexId(chainId, edgeIndex)));
      }
      chains.add(chain.toArray(new String[0][]));
    }
    assertEquals(lastEdge, shape.numEdges());
    String[][][] actual = chains.toArray(new String[0][][]);
    return Arrays.deepToString(actual);
  }

  public void checkEncodeDecode(S2LaxPolylineShape expected) throws IOException {
    for (S2Coder<S2LaxPolylineShape> coder :
        Lists.newArrayList(S2LaxPolylineShape.Coder.FAST, S2LaxPolylineShape.Coder.COMPACT)) {
      ByteArrayOutputStream out = new ByteArrayOutputStream();
      coder.encode(expected, out);
      Bytes data = Bytes.fromByteArray(out.toByteArray());
      Cursor cursor = data.cursor();
      S2LaxPolylineShape actual = coder.decode(data, cursor);
      assertTrue(S2ShapeUtil.equals(expected, actual));
      assertEquals(out.size(), cursor.position);
    }
  }
}
