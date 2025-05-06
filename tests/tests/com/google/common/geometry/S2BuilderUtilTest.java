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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertThrows;

import com.google.common.geometry.S2BuilderGraph.EdgeList;
import com.google.common.geometry.S2BuilderUtil.GraphShape;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.util.Comparator;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/**
 * Tests for S2BuilderUtil.GraphShape and static methods in S2BuilderUtil.
 *
 * <p>S2BuilderUtil.GraphClone and S2BuilderUtil.GraphCloningLayer are tested as part of the
 * S2Builder tests.
 */
@RunWith(JUnit4.class)
public final class S2BuilderUtilTest extends GeometryTestCase {

  private final List<S2Point> vertices = S2TextFormat.parsePointsOrDie("0:0, 2:0, 2:2, 0:2");

  /** Sorts S2Points by their distance from lat,lng 0:0 */
  private static class FromOrigin implements Comparator<S2Point> {
    private final S2Point origin = S2TextFormat.makePointOrDie("0:0");
    @Override
    public int compare(S2Point left, S2Point right) {
      return S2Predicates.compareDistances(origin, left, right);
    }
  }

  @Test
  public void testDeduplicateSortedList() {
    List<S2Point> input = S2TextFormat.parsePointsOrDie(
        "2:0, 2:0, 0:0, 1:0, 1:0, 2:2, 1:0, 1:2, 2:2, 2:2");
    List<S2Point> expected = S2TextFormat.parsePointsOrDie(
        "0:0, 1:0, 2:0, 1:2, 2:2");
    input.sort(new FromOrigin());
    S2BuilderUtil.deduplicateSortedList(input);
    assertEquals(expected, input);
  }

  @Test
  public void testGraphShapeBasic() {
    // Set up a GraphShape on the 'vertices' above that contains two edges.
    EdgeList edges = new EdgeList();
    edges.add(0, 1); // 0:0 to 2:0
    edges.add(2, 3); // 2:2 to 0:2
    final GraphShape vShape = new GraphShape(edges, vertices);

    assertEquals(1, vShape.dimension());
    assertEquals(2, vShape.numEdges());
    assertFalse(vShape.hasInterior());
    assertFalse(vShape.containsOrigin());
    assertEquals(2, vShape.numChains());

    MutableEdge edge = new MutableEdge();
    vShape.getEdge(0, edge);
    assertSame(vertices.get(0), edge.getStart());
    assertSame(vertices.get(1), edge.getEnd());
    vShape.getEdge(1, edge);
    assertSame(vertices.get(2), edge.getStart());
    assertSame(vertices.get(3), edge.getEnd());

    // There are two chains, each one edge.
    assertEquals(0, vShape.getChainStart(0));
    assertEquals(1, vShape.getChainStart(1));

    assertEquals(1, vShape.getChainLength(0));
    assertEquals(1, vShape.getChainLength(1));

    MutableEdge mutableEdge = new MutableEdge();
    vShape.getChainEdge(0, 0, mutableEdge);
    assertSame(vertices.get(0), mutableEdge.getStart());
    assertSame(vertices.get(1), mutableEdge.getEnd());
    vShape.getChainEdge(1, 0, mutableEdge);
    assertSame(vertices.get(2), mutableEdge.getStart());
    assertSame(vertices.get(3), mutableEdge.getEnd());

    // getChainVertex can get the start and end of each of the two edges.
    assertSame(vertices.get(0), vShape.getChainVertex(0, 0));
    assertSame(vertices.get(1), vShape.getChainVertex(0, 1));
    assertSame(vertices.get(2), vShape.getChainVertex(1, 0));
    assertSame(vertices.get(3), vShape.getChainVertex(1, 1));
  }

  @Test
  public void testGetChainArgumentValidation() {
    EdgeList edges = new EdgeList();
    edges.add(0, 1); // 0:0 to 2:0
    edges.add(2, 3); // 2:2 to 0:2
    final GraphShape vShape = new GraphShape(edges, vertices);

    // VEdge has two chains of length 1. Any attept to go beyond that must throw.
    assertThrows(
        "Index out of range", IndexOutOfBoundsException.class, () -> vShape.getChainStart(2));

    assertThrows(
        "Index out of range", IndexOutOfBoundsException.class, () -> vShape.getChainLength(2));

    MutableEdge mutableEdge = new MutableEdge();
    assertThrows(
        "Index out of range",
        IndexOutOfBoundsException.class,
        () -> vShape.getChainEdge(2, 0, mutableEdge));
    // Chains have length 0 so edge offset 1 must also fail.
    assertThrows(
        "Index out of range",
        IndexOutOfBoundsException.class,
        () -> vShape.getChainEdge(0, 1, mutableEdge));

    // Offsets 0 and 1 are valid vertex offsets, anything higher is invalid.
    assertThrows(
        "Index out of range",
        IndexOutOfBoundsException.class,
        () -> assertEquals(vertices.get(0), vShape.getChainVertex(2, 0)));
    assertThrows(
        "Index out of range",
        IndexOutOfBoundsException.class,
        () -> assertEquals(vertices.get(0), vShape.getChainVertex(0, 2)));
  }
}
