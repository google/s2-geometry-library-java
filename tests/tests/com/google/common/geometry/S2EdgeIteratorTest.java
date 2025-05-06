/*
 * Copyright 2025 Google Inc.
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

import static com.google.common.geometry.S2TextFormat.makeIndexOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.geometry.S2Shape.MutableEdge;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for S2EdgeIterator. */
@RunWith(JUnit4.class)
public final class S2EdgeIteratorTest {

  /** Returns the full list of edges in the given S2ShapeIndex. */
  private List<MutableEdge> getEdges(S2ShapeIndex index) {
    List<MutableEdge> result = new ArrayList<>();
    for (S2Shape shape : index.getShapes()) {
      if (shape == null) {
        continue;
      }
      for (int j = 0; j < shape.numEdges(); ++j) {
        MutableEdge edge = new MutableEdge();
        shape.getEdge(j, edge);
        result.add(edge);
      }
    }
    return result;
  }

  /** Verifies that the edges produced by an S2EdgeIterator matches getEdges. */
  private void verify(S2ShapeIndex index) {
    List<MutableEdge> expected = getEdges(index);

    int i = 0;
    int shapeId = -1;
    int edgeId = -1;
    MutableEdge actualEdge = new MutableEdge();
    for (S2EdgeIterator it = new S2EdgeIterator(index); !it.done(); it.next(), ++edgeId, ++i) {
      // The iterator visits the edges of each shape in order.  When we see a new shape id, reset
      // the edgeId count.
      if (it.shapeId() != shapeId) {
        shapeId = it.shapeId();
        edgeId = 0;
      }

      assertTrue(i < expected.size());
      it.getEdge(actualEdge);
      assertTrue(expected.get(i).isEqualTo(actualEdge));
      assertEquals(it.edgeId(), edgeId);
      // assertEquals(it.shape_edgeId(), ShapeEdgeId(shapeId, edgeId));
    }
  }

  @Test
  public void testEmpty() {
    S2ShapeIndex index = makeIndexOrDie("##");
    verify(index);
  }

  @Test
  public void testPoints() {
    S2ShapeIndex index = makeIndexOrDie("0:0|1:1##");
    verify(index);
  }

  @Test
  public void testLines() {
    S2ShapeIndex index = makeIndexOrDie("#0:0,10:10|5:5,5:10|1:2,2:1#");
    verify(index);
  }

  @Test
  public void testPolygons() {
    S2ShapeIndex index = makeIndexOrDie("##10:10,10:0,0:0|-10:-10,-10:0,0:0,0:-10");
    verify(index);
  }

  @Test
  public void testCollection() {
    S2ShapeIndex index =
        makeIndexOrDie(
            "1:1|7:2#1:1,2:2,3:3|2:2,1:7#10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");
    verify(index);
  }

  @Test
  public void testAssignmentAndEquality() {
    S2ShapeIndex index1 =
        makeIndexOrDie(
            "1:1|7:2#1:1,2:2,3:3|2:2,1:7#10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");

    S2ShapeIndex index2 =
        makeIndexOrDie(
            "1:1|7:2#1:1,2:2,3:3|2:2,1:7#10:10,10:0,0:0;20:20,20:10,10:10|15:15,15:0,0:0");

    S2EdgeIterator it1 = new S2EdgeIterator(index1);
    S2EdgeIterator it2 = new S2EdgeIterator(index2);

    // Different indices, so different iterators.
    assertFalse(it1.isEqualTo(it2));

    // Make it2 a copy of it1, on the same index and at the same position.
    it1 = new S2EdgeIterator(it2);
    assertTrue(it1.isEqualTo(it2));

    // Different positions, so not equal.
    it1.next();
    assertFalse(it1.isEqualTo(it2));

    // Same positions again.
    it2.next();
    assertTrue(it1.isEqualTo(it2));
  }
}
