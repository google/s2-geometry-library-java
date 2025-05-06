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

import static com.google.common.geometry.S2TextFormat.makeLaxPolylineOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2LaxPolylineLayer}. */
@RunWith(JUnit4.class)
public final class S2LaxPolylineLayerTest extends GeometryTestCase {

  @Test
  public void testNoEdges() {
    checkS2LaxPolylineShape(ImmutableList.of(), "", new S2Builder.Builder());
  }

  /**
   * Even with undirected edges, S2LaxPolylineLayer prefers to reconstruct edges in their original
   * direction.
   */
  @Test
  public void testOneEdge() {
    checkS2LaxPolylineShapeUnchanged("3:4, 1:1");
    checkS2LaxPolylineShapeUnchanged("1:1, 3:4");
  }

  /** A line with backtracking is reconstructed unchanged. */
  @Test
  public void testStraightLineWithBacktracking() {
    checkS2LaxPolylineShapeUnchanged("0:0, 1:0, 2:0, 3:0, 2:0, 1:0, 2:0, 3:0, 4:0");
  }

  /**
   * Test that the "early walk termination" code (which is needed by S2LaxPolylineLayer in order to
   * implement idempotency) does not create two polylines when it is possible to assemble the edges
   * into one.
   */
  @Test
  public void testEarlyWalkTerminationWithEndLoop1() {
    // This example tests a code path where the early walk termination code should not be triggered
    // at all (but was at one point due to a bug).
    S2Builder.Builder options = new S2Builder.Builder();
    options.setSnapFunction(new S2BuilderSnapFunctions.IntLatLngSnapFunction(2));
    checkS2LaxPolylineShape(ImmutableList.of("0:0, 0:2, 0:1"), "0:0, 0:1, 0:2, 0:1", options);
  }

  /**
   * This tests a different code path where the walk is terminated early (yield a polyline with one
   * edge), and then the walk is "maximized" by appending a two-edge loop to the end.
   */
  @Test
  public void testEarlyWalkTerminationWithEndLoop2() {
    checkS2LaxPolylineShape(
        ImmutableList.of("0:0, 0:1", "0:2, 0:1", "0:1, 0:2"),
        "0:0, 0:1, 0:2, 0:1",
        new S2Builder.Builder());
  }

  @Test
  public void testSimpleLoop() {
    checkS2LaxPolylineShapeUnchanged("0:0, 0:5, 5:5, 5:0, 0:0");
  }

  /**
   * This polyline consists of many overlapping loops that keep returning to the same starting
   * vertex (2:2). This tests whether the implementation is able to assemble the polyline in the
   * original order.
   */
  @Test
  public void testManyLoops() {
    checkS2LaxPolylineShapeUnchanged(
        "0:0, 2:2, 2:4, 2:2, 2:4, 4:4, 4:2, 2:2, 4:4, 4:2, 2:2, 2:0, 2:2, "
            + "2:0, 4:0, 2:2, 4:2, 2:2, 0:2, 0:4, 2:2, 0:4, 0:2, 2:2, 0:4, 2:2, "
            + "0:2, 2:2, 0:0, 0:2, 2:2, 0:0");
  }

  /**
   * This test consists of 5 squares that touch diagonally, similar to the 5 white squares of a 3x3
   * chessboard. The edges of these squares need to be reordered to assemble them into a single
   * unbroken polyline.
   */
  @Test
  public void testUnorderedLoops() {
    checkS2LaxPolylineShape(
        ImmutableList.of(
            "3:3, 3:2, 2:2, 2:3, 3:3",
            "1:0, 0:0, 0:1, 1:1, 1:0",
            "3:1, 3:0, 2:0, 2:1, 3:1",
            "1:3, 1:2, 0:2, 0:1, 1:3",
            "1:1, 1:2, 2:2, 2:1, 1:1" // Central square
            ),
        "3:3, 3:2, 2:2, 2:1, 3:1, 3:0, 2:0, 2:1, 1:1, 1:0, 0:0, "
            + "0:1, 1:1, 1:2, 0:2, 0:1, 1:3, 1:2, 2:2, 2:3, 3:3",
        new S2Builder.Builder());
  }

  /**
   * Test reconstruction of a polyline where two edges have been split into many pieces by crossing
   * edges. This example is particularly challenging because (1) the edges form a loop, and (2) the
   * first and last edges are identical (but reversed). This is designed to test the heuristics that
   * attempt to find the first edge of the input polyline.
   */
  @Test
  public void testSplitEdges() {
    S2Builder.Builder options = new S2Builder.Builder();
    options.setSplitCrossingEdges(true);
    options.setSnapFunction(new S2BuilderSnapFunctions.IntLatLngSnapFunction(7));
    checkS2LaxPolylineShape(
        ImmutableList.of("0:10, 0:0, 1:0, -1:2, 1:4, -1:6, 1:8, -1:10, -5:0, 0:0, 0:10"),
        "0:10, 0:9, 0:7, 0:5, 0:3, 0:1, 0:0, 1:0, 0:1, -1:2, 0:3, 1:4, 0:5, "
            + "-1:6, 0:7, 1:8, 0:9, -1:10, -5:0, 0:0, 0:1, 0:3, 0:5, 0:7, 0:9, 0:10",
        options);
  }

  @Test
  public void testSimpleEdgeLabels() {
    S2Builder builder = new S2Builder.Builder().build();
    IntVector labelSetIds = new IntVector();
    IdSetLexicon labelSetLexicon = new IdSetLexicon();
    S2LaxPolylineLayer layer =
        new S2LaxPolylineLayer(
            new S2LaxPolylineLayer.Options(EdgeType.UNDIRECTED), labelSetLexicon, labelSetIds);
    builder.startLayer(layer);
    builder.setLabel(5);
    builder.addShape(makeLaxPolylineOrDie("0:0, 0:1, 0:2"));
    builder.pushLabel(7);
    builder.addShape(makeLaxPolylineOrDie("0:3, 0:2"));
    builder.clearLabels();
    builder.addShape(makeLaxPolylineOrDie("0:3, 0:4, 0:5"));
    builder.setLabel(11);
    builder.addShape(makeLaxPolylineOrDie("0:6, 0:5"));
    S2Error error = new S2Error();
    assertTrue(builder.build(error));

    ImmutableList<ImmutableList<Integer>> expected =
        ImmutableList.of(
            ImmutableList.of(5),
            ImmutableList.of(5),
            ImmutableList.of(5, 7),
            ImmutableList.of(),
            ImmutableList.of(),
            ImmutableList.of(11));
    assertEquals(expected.size(), labelSetIds.size());
    for (int i = 0; i < expected.size(); ++i) {
      assertEquals(expected.get(i).size(), labelSetLexicon.idSet(labelSetIds.get(i)).size());
      int j = 0;
      for (OfInt iter = labelSetLexicon.idSet(labelSetIds.get(i)).intIterator(); iter.hasNext(); ) {
        int label = iter.nextInt();
        assertEquals((int) expected.get(i).get(j++), label);
      }
    }
  }

  @Test
  public void testAntipodalVertices() {
    S2Builder builder = new S2Builder.Builder().build();
    S2LaxPolylineLayer layer = new S2LaxPolylineLayer();
    builder.startLayer(layer);
    builder.addEdge(new S2Point(1, 0, 0), new S2Point(-1, 0, 0));
    S2Error error = new S2Error();
    assertTrue(builder.build(error));
    S2LaxPolylineShape output = layer.getPolyline();
    assertEquals(2, output.numVertices());
    assertEquals(output.vertex(0), new S2Point(1, 0, 0));
    assertEquals(output.vertex(1), new S2Point(-1, 0, 0));
  }

  private void checkS2LaxPolylineShape(
      List<String> inputStrs, String expectedStr, EdgeType edgeType, S2Builder.Builder options) {
    S2Builder builder = options.build();
    S2LaxPolylineLayer layer = new S2LaxPolylineLayer(new S2LaxPolylineLayer.Options(edgeType));
    builder.startLayer(layer);
    for (String inputStr : inputStrs) {
      builder.addShape(makeLaxPolylineOrDie(inputStr));
    }
    S2Error error = new S2Error();
    assertTrue(builder.build(error));
    assertEquals(expectedStr, S2TextFormat.toString(layer.getPolyline()));
  }

  // Convenience function that tests both directed and undirected edges.
  private void checkS2LaxPolylineShape(
      List<String> inputStrs, String expectedStr, S2Builder.Builder options) {
    checkS2LaxPolylineShape(inputStrs, expectedStr, EdgeType.DIRECTED, options);
    checkS2LaxPolylineShape(inputStrs, expectedStr, EdgeType.UNDIRECTED, options);
  }

  private void checkS2LaxPolylineShapeUnchanged(String inputStr) {
    checkS2LaxPolylineShape(ImmutableList.of(inputStr), inputStr, new S2Builder.Builder());
  }
}
