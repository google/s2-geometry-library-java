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

import static com.google.common.geometry.S2TextFormat.makePolylineOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2BuilderSnapFunctions.IntLatLngSnapFunction;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2PolylineLayer}. */
@RunWith(JUnit4.class)
public final class S2PolylineLayerTest extends GeometryTestCase {

  @Test
  public void testNoEdges() {
    checkS2Polyline(ImmutableList.of(), "", null);
  }

  /**
   * Even with undirected edges, S2PolylineLayer prefers to reconstruct edges in their original
   * direction.
   */
  @Test
  public void testOneEdge() {
    checkS2PolylineUnchanged("3:4, 1:1");
    checkS2PolylineUnchanged("1:1, 3:4");
  }

  /**
   * Check that a polyline can be reconstructed with edges that backtrack between the same vertices,
   * for walks of increasing complexity.
   */
  @Test
  public void testStraightLineWithBacktracking() {
    // fwd, back, fwd
    checkS2PolylineUnchanged("0:0, 1:0, 0:0, 1:0");
    // fwd, fwd, back
    checkS2PolylineUnchanged("0:0, 1:0, 2:0, 1:0");
    // fwd, fwd, back, fwd
    checkS2PolylineUnchanged("0:0, 1:0, 2:0, 1:0, 2:0");
    // fwd, fwd, back, fwd, fwd
    checkS2PolylineUnchanged("0:0, 1:0, 2:0, 1:0, 2:0, 3:0");
    // fwd, fwd, back, back, fwd, fwd
    checkS2PolylineUnchanged("0:0, 1:0, 2:0, 1:0, 0:0, 1:0, 2:0");
    // A loop from the second vertex of the line: fwd, fwd, fwd, back, back, fwd, fwd
    checkS2PolylineUnchanged("0:0, 1:0, 2:0, 3:0, 2:0, 1:0, 2:0, 3:0");
    // A loop from the first vertex of the line: fwd, fwd, back, back, fwd, fwd, fwd
    checkS2PolylineUnchanged("0:0, 1:0, 2:0, 1:0, 0:0, 1:0, 2:0, 3:0");
    checkS2PolylineUnchanged("0:0, 1:0, 2:0, 3:0, 2:0, 1:0, 2:0, 3:0");
    checkS2PolylineUnchanged("0:0, 1:0, 2:0, 3:0, 2:0, 1:0, 2:0, 3:0, 4:0");
  }

  /**
   * Test that the "early walk termination" code (which is needed by S2PolylineVectorLayer in order
   * to implement idempotency) does not create two polylines when it is possible to assemble the
   * edges into one.
   *
   * <p>This example tests a code path where the early walk termination code should not be triggered
   * at all.
   */
  @Test
  public void testEarlyWalkTerminationWithEndLoop1() {
    S2Builder.Builder options = new S2Builder.Builder();
    options.setSnapFunction(new IntLatLngSnapFunction(2));
    checkS2Polyline(ImmutableList.of("0:0, 0:2, 0:1"), "0:0, 0:1, 0:2, 0:1", options);
  }

  /**
   * This tests a different code path where the walk is terminated early (yield a polyline with one
   * edge), and then the walk is "maximized" by appending a two-edge loop to the end.
   */
  @Test
  public void testEarlyWalkTerminationWithEndLoop2() {
    checkS2Polyline(
        ImmutableList.of("0:0, 0:1", "0:2, 0:1", "0:1, 0:2"), "0:0, 0:1, 0:2, 0:1", null);
  }

  @Test
  public void testSimpleLoop() {
    checkS2PolylineUnchanged("0:0, 0:5, 5:5, 5:0, 0:0");
  }

  /**
   * This polyline consists of many overlapping loops that keep returning to the same starting
   * vertex (2:2). This tests whether the implementation is able to assemble the polyline in the
   * original order.
   */
  @Test
  public void testManyLoops() {
    checkS2PolylineUnchanged(
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
    checkS2Polyline(
        ImmutableList.of(
            "3:3, 3:2, 2:2, 2:3, 3:3",
            "1:0, 0:0, 0:1, 1:1, 1:0",
            "3:1, 3:0, 2:0, 2:1, 3:1",
            "1:3, 1:2, 0:2, 0:1, 1:3",
            "1:1, 1:2, 2:2, 2:1, 1:1" // Central square
            ),
        "3:3, 3:2, 2:2, 2:1, 3:1, 3:0, 2:0, 2:1, 1:1, 1:0, 0:0, "
            + "0:1, 1:1, 1:2, 0:2, 0:1, 1:3, 1:2, 2:2, 2:3, 3:3",
        null);
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
    options.setSnapFunction(new IntLatLngSnapFunction(7));
    checkS2Polyline(
        ImmutableList.of("0:10, 0:0, 1:0, -1:2, 1:4, -1:6, 1:8, -1:10, -5:0, 0:0, 0:10"),
        "0:10, 0:9, 0:7, 0:5, 0:3, 0:1, 0:0, 1:0, 0:1, -1:2, 0:3, 1:4, 0:5, "
            + "-1:6, 0:7, 1:8, 0:9, -1:10, -5:0, 0:0, 0:1, 0:3, 0:5, 0:7, 0:9, 0:10",
        options);
  }

  /** Tests adding and retrieving labels from polyline edges. */
  @Test
  public void testSimpleEdgeLabels() {
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();

    IntVector labelSetIds = new IntVector();
    IdSetLexicon labelSetLexicon = new IdSetLexicon();
    S2PolylineLayer.Options layerOptions = new S2PolylineLayer.Options(EdgeType.UNDIRECTED);
    S2PolylineLayer layer = new S2PolylineLayer(layerOptions, labelSetLexicon, labelSetIds);
    builder.startLayer(layer);
    builder.setLabel(5);
    builder.addPolyline(makePolylineOrDie("0:0, 0:1, 0:2"));
    builder.pushLabel(7);
    builder.addPolyline(makePolylineOrDie("0:3, 0:2"));
    builder.clearLabels();
    builder.addPolyline(makePolylineOrDie("0:3, 0:4, 0:5"));
    builder.setLabel(11);
    builder.addPolyline(makePolylineOrDie("0:6, 0:5"));
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
  public void testInvalidPolyline() {
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();

    S2PolylineLayer.Options layerOptions = new S2PolylineLayer.Options();
    layerOptions.setValidate(true);
    S2PolylineLayer layer = new S2PolylineLayer(layerOptions);
    builder.startLayer(layer);
    List<S2Point> vertices = new ArrayList<>();
    vertices.add(new S2Point(1, 0, 0));
    vertices.add(new S2Point(-1, 0, 0));
    // Internal assertions must be disabled to allow invalid polylines to be built.
    S2Polyline input = uncheckedCreate(() -> new S2Polyline(vertices));
    builder.addPolyline(input);
    S2Error error = new S2Error();
    assertFalse(uncheckedCreate(() -> builder.build(error)));
    assertEquals(S2Error.Code.ANTIPODAL_VERTICES, error.code());
  }

  private void checkS2Polyline(
      List<String> inputStrs, String expectedStr, EdgeType edgeType, S2Builder.Builder options) {
    if (options == null) {
      options = new S2Builder.Builder();
    }
    S2Builder builder = options.build();
    S2PolylineLayer.Options layerOptions = new S2PolylineLayer.Options(edgeType);
    S2PolylineLayer layer = new S2PolylineLayer(layerOptions);
    builder.startLayer(layer);
    for (String inputStr : inputStrs) {
      builder.addPolyline(makePolylineOrDie(inputStr));
    }
    S2Error error = new S2Error();
    assertTrue(builder.build(error));
    assertEquals(expectedStr, S2TextFormat.toString(layer.getPolyline()));
  }

  /** Convenience function that tests both directed and undirected edges. */
  private void checkS2Polyline(
      List<String> inputStrs, String expectedStr, S2Builder.Builder options) {
    if (options == null) {
      options = new S2Builder.Builder();
    }
    checkS2Polyline(inputStrs, expectedStr, EdgeType.DIRECTED, options);
    checkS2Polyline(inputStrs, expectedStr, EdgeType.UNDIRECTED, options);
  }

  /**
   * Checks that building the given S2TextFormat polyline into an S2PolylineLayer results in an
   * unchanged output polyline, for both directed and undirected edges.
   */
  private void checkS2PolylineUnchanged(String inputStr) {
    checkS2Polyline(ImmutableList.of(inputStr), inputStr, null);
  }
}
