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

import static com.google.common.geometry.S2TextFormat.makeLoopOrDie;
import static com.google.common.geometry.S2TextFormat.makePolylineOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2PolygonLayer}. */
@RunWith(JUnit4.class)
public final class S2PolygonLayerTest extends GeometryTestCase {

  @Test
  public void testEmpty() {
    checkS2PolygonUnchanged("");
  }

  @Test
  public void testFull() {
    checkS2PolygonUnchanged("full");
  }

  @Test
  public void testSmallLoop() {
    checkS2PolygonUnchanged("0:0, 0:1, 1:1");
  }

  @Test
  public void testThreeLoops() {
    // The second two loops are nested.
    checkS2PolygonUnchanged("0:1, 1:1, 0:0; " + "3:3, 3:6, 6:6, 6:3; " + "4:4, 4:5, 5:5, 5:4");
  }

  @Test
  public void testPartialLoop() {
    checkS2PolygonError(
        ImmutableList.of("0:1, 2:3, 4:5"), S2Error.Code.BUILDER_EDGES_DO_NOT_FORM_LOOPS);
  }

  @Test
  public void testInvalidPolygon() {
    // TODO(user): When S2Polygon reworking is done, only one of these errors will occur.
    checkS2PolygonError(
        ImmutableList.of("0:0, 0:10, 10:0, 10:10, 0:0"),
        ImmutableSet.of(S2Error.Code.LOOP_SELF_INTERSECTION, S2Error.Code.OVERLAPPING_GEOMETRY));
  }

  /**
   * Check that S2PolygonLayer can assemble polygons even when there are duplicate edges (after
   * sibling pairs are removed), and then report the duplicate edges as an error.
   */
  @Test
  public void testDuplicateInputEdges() {
    S2Builder.Builder builderOptions = new S2Builder.Builder();
    S2Builder builder = builderOptions.build();

    S2PolygonLayer.Options layerOptions = new S2PolygonLayer.Options();
    layerOptions.setValidate(true);
    S2PolygonLayer layer = new S2PolygonLayer(layerOptions);
    builder.startLayer(layer);

    builder.addPolyline(makePolylineOrDie("0:0, 0:2, 2:2, 1:1, 0:2, 2:2, 2:0, 0:0"));
    S2Error error = new S2Error();
    // Internal assertions must be disabled to allow invalid polygons to be built.
    boolean result = uncheckedCreate(() -> builder.build(error));
    assertFalse(result);

    // Different implementations of S2Polygon validation find one or the other. Both are reasonable.
    // TODO(user): When S2Polygon reworking is done, only one of these errors will occur.
    ImmutableSet<S2Error.Code> expectedCodes =
        ImmutableSet.of(
            S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS,
            S2Error.Code.POLYGON_LOOPS_SHARE_EDGE);

    assertTrue("Got " + error.code(), expectedCodes.contains(error.code()));
    S2Polygon output = layer.getPolygon();
    assertEquals(2, output.numLoops());
    S2Loop loop0 = makeLoopOrDie("0:0, 0:2, 2:2, 2:0");
    S2Loop loop1 = makeLoopOrDie("0:2, 2:2, 1:1");
    assertEquals(output.loop(0), loop0);
    assertEquals(output.loop(1), loop1);
  }

  /**
   * Tests fetching edge labels on directed edges. Uses a polygon consisting of 3 loops. The loops
   * are reordered and some of the loops are inverted during S2Polygon construction.
   */
  @Test
  public void testDirectedEdgeLabels() {
    checkEdgeLabels(EdgeType.DIRECTED);
  }

  /**
   * Tests fetching edge labels on undirected edges. Uses a polygon consisting of 3 loops. The loops
   * are reordered and some of the loops are inverted during S2Polygon construction.
   */
  @Test
  public void testUndirectedEdgeLabels() {
    checkEdgeLabels(EdgeType.UNDIRECTED);
  }

  /** Tests the situation where labels are requested but none were provided. */
  @Test
  public void testLabelsRequestedButNotProvided() {
    S2Builder.Builder builderOptions = new S2Builder.Builder();
    S2Builder builder = builderOptions.build();
    IntVector labelSetIds = new IntVector();
    IdSetLexicon labelSetLexicon = new IdSetLexicon();
    S2PolygonLayer layer =
        new S2PolygonLayer(new S2PolygonLayer.Options(), labelSetLexicon, labelSetIds);
    builder.startLayer(layer);
    builder.addPolyline(makePolylineOrDie("0:0, 0:1, 1:0, 0:0"));
    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));

    ArrayList<IntVector> labelSetIdsForLoops = layer.getLabelSetIdsForLoops();
    assertEquals(1, labelSetIdsForLoops.size()); // One loop.
    assertEquals(3, labelSetIdsForLoops.get(0).size()); // Three edges.
    labelSetIdsForLoops
        .get(0)
        .forEach(labelSetId -> assertEquals(IdSetLexicon.EMPTY_SET_ID, labelSetId));

    S2Polygon output = layer.getPolygon();
    assertEquals(1, output.numLoops());
    assertEquals(3, labelSetIds.size()); // Three edges.
    labelSetIds.forEach(labelSetId -> assertEquals(IdSetLexicon.EMPTY_SET_ID, labelSetId));
  }

  /** Three loops (two shells and one hole) that combine into one. */
  @Test
  public void testThreeLoopsIntoOne() {
    checkS2Polygon(
        ImmutableList.of(
            "10:0, 0:0, 0:10, 5:10, 10:10, 10:5",
            "0:10, 0:15, 5:15, 5:10",
            "10:10, 5:10, 5:5, 10:5"),
        "10:5, 10:0, 0:0, 0:10, 0:15, 5:15, 5:10, 5:5");
  }

  /**
   * A big CCW triangle containing 3 CW triangular holes. The whole thing looks like a pyramid of
   * nine triangles. The output consists of 6 positive triangles with no holes.
   */
  @Test
  public void testTrianglePyramid() {
    checkS2Polygon(
        ImmutableList.of(
            "0:0, 0:2, 0:4, 0:6, 1:5, 2:4, 3:3, 2:2, 1:1",
            "0:2, 1:1, 1:3",
            "0:4, 1:3, 1:5",
            "1:3, 2:2, 2:4"),
        "0:4, 0:6, 1:5; 2:4, 3:3, 2:2; 2:2, 1:1, 1:3; "
            + "1:1, 0:0, 0:2; 1:3, 0:2, 0:4; 1:3, 1:5, 2:4");
  }

  /**
   * A complex set of nested polygons, with the loops in random order and the vertices in random
   * cyclic order within each loop. This test checks that the order (after S2Polygon.initNested is
   * called) is preserved exactly, whether directed or undirected edges are used.
   */
  @Test
  public void testComplexNesting() {
    checkS2PolygonUnchanged(
        "47:15, 47:5, 5:5, 5:15; "
            + "35:12, 35:7, 27:7, 27:12; "
            + "1:50, 50:50, 50:1, 1:1; "
            + "42:22, 10:22, 10:25, 42:25; "
            + "47:30, 47:17, 5:17, 5:30; "
            + "7:27, 45:27, 45:20, 7:20; "
            + "37:7, 37:12, 45:12, 45:7; "
            + "47:47, 47:32, 5:32, 5:47; "
            + "50:60, 50:55, 1:55, 1:60; "
            + "25:7, 17:7, 17:12, 25:12; "
            + "7:7, 7:12, 15:12, 15:7");
  }

  /** Five nested loops that touch at one common point. */
  @Test
  public void testFiveLoopsTouchingAtOneCommonPoint() {
    checkS2PolygonUnchanged(
        "0:0, 0:10, 10:10, 10:0; "
            + "0:0, 1:9, 9:9, 9:1; "
            + "0:0, 2:8, 8:8, 8:2; "
            + "0:0, 3:7, 7:7, 7:3; "
            + "0:0, 4:6, 6:6, 6:4");
  }

  /**
   * Four diamonds nested inside each other, where each diamond shares two vertices with the diamond
   * inside it and shares its other two vertices with the diamond that contains it. The resulting
   * shape looks vaguely like an eye made out of chevrons.
   */
  @Test
  public void testFourNestedDiamondsTouchingAtTwoPointsPerPair() {
    checkS2Polygon(
        ImmutableList.of(
            "0:10, -10:0, 0:-10, 10:0",
            "0:-20, -10:0, 0:20, 10:0",
            "0:-10, -5:0, 0:10, 5:0",
            "0:5, -5:0, 0:-5, 5:0"),
        "10:0, 0:10, -10:0, 0:20; "
            + "0:-20, -10:0, 0:-10, 10:0; "
            + "5:0, 0:-10, -5:0, 0:-5; "
            + "0:5, -5:0, 0:10, 5:0");
  }

  /** Seven diamonds nested within each other touching at one point between each nested pair. */
  @Test
  public void testSevenDiamondsTouchingAtOnePointPerPair() {
    checkS2PolygonUnchanged(
        "0:-70, -70:0, 0:70, 70:0; "
            + "0:-70, -60:0, 0:60, 60:0; "
            + "0:-50, -60:0, 0:50, 50:0; "
            + "0:-40, -40:0, 0:50, 40:0; "
            + "0:-30, -30:0, 0:30, 40:0; "
            + "0:-20, -20:0, 0:30, 20:0; "
            + "0:-10, -20:0, 0:10, 10:0");
  }

  /**
   * Checks that building the given S2TextFormat polygons in 'inputStrs' in an S2PolygonLayer with
   * the given EdgeType produces an S2Polygon with the same S2TextFormat as 'expectedStr'.
   */
  private void checkS2Polygon(List<String> inputStrs, String expectedStr, EdgeType edgeType) {
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    S2PolygonLayer layer = new S2PolygonLayer(new S2PolygonLayer.Options(edgeType));
    builder.startLayer(layer);
    boolean isFull = false;
    for (String inputStr : inputStrs) {
      if (inputStr.equals("full")) {
        isFull = true;
      }
      builder.addPolygon(makeVerbatimPolygonOrDie(inputStr));
    }
    builder.addIsFullPolygonPredicate(S2Builder.isFullPolygon(isFull));
    S2Error error = new S2Error();
    assertTrue(builder.build(error));
    // The input strings in tests may not be in normalized form, so we build an S2Polygon and
    // convert it back to a string.
    S2Polygon expected = S2TextFormat.makePolygonOrDie(expectedStr);
    assertEquals(S2TextFormat.toString(expected), S2TextFormat.toString(layer.getPolygon()));
  }

  /**
   * Checks that building the given S2TextFormat polygons in 'inputStrs' in an S2PolygonLayer with
   * either of DIRECTED or UNDIRECTED edges produces an S2Polygon with the same S2TextFormat as
   * 'expectedStr'.
   */
  private void checkS2Polygon(List<String> inputStrs, String expectedStr) {
    checkS2Polygon(inputStrs, expectedStr, EdgeType.DIRECTED);
    checkS2Polygon(inputStrs, expectedStr, EdgeType.UNDIRECTED);
  }

  /**
   * Checks that building the given S2TextFormat polygon in 'inputStr' in an S2PolygonLayer with
   * either of DIRECTED or UNDIRECTED edges produces the same polygon as output.
   */
  private void checkS2PolygonUnchanged(String inputStr) {
    checkS2Polygon(ImmutableList.of(inputStr), inputStr);
  }

  /**
   * Checks that building the given S2TextFormat *polylines* in 'inputStrs' into an S2PolygonLayer
   * with the given EdgeType results in an error with one of the 'expectedError' codes.
   */
  private void checkS2PolygonError(
      List<String> inputStrs, Set<S2Error.Code> expectedError, EdgeType edgeType) {
    S2Builder.Builder builderOptions = new S2Builder.Builder();
    S2Builder builder = builderOptions.build();

    S2PolygonLayer.Options layerOptions = new S2PolygonLayer.Options(edgeType);
    layerOptions.setValidate(true);

    builder.startLayer(new S2PolygonLayer(layerOptions));
    for (String inputStr : inputStrs) {
      builder.addPolyline(makePolylineOrDie(inputStr));
    }
    S2Error error = new S2Error();
    // Internal assertions must be disabled to allow invalid polygons to be built.
    boolean result = uncheckedCreate(() -> builder.build(error));
    assertFalse(result);
    assertTrue(expectedError.contains(error.code()));
  }

  private void checkS2PolygonError(List<String> inputStrs, S2Error.Code expectedError) {
    checkS2PolygonError(inputStrs, ImmutableSet.of(expectedError), EdgeType.DIRECTED);
    checkS2PolygonError(inputStrs, ImmutableSet.of(expectedError), EdgeType.UNDIRECTED);
  }

  private void checkS2PolygonError(List<String> inputStrs, Set<S2Error.Code> expectedError) {
    checkS2PolygonError(inputStrs, expectedError, EdgeType.DIRECTED);
    checkS2PolygonError(inputStrs, expectedError, EdgeType.UNDIRECTED);
  }

  /**
   * A map from a point representing an edge, to the set of labels on that edge.
   *
   * <p>Since we don't expect to have any crossing edges, the key for each edge is simply the
   * S2Point vector sum of its endpoints. This key has the advantage of being unchanged when the
   * endpoints of an edge are swapped.
   */
  private static class EdgeLabelMap extends HashMap<S2Point, HashSet<Integer>> {
    /** Gets the edge labels for the given key as an IdSet. */
    public IdSetLexicon.IdSet getIdSet(S2Point key) {
      return new IdSetLexicon.TestIdSet(get(key));
    }
  }

  /**
   * Adds the given polyline to the given builder, and also updates the provided EdgeLabelMap with
   * the labels for each added edge. Edges are identified by the combination of their (unordered)
   * endpoints.
   */
  private static void addPolylineWithLabels(
      S2Polyline polyline,
      EdgeType edgeType,
      int labelBegin,
      S2Builder builder,
      EdgeLabelMap edgeLabelMap) {
    for (int i = 0; i + 1 < polyline.numVertices(); ++i) {
      int label = labelBegin + i;
      builder.setLabel(label);
      // With undirected edges, reverse the direction of every other input edge.
      int dir = edgeType == EdgeType.DIRECTED ? 1 : (i & 1);
      S2Point a = polyline.vertex(i + (1 - dir));
      S2Point b = polyline.vertex(i + dir);
      builder.addEdge(a, b);
      System.err.println(
          "Edge "
              + a.toDegreesString()
              + " - "
              + b.toDegreesString()
              + " added with label "
              + label);
      S2Point key = S2Point.add(polyline.vertex(i), polyline.vertex(i + 1));
      edgeLabelMap.computeIfAbsent(key, (S2Point k) -> new HashSet<Integer>());
      edgeLabelMap.get(key).add(label);
    }
  }

  private static void checkEdgeLabels(EdgeType edgeType) {
    S2Builder.Builder builderOptions = new S2Builder.Builder();
    S2Builder builder = builderOptions.build();
    IntVector labelSetIds = new IntVector();
    IdSetLexicon labelSetLexicon = new IdSetLexicon();
    S2PolygonLayer layer =
        new S2PolygonLayer(new S2PolygonLayer.Options(edgeType), labelSetLexicon, labelSetIds);
    builder.startLayer(layer);

    // We use a polygon consisting of 3 loops, two of which are triangles and the other is a square.
    // The loops are reordered and some of the loops are inverted during S2Polygon construction.
    EdgeLabelMap edgeLabelMap = new EdgeLabelMap();
    addPolylineWithLabels(
        makePolylineOrDie("0:0, 9:1, 1:9, 0:0, 2:8, 8:2, 0:0, 0:10, 10:10, 10:0, 0:0"),
        edgeType,
        0,
        builder,
        edgeLabelMap);
    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));
    S2Polygon output = layer.getPolygon();
    S2Shape polygon = output.shape();

    S2Shape.MutableEdge edge = new S2Shape.MutableEdge();

    // The square is the first loop, as it encloses the other two loops.
    IntVector expectedLoopSizes = IntVector.of(4, 3, 3);

    // 10 edges in three loops.
    assertEquals(10, labelSetIds.size());
    assertEquals(expectedLoopSizes.size(), output.numLoops());
    for (int loopId = 0; loopId < output.numLoops(); ++loopId) {
      assertEquals(expectedLoopSizes.get(loopId), output.loop(loopId).numEdges());
    }

    // Check the label sets on each loop edge. This is visible only for testing.
    ArrayList<IntVector> labelSetIdsForLoops = layer.getLabelSetIdsForLoops();
    assertEquals(expectedLoopSizes.size(), labelSetIdsForLoops.size());

    for (int loopNum = 0; loopNum < expectedLoopSizes.size(); ++loopNum) {
      assertEquals(expectedLoopSizes.get(loopNum), labelSetIdsForLoops.get(loopNum).size());
      for (int loopOffset = 0; loopOffset < labelSetIdsForLoops.get(loopNum).size(); ++loopOffset) {
        // Note that even if the loop was reversed when constructing the polygon, we're still
        // iterating the edges in the forward order here, not using orientedVertex().
        S2Point a = output.loop(loopNum).vertex(loopOffset);
        S2Point b = output.loop(loopNum).vertex(loopOffset + 1);
        S2Point key = S2Point.add(a, b);
        IdSetLexicon.IdSet expectedLabelIds = edgeLabelMap.getIdSet(key);

        int idSetId = labelSetIdsForLoops.get(loopNum).get(loopOffset);
        IdSetLexicon.IdSet actualLabelIds = labelSetLexicon.idSet(idSetId);

        assertEquals(expectedLabelIds, actualLabelIds);
      }
    }

    // Check the label sets on each shape edge.
    for (int shapeEdgeId = 0; shapeEdgeId < polygon.numEdges(); ++shapeEdgeId) {
      int idSetId = labelSetIds.get(shapeEdgeId);
      IdSetLexicon.IdSet actualLabelIds = labelSetLexicon.idSet(idSetId);

      polygon.getEdge(shapeEdgeId, edge);
      S2Point key = S2Point.add(edge.getStart(), edge.getEnd());
      IdSetLexicon.IdSet expectedLabelIds = edgeLabelMap.getIdSet(key);

      assertEquals(expectedLabelIds, actualLabelIds);
    }
  }
}
