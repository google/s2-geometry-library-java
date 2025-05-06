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
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2BuilderGraph.PolylineType;
import com.google.common.geometry.S2BuilderSnapFunctions.IntLatLngSnapFunction;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2PolylineVectorLayer}. */
@RunWith(JUnit4.class)
public final class S2PolylineVectorLayerTest extends GeometryTestCase {

  @Test
  public void testNoEdges() {
    checkS2PolylineVectorUnchanged(ImmutableList.of());
  }

  @Test
  public void testSingleVertexPolylines() {
    // One polyline with a single vertex.
    checkS2PolylineVectorUnchanged(ImmutableList.of("2:2"));

    // One polyline with two vertices, one with a single vertex, not touching.
    checkS2PolylineVectorUnchanged(ImmutableList.of("2:2", "1:1, 1:0"));
  }

  /**
   * Check that polylines are joined together when possible, even if they were not adjacent in the
   * input. For undirected edges, the polyline direction should be chosen such that the first edge
   * of the polyline was added to S2Builder before the last edge of the polyline.
   */
  @Test
  public void testJoiningPolylines() {
    checkS2PolylineVector(
        ImmutableList.of("1:1, 2:2", "3:3, 2:2", "0:0, 1:1"),
        ImmutableList.of("3:3, 2:2", "0:0, 1:1, 2:2"),
        EdgeType.DIRECTED,
        null,
        null);
    checkS2PolylineVector(
        ImmutableList.of("1:1, 2:2", "3:3, 2:2", "0:0, 1:1"),
        ImmutableList.of("3:3, 2:2, 1:1, 0:0"),
        EdgeType.UNDIRECTED,
        null,
        null);
    // A case where one of the input polylines is degenerate
    checkS2PolylineVector(
        ImmutableList.of("1:1, 2:2", "3:3, 2:2", "1:1"),
        ImmutableList.of("1:1, 2:2, 3:3"),
        EdgeType.UNDIRECTED,
        null,
        null);
  }

  /** Test a complex network of polylines that meet at shared vertices. */
  @Test
  public void testSegmentNetwork() {
    checkS2PolylineVectorUnchanged(
        ImmutableList.of(
            "0:0, 1:1, 2:2",
            "2:2, 2:3, 2:4",
            "2:4, 3:4, 4:4",
            "2:2, 3:2, 4:2",
            "4:2, 4:3, 4:4",
            "1:0, 2:2",
            "0:1, 2:2",
            "5:4, 4:4",
            "4:5, 4:4",
            "2:4, 2:5, 1:5, 1:4, 2:4",
            "4:2, 6:1, 5:0", // Two nested loops
            "4:2, 7:0, 6:-1",
            "11:1, 11:0, 10:0, 10:1, 11:1" // Isolated loop
            ));
  }

  /**
   * This checks idempotency for directed edges in the case of several polylines that share edges
   * (and that even share loops). The test happens to pass for undirected edges as well.
   */
  @Test
  public void testMultipleIntersectingWalks() {
    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setPolylineType(PolylineType.WALK);
    ImmutableList<String> input =
        ImmutableList.of(
            "5:5, 5:6, 6:5, 5:5, 5:4, 5:3",
            "4:4, 5:5, 6:5, 5:6, 5:5, 5:6, 6:5, 5:5, 4:5",
            "3:5, 5:5, 5:6, 6:5, 5:5, 5:6, 6:6, 7:7");
    checkS2PolylineVector(input, input, layerOptions, null);
  }

  /**
   * This checks idempotency for cases where earlier polylines in the input happen to terminate in
   * the middle of later polylines. This requires building non-maximal polylines.
   */
  @Test
  public void testEarlyWalkTermination() {
    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setPolylineType(PolylineType.WALK);
    ImmutableList<String> input =
        ImmutableList.of("0:1, 1:1", "1:0, 1:1, 1:2", "0:2, 1:2, 2:2", "2:1, 2:2, 2:3");
    checkS2PolylineVector(input, input, layerOptions, null);
  }

  /**
   * A single input edge is split into several segments by removing portions of it, and then each of
   * those segments becomes one edge of a loop.
   */
  @Test
  public void testInputEdgeStartsMultipleLoops() {
    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setPolylineType(PolylineType.WALK).setSiblingPairs(SiblingPairs.DISCARD);

    S2Builder.Builder builderOptions = new S2Builder.Builder();
    builderOptions.setSnapFunction(new IntLatLngSnapFunction(7));
    ImmutableList<String> input =
        ImmutableList.of(
            "0:10, 0:0",
            "0:6, 1:6, 1:7, 0:7, 0:8",
            "0:8, 1:8, 1:9, 0:9, 0:10",
            "0:2, 1:2, 1:3, 0:3, 0:4",
            "0:0, 1:0, 1:1, 0:1, 0:2",
            "0:4, 1:4, 1:5, 0:5, 0:6");
    ImmutableList<String> expected =
        ImmutableList.of(
            "0:1, 0:0, 1:0, 1:1, 0:1",
            "0:3, 0:2, 1:2, 1:3, 0:3",
            "0:5, 0:4, 1:4, 1:5, 0:5",
            "0:7, 0:6, 1:6, 1:7, 0:7",
            "0:9, 0:8, 1:8, 1:9, 0:9");
    checkS2PolylineVector(input, expected, layerOptions, builderOptions);
  }

  /** Verifies that the setValidate() option works. */
  @Test
  public void testValidateTrue() {
    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setValidate(true);

    S2Builder.Builder builderOptions = new S2Builder.Builder();
    S2Builder builder = builderOptions.build();

    S2PolylineVectorLayer layer = new S2PolylineVectorLayer(layerOptions);
    builder.startLayer(layer);
    builder.addEdge(new S2Point(1, 0, 0), new S2Point(-1, 0, 0));
    S2Error error = new S2Error();
    assertFalse(builder.build(error));
    assertEquals(S2Error.Code.ANTIPODAL_VERTICES, error.code());
  }

  /**
   * Tests setting labels on edges, building an S2PolylineVectorLayer with setWithEdgeLabels(true),
   * and retrieving the labels
   */
  @Test
  public void testSimpleEdgeLabels() {
    S2Builder.Builder builderOptions = new S2Builder.Builder();
    S2Builder builder = builderOptions.build();

    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setEdgeType(EdgeType.UNDIRECTED).setDuplicateEdges(DuplicateEdges.MERGE);
    List<IntVector> labelSetIds = new ArrayList<>();
    IdSetLexicon labelSetLexicon = new IdSetLexicon();

    S2PolylineVectorLayer layer =
        new S2PolylineVectorLayer(layerOptions, labelSetLexicon, labelSetIds);

    builder.startLayer(layer);
    builder.setLabel(1);
    builder.addPolyline(makePolylineOrDie("0:0, 0:1, 0:2"));
    builder.setLabel(2);
    builder.addPolyline(makePolylineOrDie("0:3, 0:2, 0:1"));
    builder.clearLabels();
    builder.addPolyline(makePolylineOrDie("0:4, 0:5"));
    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));

    ImmutableList<List<List<Integer>>> expected =
        ImmutableList.of(
            ImmutableList.of(ImmutableList.of(1), ImmutableList.of(1, 2), ImmutableList.of(2)),
            ImmutableList.of(ImmutableList.of()));
    assertEquals(expected.size(), labelSetIds.size());
    for (int i = 0; i < expected.size(); ++i) {
      assertEquals(expected.get(i).size(), labelSetIds.get(i).size());
      for (int j = 0; j < expected.get(i).size(); ++j) {
        assertEquals(
            expected.get(i).get(j).size(), labelSetLexicon.idSet(labelSetIds.get(i).get(j)).size());
        int k = 0;
        for (OfInt iter = labelSetLexicon.idSet(labelSetIds.get(i).get(j)).intIterator();
            iter.hasNext(); ) {
          int label = iter.nextInt();
          assertEquals((int) expected.get(i).get(j).get(k++), label);
        }
      }
    }
  }

  private void checkS2PolylineVector(
      List<String> inputStrs,
      List<String> expectedStrs,
      EdgeType edgeType,
      S2PolylineVectorLayer.Options layerOptions,
      S2Builder.Builder builderOptions) {
    if (layerOptions == null) {
      layerOptions = new S2PolylineVectorLayer.Options();
    }
    if (builderOptions == null) {
      builderOptions = new S2Builder.Builder();
    }

    layerOptions.setEdgeType(edgeType);
    S2Builder builder = builderOptions.build();
    S2PolylineVectorLayer layer = new S2PolylineVectorLayer(layerOptions);
    builder.startLayer(layer);
    for (String inputStr : inputStrs) {
      builder.addPolyline(makePolylineOrDie(inputStr));
    }
    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));
    List<S2Polyline> output = layer.getPolylines();
    assertEquals(expectedStrs.size(), output.size());

    for (int i = 0; i < output.size(); i++) {
      // The input strings in tests may not be in normalized form, so we build an S2Polyline and
      // convert it back to a string.
      S2Polyline expected = S2TextFormat.makePolylineOrDie(expectedStrs.get(i));
      S2Polyline actual = output.get(i);
      assertEquals(S2TextFormat.toString(expected), S2TextFormat.toString(actual));
    }
  }

  // Convenience function that tests both directed and undirected edges.
  private void checkS2PolylineVector(
      List<String> inputStrs,
      List<String> expectedStrs,
      S2PolylineVectorLayer.Options layerOptions,
      S2Builder.Builder builderOptions) {
    checkS2PolylineVector(inputStrs, expectedStrs, EdgeType.DIRECTED, layerOptions, builderOptions);
    checkS2PolylineVector(
        inputStrs, expectedStrs, EdgeType.UNDIRECTED, layerOptions, builderOptions);
  }

  private void checkS2PolylineVectorUnchanged(List<String> inputStrs) {
    checkS2PolylineVector(inputStrs, inputStrs, null, null);
  }
}
