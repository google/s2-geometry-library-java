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

import static com.google.common.geometry.S2TextFormat.makeLaxPolygonOrDie;
import static com.google.common.geometry.S2TextFormat.makePointOrDie;
import static com.google.common.geometry.S2TextFormat.makePolylineOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2LaxPolygonLayer.DegenerateBoundaries;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import com.google.common.geometry.primitives.Ints.ImmutableIntSequence;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2LaxPolygonLayer}. */
@RunWith(JUnit4.class)
public final class S2LaxPolygonLayerTest extends GeometryTestCase {
  /** Returns a predicate that throws an exception if it is called. */
  private static S2Builder.IsFullPolygonPredicate throwsIfCalled() {
    return graph -> {
      throw new AssertionError("IsFullPolygonPredicate was called, should not be.");
    };
  }

  private void checkLaxPolygon(
      String inputStr,
      String expectedStr,
      EdgeType edgeType,
      DegenerateBoundaries degenerateBoundaries) {
    S2Builder builder = new S2Builder.Builder().build();
    S2LaxPolygonLayer.Options options = new S2LaxPolygonLayer.Options();
    options.setEdgeType(edgeType);
    options.setDegenerateBoundaries(degenerateBoundaries);
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer(options);
    builder.startLayer(layer);

    S2LaxPolygonShape polygon = makeLaxPolygonOrDie(inputStr);
    builder.addShape(polygon);

    // In order to construct polygons that are full except possibly for a collection of degenerate
    // holes, we must supply S2Builder with a predicate that distinguishes empty polygons from full
    // ones (modulo degeneracies).
    boolean hasFullLoop = false;
    for (int i = 0; i < polygon.numChains(); ++i) {
      if (polygon.getChainLength(i) == 0) {
        hasFullLoop = true;
      }
    }
    builder.addIsFullPolygonPredicate(S2Builder.isFullPolygon(hasFullLoop));
    S2Error error = new S2Error();
    assertTrue(builder.build(error));
    S2LaxPolygonShape output = layer.getPolygon();
    String actualStr = S2TextFormat.toString(output, "; ");
    assertEquals(expectedStr, actualStr);
  }

  private void checkLaxPolygon(
      String inputStr, String expectedStr, DegenerateBoundaries degenerateBoundaries) {
    checkLaxPolygon(inputStr, expectedStr, EdgeType.DIRECTED, degenerateBoundaries);

    // TODO(user): Implement.
    // checkLaxPolygon(inputStr, expectedStr, EdgeType.UNDIRECTED, degenerateBoundaries);
  }

  void checkLaxPolygonUnchanged(String inputStr, DegenerateBoundaries degenerateBoundaries) {
    checkLaxPolygon(inputStr, inputStr, degenerateBoundaries);
  }

  private static final ImmutableList<DegenerateBoundaries> ALL_DEGENERATE_BOUNDARIES =
      ImmutableList.of(
          DegenerateBoundaries.DISCARD,
          DegenerateBoundaries.DISCARD_HOLES,
          DegenerateBoundaries.DISCARD_SHELLS,
          DegenerateBoundaries.KEEP);

  @Test
  public void testEmpty() {
    for (DegenerateBoundaries db : ALL_DEGENERATE_BOUNDARIES) {
      checkLaxPolygonUnchanged("", db);
    }
  }

  @Test
  public void testFull() {
    for (DegenerateBoundaries db : ALL_DEGENERATE_BOUNDARIES) {
      checkLaxPolygonUnchanged("full", db);
    }
  }

  @Test
  public void testOneNormalShell() {
    for (DegenerateBoundaries db : ALL_DEGENERATE_BOUNDARIES) {
      checkLaxPolygonUnchanged("0:0, 0:1, 1:1", db);
    }
  }

  @Test
  public void testIsFullPolygonPredicateNotCalled() {
    // Test that the IsFullPolygonPredicate is not called when at least one non-degenerate loop is
    // present.
    for (DegenerateBoundaries degenerateBoundaries : ALL_DEGENERATE_BOUNDARIES) {
      S2Builder builder = new S2Builder.Builder().build();
      S2LaxPolygonLayer.Options options = new S2LaxPolygonLayer.Options();
      options.setEdgeType(EdgeType.DIRECTED);
      options.setDegenerateBoundaries(degenerateBoundaries);
      S2LaxPolygonLayer layer = new S2LaxPolygonLayer(options);
      builder.startLayer(layer);
      S2LaxPolygonShape polygon = makeLaxPolygonOrDie("0:0, 0:1, 1:1");
      builder.addShape(polygon);
      // If the predicate is called, it will throw an exception.
      builder.addIsFullPolygonPredicate(throwsIfCalled());
      S2Error error = new S2Error();
      assertTrue(builder.build(error));
    }
  }

  @Test
  public void testSingleDegenerateEdgeBecomesShell() {
    S2Builder builder = new S2Builder.Builder().build();
    S2LaxPolygonLayer.Options options = new S2LaxPolygonLayer.Options();
    options.setEdgeType(EdgeType.DIRECTED);
    options.setDegenerateBoundaries(DegenerateBoundaries.KEEP);
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer(options);
    builder.startLayer(layer);
    // A loop with one vertex defines a single degenerate edge.
    S2Point p = makePointOrDie("1:1");
    builder.addPoint(p);

    // The S2LaxPolygonLayer needs to determine if the degenerate edge is intended to represent a
    // degenerate hole in the full polygon, or simply a degenerate shell. It calls the predicate to
    // decide that, which we set to be "false". Without a full loop to contain it, the edge is a
    // shell.
    builder.addIsFullPolygonPredicate(g -> false);

    S2Error error = new S2Error();
    boolean result = builder.build(error);
    assertTrue(error.text(), result);

    // The point should have become a degenerate shell.
    S2LaxPolygonShape polygon = layer.getPolygon();
    assertEquals(1, polygon.numChains());
    assertEquals(1, polygon.chain(0).size());
    assertTrue(polygon.chain(0).get(0).equalsPoint(p));
  }

  @Test
  public void testDegenerateEdgeBecomesHoleInSphere() {
    S2Builder builder = new S2Builder.Builder().build();
    S2LaxPolygonLayer.Options options = new S2LaxPolygonLayer.Options();
    options.setEdgeType(EdgeType.DIRECTED);
    options.setDegenerateBoundaries(DegenerateBoundaries.KEEP);
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer(options);
    builder.startLayer(layer);
    S2Point p = makePointOrDie("1:1");
    builder.addPoint(p);

    // Just like the test case above, but now the degenerate edge becomes a hole because it is
    // contained by the full polygon.
    builder.addIsFullPolygonPredicate(g -> true);

    S2Error error = new S2Error();
    boolean result = builder.build(error);
    assertTrue(error.text(), result);

    // The result has two loops: an empty one to define the full sphere, and the degenerate hole.
    S2LaxPolygonShape polygon = layer.getPolygon();
    assertEquals(2, polygon.numChains());
    assertEquals(0, polygon.getChainLength(0)); // empty chain for full loop
    assertEquals(1, polygon.getChainLength(1)); // degenerate loop
    String s = S2TextFormat.toString(polygon, "; ");
    assertTrue(s.equals("full; 1:1"));
  }

  @Test
  public void testDegenerateEdgeBecomesHoleInLoop() {
    S2Builder builder = new S2Builder.Builder().build();
    S2LaxPolygonLayer.Options options = new S2LaxPolygonLayer.Options();
    options.setEdgeType(EdgeType.DIRECTED);
    options.setDegenerateBoundaries(DegenerateBoundaries.KEEP);
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer(options);
    builder.startLayer(layer);
    S2Loop loop = S2TextFormat.makeLoopOrDie("0:0, 0:2, 2:2, 2:0");
    builder.addShape(loop);
    S2Point p = makePointOrDie("1:1");
    builder.addPoint(p);

    // In this case the layer doesn't need to check the full polygon predicate.
    builder.addIsFullPolygonPredicate(throwsIfCalled());

    S2Error error = new S2Error();
    boolean result = builder.build(error);
    assertTrue(error.text(), result);

    S2LaxPolygonShape polygon = layer.getPolygon();
    assertEquals(2, polygon.numChains());
    assertEquals(4, polygon.getChainLength(0));
    assertEquals(1, polygon.getChainLength(1));
    String s = S2TextFormat.toString(polygon, "; ");
    assertTrue(s.equals("0:0, 0:2, 2:2, 2:0; 1:1"));
  }

  @Test
  public void testTwoNormalShellsOneNormalHole() {
    // The second two loops are nested. Note that S2LaxPolygon and S2Polygon require opposite vertex
    // orderings for holes.
    for (DegenerateBoundaries db : ALL_DEGENERATE_BOUNDARIES) {
      checkLaxPolygonUnchanged(
          "0:1, 1:1, 0:0; " + "3:3, 3:6, 6:6, 6:3; " + "4:4, 5:4, 5:5, 4:5", db);
    }
  }

  @Test
  public void testAllDegenerateShells() {
    for (DegenerateBoundaries db :
        ImmutableList.of(DegenerateBoundaries.KEEP, DegenerateBoundaries.DISCARD_HOLES)) {
      checkLaxPolygonUnchanged("1:1; 2:2, 3:3", db);
    }
    for (DegenerateBoundaries db :
        ImmutableList.of(DegenerateBoundaries.DISCARD, DegenerateBoundaries.DISCARD_SHELLS)) {
      checkLaxPolygon("1:1; 2:2, 3:3", "", db);
    }
  }

  @Test
  public void testAllDegenerateHoles() {
    for (DegenerateBoundaries db :
        ImmutableList.of(DegenerateBoundaries.KEEP, DegenerateBoundaries.DISCARD_SHELLS)) {
      checkLaxPolygonUnchanged("full; 1:1; 2:2, 3:3", db);
    }
    for (DegenerateBoundaries db :
        ImmutableList.of(DegenerateBoundaries.DISCARD, DegenerateBoundaries.DISCARD_HOLES)) {
      checkLaxPolygon("full; 1:1; 2:2, 3:3", "full", db);
    }
  }

  @Test
  public void testSomeDegenerateShells() {
    String kNormal = "0:0, 0:9, 9:0; 1:1, 7:1, 1:7";
    String kInput = kNormal + "; 3:2; 2:2, 2:3";
    checkLaxPolygonUnchanged(kInput, DegenerateBoundaries.KEEP);
    checkLaxPolygonUnchanged(kInput, DegenerateBoundaries.DISCARD_HOLES);
    checkLaxPolygon(kInput, kNormal, DegenerateBoundaries.DISCARD);
    checkLaxPolygon(kInput, kNormal, DegenerateBoundaries.DISCARD_SHELLS);
  }

  @Test
  public void testSomeDegenerateHoles() {
    for (DegenerateBoundaries db :
        ImmutableList.of(DegenerateBoundaries.KEEP, DegenerateBoundaries.DISCARD_SHELLS)) {
      checkLaxPolygonUnchanged("0:0, 0:9, 9:0; 1:1; 2:2, 3:3", db);
    }
    for (DegenerateBoundaries db :
        ImmutableList.of(DegenerateBoundaries.DISCARD, DegenerateBoundaries.DISCARD_HOLES)) {
      checkLaxPolygon("0:0, 0:9, 9:0; 1:1; 2:2, 3:3", "0:0, 0:9, 9:0", db);
    }
  }

  @Test
  public void testNormalAndDegenerateShellsAndHoles() {
    // We start with two normal shells and one normal hole.
    String kNormal = "0:0, 0:9, 9:9, 9:0; " + "0:10, 0:19, 9:19, 9:10; 1:11, 8:11, 8:18, 1:18";
    // These are the same loops augmented with degenerate interior filaments (holes). Note that one
    // filament connects the second shell and hole above, transforming them into a single loop.
    String kNormalWithDegenHoles =
        "0:0, 0:9, 1:8, 1:7, 1:8, 0:9, 9:9, 9:0; "
            + "0:10, 0:19, 9:19, 9:10, 0:10, 1:11, 8:11, 8:18, 1:18, 1:11";
    // Then we add other degenerate shells and holes, including a sibling pair that connects the two
    // shells above.
    String kDegenShells = "0:9, 0:10; 2:12; 3:13, 3:14; 20:20; 10:0, 10:1";
    String kDegenHoles = "2:5; 3:6, 3:7; 8:8";
    String kInput = kNormalWithDegenHoles + "; " + kDegenShells + "; " + kDegenHoles;
    checkLaxPolygon(kInput, kNormal, DegenerateBoundaries.DISCARD);
    checkLaxPolygon(kInput, kNormal + "; " + kDegenShells, DegenerateBoundaries.DISCARD_HOLES);
    checkLaxPolygon(
        kInput, kNormalWithDegenHoles + "; " + kDegenHoles, DegenerateBoundaries.DISCARD_SHELLS);
    checkLaxPolygon(kInput, kInput, DegenerateBoundaries.KEEP);
  }

  @Test
  public void testPartialLoop() {
    S2Builder builder = new S2Builder.Builder().build();
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer();
    builder.startLayer(layer);
    builder.addPolyline(makePolylineOrDie("0:1, 2:3, 4:5"));
    S2Error error = new S2Error();
    assertFalse(builder.build(error));
    assertEquals(S2Error.Code.BUILDER_EDGES_DO_NOT_FORM_LOOPS, error.code());
    try {
      S2LaxPolygonShape unused = layer.getPolygon();
      fail("Expected exception was not thrown.");
    } catch (NullPointerException expected) {
      // test passes
    }
  }

  // TODO(user): Implement validation of S2LaxPolygonShape.
  /*
   @Test
  public void testInvalidPolygon() {
    S2Builder builder = new S2Builder.Builder().build();
    S2LaxPolygonLayer.Options options = new S2LaxPolygonLayer.Options();
    // options.setValidate(true);
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer(options);
    builder.startLayer(layer);
    builder.addPolyline(makePolylineOrDie("0:0, 0:10, 10:0, 10:10, 0:0"));
    S2Error error = new S2Error();
    assertFalse(builder.build(error));
    assertEquals(S2Error.Code.LOOP_SELF_INTERSECTION, error.code());
  }
  */

  /**
   * Check that S2LaxPolygonLayer removes duplicate edges in such a way that degeneracies are not
   * lost.
   */
  @Test
  public void testDuplicateInputEdges() {
    S2Builder builder = new S2Builder.Builder().build();
    S2LaxPolygonLayer.Options options = new S2LaxPolygonLayer.Options();
    options.setDegenerateBoundaries(DegenerateBoundaries.KEEP);
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer(options);
    builder.startLayer(layer);
    builder.addShape(makeLaxPolygonOrDie("0:0, 0:5, 5:5, 5:0"));
    builder.addPoint(makePointOrDie("0:0"));
    builder.addPoint(makePointOrDie("1:1"));
    builder.addPoint(makePointOrDie("1:1"));
    builder.addShape(makeLaxPolygonOrDie("2:2, 2:3"));
    builder.addShape(makeLaxPolygonOrDie("2:2, 2:3"));
    S2Error error = new S2Error();
    assertTrue(builder.build(error));
    S2LaxPolygonShape output = layer.getPolygon();
    assertEquals("0:0, 0:5, 5:5, 5:0; 1:1; 2:2, 2:3", S2TextFormat.toString(output, "; "));
  }

  /** For undirected edges, sort the vertices in lexicographic order. */
  private static ImmutableEdge getKey(S2Shape.MutableEdge edge, EdgeType edgeType) {
    if (edgeType == EdgeType.UNDIRECTED && edge.b.lessThan(edge.a)) {
      edge.reverse();
    }
    return new ImmutableEdge(edge);
  }

  private static void addShapeWithLabels(
      S2Shape shape,
      EdgeType edgeType,
      S2Builder builder,
      Map<ImmutableEdge, IntVector> edgeLabelMap) {
    int kLabelBegin = 1234; // Arbitrary.

    for (int e = 0; e < shape.numEdges(); ++e) {
      int label = kLabelBegin + e;
      builder.setLabel(label);
      // For undirected edges, reverse the direction of every other input edge.
      S2Shape.MutableEdge edge = new S2Shape.MutableEdge();
      shape.getEdge(e, edge);
      if (edgeType == EdgeType.UNDIRECTED && ((e & 1) != 0)) {
        edge.reverse();
      }
      builder.addEdge(edge.a, edge.b);
      ImmutableEdge key = getKey(edge, edgeType);
      edgeLabelMap.putIfAbsent(key, new IntVector());
      edgeLabelMap.get(key).add(label);
    }
  }

  /**
   * Converts "inputStr" to an S2LaxPolygonShape, assigns labels to its edges, then uses
   * S2LaxPolygonLayer with the given arguments to build a new S2LaxPolygonShape and verifies that
   * all edges have the expected labels. (This function does not test whether the output edges are
   * correct.)
   */
  private static void checkEdgeLabels(
      String inputStr, EdgeType edgeType, DegenerateBoundaries degenerateBoundaries) {
    S2Builder builder = new S2Builder.Builder().build();
    ArrayList<IntVector> labelSetIds = new ArrayList<>();
    IdSetLexicon labelSetLexicon = new IdSetLexicon();
    S2LaxPolygonLayer.Options options = new S2LaxPolygonLayer.Options();
    options.setEdgeType(edgeType);
    options.setDegenerateBoundaries(degenerateBoundaries);
    S2LaxPolygonLayer layer = new S2LaxPolygonLayer(options, labelSetLexicon, labelSetIds);
    builder.startLayer(layer);

    HashMap<ImmutableEdge, IntVector> edgeLabelMap = new HashMap<>();
    addShapeWithLabels(makeLaxPolygonOrDie(inputStr), edgeType, builder, edgeLabelMap);
    // Sort and deduplicate all the IntVectors in the edgeLabelMap
    for (IntVector v : edgeLabelMap.values()) {
      v.sort();
      v.unique();
    }

    S2Error error = new S2Error();
    assertTrue(builder.build(error));
    S2LaxPolygonShape output = layer.getPolygon();
    S2Shape.MutableEdge edge = new S2Shape.MutableEdge();

    for (int i = 0; i < output.numChains(); ++i) {
      for (int j = 0; j < output.chain(i).size(); ++j) {
        output.getChainEdge(i, j, edge);
        ImmutableEdge key = getKey(edge, edgeType);
        ImmutableIntSequence expectedLabels = ImmutableIntSequence.viewOf(edgeLabelMap.get(key));
        ImmutableIntSequence actualLabels =
            ImmutableIntSequence.viewOf(labelSetLexicon.idSet(labelSetIds.get(i).get(j)));
        assertEquals(
            Platform.formatString(
                "Chain %d, chainEdge %d is edge %s. EdgeType %s for key %s, expected labels %s.\n",
                i, j, edge, edgeType, key, expectedLabels),
            expectedLabels,
            actualLabels);
      }
    }
  }

  @Test
  public void testEdgeLabels() {
    // TODO(user): Implement EdgeType.UNDIRECTED.
    for (EdgeType edgeType : ImmutableList.of(EdgeType.DIRECTED)) {
      for (DegenerateBoundaries db : ALL_DEGENERATE_BOUNDARIES) {
        // Test a polygon with normal and degenerate shells and holes. Note that this
        // S2LaxPolygonShape has duplicate edges and is therefore not valid in most contexts.
        checkEdgeLabels(
            "1:1, 1:2; 0:0, 0:9, 9:9, 9:0; 1:2, 1:1; 3:3, 8:3, 8:8, 3:8; 4:4; 4:5, 5:5; 4:4",
            edgeType,
            db);
      }
    }
  }
}
