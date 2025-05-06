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

import static com.google.common.geometry.S2TextFormat.makeLaxPolylineOrDie;
import static com.google.common.geometry.S2TextFormat.makePolylineOrDie;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2BuilderGraph.DegenerateBoundaries;
import com.google.common.geometry.S2BuilderGraph.DirectedComponent;
import com.google.common.geometry.S2BuilderGraph.Edge;
import com.google.common.geometry.S2BuilderGraph.EdgeList;
import com.google.common.geometry.S2BuilderGraph.LoopType;
import com.google.common.geometry.S2BuilderGraph.PolylineBuilder;
import com.google.common.geometry.S2BuilderGraph.UndirectedComponent;
import com.google.common.geometry.S2BuilderGraph.VertexOutEdgeIds;
import com.google.common.geometry.S2BuilderGraph.VertexOutEdges;
import com.google.common.geometry.S2BuilderGraph.VertexOutMap;
import com.google.common.geometry.S2BuilderUtil.GraphClone;
import com.google.common.geometry.S2BuilderUtil.GraphCloningLayer;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import com.google.common.geometry.primitives.Ints.ImmutableIntSequence;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/**
 * Unit tests for S2BuilderGraph. More extensive testing of S2BuilderGraph is via the Layer
 * implementation unit tests.
 */
@RunWith(JUnit4.class)
public final class S2BuilderGraphTest extends GeometryTestCase {

  /** Test default values and accessors of GraphOptions. */
  @Test
  public void testGraphOptions() {
    GraphOptions options = new GraphOptions();
    // Default values
    assertEquals(EdgeType.DIRECTED, options.edgeType());
    assertEquals(DegenerateEdges.KEEP, options.degenerateEdges());
    assertEquals(DuplicateEdges.KEEP, options.duplicateEdges());
    assertEquals(SiblingPairs.KEEP, options.siblingPairs());
    assertTrue(options.allowVertexFiltering());

    // Not the defaults.
    options =
        new GraphOptions(
            EdgeType.UNDIRECTED,
            DegenerateEdges.DISCARD,
            DuplicateEdges.MERGE,
            SiblingPairs.DISCARD_EXCESS);
    options.setAllowVertexFiltering(false);
    assertEquals(EdgeType.UNDIRECTED, options.edgeType());
    assertEquals(DegenerateEdges.DISCARD, options.degenerateEdges());
    assertEquals(DuplicateEdges.MERGE, options.duplicateEdges());
    assertEquals(SiblingPairs.DISCARD_EXCESS, options.siblingPairs());
    assertFalse(options.allowVertexFiltering());
  }

  /** Construct a graph and check basic accessors. */
  @Test
  public void testConstructAndAccessGraph() {
    GraphOptions options = new GraphOptions();
    ArrayList<S2Point> vertices = new ArrayList<>();
    vertices.add(new S2Point(1, 0, 0));
    EdgeList edges = new EdgeList();
    edges.add(37, 42);
    IntVector inputEdgeIdSetIds = IntVector.empty();
    inputEdgeIdSetIds.add(0);
    IdSetLexicon inputEdgeIdSetLexicon = new IdSetLexicon();
    IdSetLexicon labelSetLexicon = new IdSetLexicon();
    IntVector labelSetIds = IntVector.empty(); // Empty means no labels are present.
    S2BuilderGraph g =
        new S2BuilderGraph(
            options,
            vertices,
            edges,
            inputEdgeIdSetIds,
            inputEdgeIdSetLexicon,
            labelSetIds,
            labelSetLexicon,
            null);

    assertEquals(1, g.numVertices());
    assertEquals(1, g.numEdges());
    S2Point p = g.vertex(0);
    assertEquals(p, new S2Point(1, 0, 0));

    assertEquals(37, g.edges().getSrcId(0));
    assertEquals(42, g.edges().getDstId(0));

    // TODO(torrey): Add more accessor checks when possible.
  }

  /** Tests the situation where labels are requested but none were provided. */
  @Test
  public void testLabelsRequestedButNotProvided() {
    GraphOptions options =
        new GraphOptions(
            EdgeType.DIRECTED,
            DegenerateEdges.KEEP,
            DuplicateEdges.KEEP,
            SiblingPairs.KEEP);

    // A single vertex and a degenerate edge which is a loop at that vertex.
    ImmutableList<S2Point> vertices = ImmutableList.of(new S2Point(1, 0, 0));
    EdgeList edges = new EdgeList();
    edges.add(0, 0);

    // inputEdgeIdSetIds maps an edge id to an InputEdgeIdSetId (an int) which is in turn a key in
    // inputEdgeIdSetLexicon.
    IntVector inputEdgeIdSetIds = IntVector.of(0);
    // The inputEdgeIdSetLexicon is empty, so no input edges were mapped to Edge 0.
    IdSetLexicon inputEdgeIdSetLexicon = new IdSetLexicon();
    // An empty labelSetIds means no labels are present.
    IntVector labelSetIds = new IntVector();
    // The labelSetLexicon maps a labelSetId (an int) to a set of S2Builder labels (also ints).
    IdSetLexicon labelSetLexicon = new IdSetLexicon();
    // Construct the graph.
    S2BuilderGraph g =
        new S2BuilderGraph(
            options,
            vertices,
            edges,
            inputEdgeIdSetIds,
            inputEdgeIdSetLexicon,
            labelSetIds,
            labelSetLexicon,
            null);

    assertTrue(g.labelSetIds().isEmpty());
    assertEquals(IdSetLexicon.EMPTY_SET_ID, g.labelSetId(0));
    assertEquals(0, g.labels(0).size()); // Labels for input edge 0.
    S2BuilderGraph.LabelFetcher fetcher = new S2BuilderGraph.LabelFetcher(g, EdgeType.DIRECTED);
    IntVector labels = new IntVector();
    fetcher.fetch(0, labels); // Labels for graph edge 0.
    assertTrue(labels.isEmpty());
  }

  /** Tests S2BuilderGraph.getDirectedLoops() with a simple input graph. */
  @Test
  public void testGetDirectedLoopsSimple() {
    GraphClone graphClone = new GraphClone();
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    GraphOptions graphOptions =
        new GraphOptions(
            EdgeType.DIRECTED,
            DegenerateEdges.DISCARD_EXCESS,
            DuplicateEdges.KEEP,
            SiblingPairs.KEEP);

    builder.startLayer(new GraphCloningLayer(graphOptions, graphClone));
    builder.addShape(makeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0")); // A-B-C-D-A, square

    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));

    S2BuilderGraph graph = graphClone.graph();
    ArrayList<int[]> loops = new ArrayList<>();
    assertTrue(graph.getDirectedLoops(S2BuilderGraph.LoopType.SIMPLE, loops, error));
    assertEquals(1, loops.size());
    assertEquals(4, loops.get(0).length);
  }

  /** Tests S2BuilderGraph.getDirectedLoops() with an input graph that contains degenerate edges. */
  @Test
  public void testGetDirectedLoopsDegenerateEdges() {
    GraphClone graphClone = new GraphClone();
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    GraphOptions graphOptions =
        new GraphOptions(
            EdgeType.DIRECTED,
            DegenerateEdges.DISCARD_EXCESS,
            DuplicateEdges.KEEP,
            SiblingPairs.KEEP);

    builder.startLayer(new GraphCloningLayer(graphOptions, graphClone));
    builder.addShape(makeLaxPolylineOrDie("1:1, 1:1")); // A-A, single point
    builder.addShape(makeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0")); // A-B-C-D-A, square
    builder.addShape(makeLaxPolylineOrDie("0:3, 3:3, 0:3")); // A-B-A, two points

    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));

    S2BuilderGraph graph = graphClone.graph();
    ArrayList<int[]> loops = new ArrayList<>();
    assertTrue(
        error.toString(), graph.getDirectedLoops(S2BuilderGraph.LoopType.SIMPLE, loops, error));
    assertEquals(3, loops.size());
    assertEquals(1, loops.get(0).length);
    assertEquals(4, loops.get(1).length);
    assertEquals(2, loops.get(2).length);
  }

  /** Tests S2BuilderGraph.getDirectedComponents() with a simple input graph. */
  @Test
  public void testGetDirectedComponentsSimple() {
    GraphClone graphClone = new GraphClone();
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    GraphOptions graphOptions =
        new GraphOptions(
            EdgeType.DIRECTED,
            DegenerateEdges.DISCARD_EXCESS,
            DuplicateEdges.KEEP,
            SiblingPairs.CREATE);

    builder.startLayer(new GraphCloningLayer(graphOptions, graphClone));
    builder.addShape(makeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0")); // A-B-C-D-A square

    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));

    S2BuilderGraph graph = graphClone.graph();
    ArrayList<DirectedComponent> components = new ArrayList<>();
    assertTrue(
        error.toString(),
        graph.getDirectedComponents(DegenerateBoundaries.KEEP, components, error));

    // Since sibling edges are created, the result should be a single directed component, with two
    // edge loops: one going forward around the provided ABCDA square, and the other edge loop going
    // backward around the same square.
    assertEquals(1, components.size());
    assertEquals(2, components.get(0).size());
    assertEquals(4, components.get(0).get(0).length);
    assertEquals(4, components.get(0).get(1).length);
  }

  /**
   * Tests S2BuilderGraph.getDirectedComponents() with an input graph that contains degenerate
   * edges.
   */
  @Test
  public void testGetDirectedComponentsDegenerateEdges() {
    GraphClone graphClone = new GraphClone();
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    GraphOptions graphOptions =
        new GraphOptions(
            EdgeType.DIRECTED,
            DegenerateEdges.DISCARD_EXCESS,
            DuplicateEdges.KEEP,
            SiblingPairs.CREATE);

    builder.startLayer(new GraphCloningLayer(graphOptions, graphClone));
    builder.addShape(makeLaxPolylineOrDie("1:1, 1:1")); // A degenerate edge
    builder.addShape(makeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0")); // A-B-C-D-A square

    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));

    S2BuilderGraph graph = graphClone.graph();
    ArrayList<DirectedComponent> components = new ArrayList<>();
    assertTrue(
        error.toString(),
        graph.getDirectedComponents(DegenerateBoundaries.KEEP, components, error));

    assertEquals(2, components.size());
    assertEquals(1, components.get(0).size());
    assertEquals(1, components.get(0).get(0).length);
    assertEquals(2, components.get(1).size());
    assertEquals(4, components.get(1).get(0).length);
    assertEquals(4, components.get(1).get(1).length);
  }

  @Test
  public void testBasicVertexOutApi() {
    GraphClone graphClone = new GraphClone();
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    GraphOptions graphOptions =
        new GraphOptions(
            EdgeType.DIRECTED,
            DegenerateEdges.DISCARD_EXCESS,
            DuplicateEdges.KEEP,
            SiblingPairs.KEEP);

    builder.startLayer(new GraphCloningLayer(graphOptions, graphClone));
    builder.addShape(makeLaxPolylineOrDie("0:0, 1:1, 2:0"));
    builder.addShape(makeLaxPolylineOrDie("0:0, 1:0, 2:0"));
    builder.addShape(makeLaxPolylineOrDie("0:0, 0:1, 0:2"));

    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));

    S2BuilderGraph graph = graphClone.graph();
    VertexOutMap outmap = new VertexOutMap(graphClone.graph());
    Edge edge = new Edge();

    // There should be three edges out of the first vertex, 0:0.
    VertexOutEdges edgesFrom0 = outmap.edges(0);
    assertEquals(3, edgesFrom0.size());

    // Iterate over and validate the out-edges from 0:0.
    HashSet<Integer> dstIds = new HashSet<>();
    for (int i = 0; i < edgesFrom0.size(); ++i) {
      assertEquals(0, edgesFrom0.getSrcId(i));
      int dstId = edgesFrom0.getDstId(i);
      dstIds.add(dstId);
      edgesFrom0.get(i, edge);
      assertEquals(0, edge.srcId);
      assertEquals(dstId, edge.dstId);
      edgesFrom0.getReverse(i, edge);
      assertEquals(0, edge.dstId);
      assertEquals(dstId, edge.srcId);
    }
    assertEquals(3, dstIds.size()); // out-edges are unique, no duplicates in this test.
    dstIds.clear();

    // Check that out edges from each vertex actually start at the vertex.
    for (int vertexId = 0; vertexId < graph.numVertices(); ++vertexId) {
      // Get outgoing edges from vertex 'vertexId'
      VertexOutEdges outEdges = outmap.edges(vertexId);
      for (int e = 0;  e < outEdges.size(); ++e) {
        outEdges.get(0, edge);
        assertEquals(vertexId, edge.srcId);
      }
    }

    // There should be three edge ids out of the first vertex.
    VertexOutEdgeIds edgeIdsFrom0 = outmap.edgeIds(0);
    assertEquals(3, edgeIdsFrom0.size());

    // Iterate the out edge ids and validate the edges.
    OfInt edgeIdIter = edgeIdsFrom0.intIterator();
    while (edgeIdIter.hasNext()) {
      int edgeId = edgeIdIter.nextInt();
      assertEquals(0, graph.edgeSrcId(edgeId));
      dstIds.add(graph.edgeDstId(edgeId));
    }
    assertEquals(3, dstIds.size());
    dstIds.clear();

    // There should be one edge from the first vertex to the second.
    VertexOutEdgeIds edgeIdsFrom0To1 = outmap.edgeIds(0, 1);
    assertEquals(1, edgeIdsFrom0To1.size());

    // There should be no edge from the second vertex back to the first.
    VertexOutEdgeIds edgeIdsFrom1To0 = outmap.edgeIds(1, 0);
    assertEquals(0, edgeIdsFrom1To0.size());
  }

  @Test
  public void testInEdgeIdsSortOrder() {
    GraphClone graphClone = new GraphClone();
    S2Builder.Builder options  = new S2Builder.Builder();
    S2Builder builder = options.build();
    GraphOptions graphOptions =
        new GraphOptions(
            EdgeType.UNDIRECTED,
            DegenerateEdges.DISCARD,
            DuplicateEdges.KEEP,
            SiblingPairs.KEEP);

    builder.startLayer(new GraphCloningLayer(graphOptions, graphClone));
    builder.addPolyline(makePolylineOrDie("0:0, 1:0, 2:0, 3:0, 2:0, 1:0, 2:0, 3:0"));

    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));
    S2BuilderGraph graph = graphClone.graph();
    assertTrue(graph.getInEdgeIds().isEqualTo(
                   IntVector.of(1, 0, 5, 6, 7, 2, 3, 4, 11, 12, 13, 8, 9, 10)));
  }

  @Test
  public void testGetUndirectedComponentsDegenerateEdges() {
    GraphClone graphClone = new GraphClone();
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    GraphOptions graphOptions =
        new GraphOptions(
            EdgeType.UNDIRECTED,
            DegenerateEdges.DISCARD_EXCESS,
            DuplicateEdges.KEEP,
            SiblingPairs.DISCARD_EXCESS);

    builder.startLayer(new GraphCloningLayer(graphOptions, graphClone));
    builder.addShape(makeLaxPolylineOrDie("1:1, 1:1")); // A degenerate edge
    builder.addShape(makeLaxPolylineOrDie("0:0, 0:2, 2:2, 2:0, 0:0")); // A-B-C-D-A square

    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));

    S2BuilderGraph graph = graphClone.graph();
    ArrayList<UndirectedComponent> components = new ArrayList<>();

    assertTrue(
        error.toString(),
        graph.getUndirectedComponents(LoopType.CIRCUIT, components, error));

    // The result consists of two components, each with two complements. Each complement in this
    // example has exactly one loop. The loops in both complements of the first component have 1
    // vertex, while the loops in both complements of the second component have 4 vertices.
    assertEquals(2, components.size(), 2);

    // First component is the degenerate edge.
    assertEquals(1, components.get(0).getComplement(0).size());  // Complement 0 has one loop.
    assertEquals(1, components.get(0).getComplement(0).get(0).length); // Loop has one vertex.
    assertEquals(1, components.get(0).getComplement(1).size());
    assertEquals(1, components.get(0).getComplement(1).get(0).length);

    // Second component is the ABCDA square.
    assertEquals(1, components.get(1).getComplement(0).size());
    assertEquals(4, components.get(1).getComplement(0).get(0).length);
    assertEquals(1, components.get(1).getComplement(1).size());
    assertEquals(4, components.get(1).getComplement(1).get(0).length);
  }

  @Test
  public void testGetPolylinesUndirectedDegeneratePaths() {
    GraphClone graphClone = new GraphClone();
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    GraphOptions graphOptions =
        new GraphOptions(
            EdgeType.UNDIRECTED,
            DegenerateEdges.KEEP,
            DuplicateEdges.KEEP,
            SiblingPairs.KEEP);

    builder.startLayer(new GraphCloningLayer(graphOptions, graphClone));
    builder.addShape(makeLaxPolylineOrDie("1:1, 1:1"));
    builder.addShape(makeLaxPolylineOrDie("0:0, 0:0, 0:1, 0:1, 0:2, 0:2"));
    builder.addShape(makeLaxPolylineOrDie("1:1, 1:1"));

    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));

    S2BuilderGraph graph = graphClone.graph();
    PolylineBuilder polylineBuilder = new PolylineBuilder();
    polylineBuilder.init(graph);
    List<int[]> polylines = polylineBuilder.buildPaths();
    assertEquals(7, polylines.size());
  }

  @Test
  public void testGetPolylinesUndirectedDegenerateWalks() {
    GraphClone graphClone = new GraphClone();
    S2Builder.Builder options = new S2Builder.Builder();
    S2Builder builder = options.build();
    GraphOptions graphOptions =
        new GraphOptions(
            EdgeType.UNDIRECTED,
            DegenerateEdges.KEEP,
            DuplicateEdges.KEEP,
            SiblingPairs.KEEP);

    builder.startLayer(new GraphCloningLayer(graphOptions, graphClone));
    builder.addShape(makeLaxPolylineOrDie("1:1, 1:1"));
    builder.addShape(makeLaxPolylineOrDie("0:0, 0:0, 0:1, 0:1, 0:2, 0:2"));
    builder.addShape(makeLaxPolylineOrDie("1:1, 1:1"));

    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));

    S2BuilderGraph graph = graphClone.graph();
    PolylineBuilder polylineBuilder = new PolylineBuilder();
    polylineBuilder.init(graph);
    List<int[]> polylines = polylineBuilder.buildWalks();
    assertEquals(2, polylines.size());
    assertEquals(2, polylines.get(0).length);
    assertEquals(5, polylines.get(1).length);
  }

  /** A container that defines an edge for testing with testProcessEdges below. */
  private static class TestEdge {
    final int srcId;
    final int dstId;
    final IntVector inputEdgeIds;

    TestEdge(int srcId, int dstId, IntVector inputEdgeIds) {
      this.srcId = srcId;
      this.dstId = dstId;
      this.inputEdgeIds = inputEdgeIds;
    }

    TestEdge(int srcId, int dstId) {
      this.srcId = srcId;
      this.dstId = dstId;
      this.inputEdgeIds = new IntVector();
    }
  }

  private void testProcessEdges(
      List<TestEdge> input, List<TestEdge> expected, GraphOptions options) {
    testProcessEdges(input, expected, options, S2Error.Code.NO_ERROR);
  }

  private void testProcessEdges(
      List<TestEdge> inputEdges,
      List<TestEdge> expectedEdges,
      GraphOptions options,
      S2Error.Code expectedCode) {
    // A list of Edges to be processed by the S2BuilderGraph.
    EdgeList edges = new EdgeList();
    // Each edge in the graph has a set of input edge ids that were mapped to it. Those sets are
    // stored in the IdSetLexicon, and the map from edge id to IdSet id is stored in inputIdSetIds.
    IntVector inputIdSetIds = IntVector.empty();
    IdSetLexicon idSetLexicon = new IdSetLexicon();

    for (TestEdge inputEdge : inputEdges) {
      // Add each input TestEdge to the EdgeList to be processed.
      edges.add(inputEdge.srcId, inputEdge.dstId);
      // Add each TestEdge's set of input edge ids to the lexicon, and map the edge id to those.
      inputIdSetIds.add(idSetLexicon.add(inputEdge.inputEdgeIds));
    }

    // Process the edges with the given options.
    S2Error error = new S2Error();
    S2BuilderGraph.processEdges(options, edges, inputIdSetIds, idSetLexicon, error);

    // Verify the expected error code.
    assertEquals(expectedCode, error.code());

    // Verify the processed edges.
    assertEquals(edges.size(), inputIdSetIds.size());
    for (int edgeId = 0; edgeId < expectedEdges.size(); ++edgeId) {
      TestEdge expectedEdge = expectedEdges.get(edgeId);
      assertLessThan(edgeId, edges.size());
      assertEquals("(edge " + edgeId + ") src ", expectedEdge.srcId, edges.getSrcId(edgeId));
      assertEquals("(edge " + edgeId + ") dst ", expectedEdge.dstId, edges.getDstId(edgeId));
      // Get the processed set of input edges that ended up assigned to the edge
      IdSetLexicon.IdSet actualInputEdgeIdSet = idSetLexicon.idSet(inputIdSetIds.get(edgeId));
      // Check that the IntVector and IdSet have the same input edge ids.
      assertEquals(
          "(edge " + edgeId + ") in ",
          ImmutableIntSequence.viewOf(expectedEdge.inputEdgeIds),
          actualInputEdgeIdSet);
    }

    assertEquals("Too many output edges", expectedEdges.size(), edges.size());
  }

  private TestEdge te(int srcId, int dstId) {
    return new TestEdge(srcId, dstId);
  }
  private TestEdge te(int srcId, int dstId, IntVector inputEdgeIds) {
    return new TestEdge(srcId, dstId, inputEdgeIds);
  }
  private IntVector inIds(int... ids) {
    return IntVector.of(ids);
  }
  private ImmutableList<TestEdge> testEdges(TestEdge... edges) {
    return ImmutableList.copyOf(edges);
  }

  /**
   * A graph with two degenerate edges, and no input edge ids. Both should be discarded, resulting
   * in an empty graph.
   */
  @Test
  public void testProcessEdgesDiscardDegenerateEdges() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP,
        SiblingPairs.KEEP);
    testProcessEdges(
        testEdges(te(0, 0), te(0, 0)),
        testEdges(),
        options);
  }

  /**
   * A graph with two degenerate edges, and no input edge ids. Both should be kept, resulting in an
   * unchanged graph.
   */
  @Test
  public void testProcessEdgesKeepDuplicateDegenerateEdges() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.KEEP,
        DuplicateEdges.KEEP,
        SiblingPairs.KEEP);
    testProcessEdges(
        testEdges(te(0, 0), te(0, 0)),
        testEdges(te(0, 0), te(0, 0)),
        options);
  }

  /**
   * Two degenerate edges, each having one different input edge id. The edges should be merged,
   * resulting in a graph with one edge having both input edge ids.
   */
  @Test
  public void testProcessEdgesMergeDuplicateDegenerateEdges() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.KEEP,
        DuplicateEdges.MERGE,
        SiblingPairs.KEEP);
    testProcessEdges(
        testEdges(te(0, 0, inIds(1)), te(0, 0, inIds(2))),
        testEdges(te(0, 0, inIds(1, 2))),
        options);
  }

  /**
   * Edge count should be reduced to 2 (i.e., one undirected edge), and all labels should be merged.
   */
  @Test
  public void testProcessEdgesMergeUndirectedDuplicateDegenerateEdges() {
    GraphOptions options = new GraphOptions(
        EdgeType.UNDIRECTED,
        DegenerateEdges.KEEP,
        DuplicateEdges.MERGE,
        SiblingPairs.KEEP);
    testProcessEdges(
        testEdges(te(0, 0, inIds(1)), te(0, 0), te(0, 0), te(0, 0, inIds(2))),
        testEdges(te(0, 0, inIds(1, 2)), te(0, 0, inIds(1, 2))),
        options);
  }

  /**
   * Converting from UNDIRECTED to DIRECTED cuts the edge count in half and merges any edge labels.
   */
  @Test
  public void testProcessEdgesConvertedUndirectedDegenerateEdges() {
    GraphOptions options = new GraphOptions(
        EdgeType.UNDIRECTED,
        DegenerateEdges.KEEP,
        DuplicateEdges.KEEP,
        SiblingPairs.REQUIRE);
    testProcessEdges(
        testEdges(te(0, 0, inIds(1)), te(0, 0), te(0, 0), te(0, 0, inIds(2))),
        testEdges(te(0, 0, inIds(1, 2)), te(0, 0, inIds(1, 2))),
        options);
    assertEquals(EdgeType.DIRECTED, options.edgeType());
  }

  /** Like the test above, except that we also merge duplicates. */
  @Test
  public void testProcessEdgesMergeConvertedUndirectedDuplicateDegenerateEdges() {
    GraphOptions options = new GraphOptions(
        EdgeType.UNDIRECTED,
        DegenerateEdges.KEEP,
        DuplicateEdges.MERGE,
        SiblingPairs.REQUIRE);
    testProcessEdges(
        testEdges(te(0, 0, inIds(1)), te(0, 0), te(0, 0), te(0, 0, inIds(2))),
        testEdges(te(0, 0, inIds(1, 2))),
        options);
    assertEquals(EdgeType.DIRECTED, options.edgeType());
  }

  /**
   * Test that degenerate edges are discarded if they are connected to any non-degenerate edges
   * (whether they are incoming or outgoing, and whether they are lexicographically before or after
   * the degenerate edge).
   */
  @Test
  public void testProcessEdgesDiscardExcessConnectedDegenerateEdges() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD_EXCESS,
        DuplicateEdges.KEEP,
        SiblingPairs.KEEP);
    testProcessEdges(testEdges(te(0, 0), te(0, 1)), testEdges(te(0, 1)), options);
    testProcessEdges(testEdges(te(0, 0), te(1, 0)), testEdges(te(1, 0)), options);
    testProcessEdges(testEdges(te(0, 1), te(1, 1)), testEdges(te(0, 1)), options);
    testProcessEdges(testEdges(te(1, 0), te(1, 1)), testEdges(te(1, 0)), options);
  }

  // Test that DISCARD_EXCESS merges any duplicate degenerate edges together.
  @Test
  public void testProcessEdgesDiscardExcessIsolatedDegenerateEdges() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD_EXCESS,
        DuplicateEdges.KEEP,
        SiblingPairs.KEEP);
    testProcessEdges(
        testEdges(te(0, 0, inIds(1)), te(0, 0, inIds(2))),
        testEdges(te(0, 0, inIds(1, 2))),
        options);
  }

  /**
   * Test that DISCARD_EXCESS with SiblingPairs.REQUIRE merges any duplicate edges together and
   * converts the edges from UNDIRECTED to DIRECTED.
   */
  @Test
  public void testProcessEdgesDiscardExcessConvertedUndirectedIsolatedDegenerateEdges() {
    GraphOptions options = new GraphOptions(
        EdgeType.UNDIRECTED,
        DegenerateEdges.DISCARD_EXCESS,
        DuplicateEdges.KEEP,
        SiblingPairs.REQUIRE);
    testProcessEdges(
        testEdges(te(0, 0, inIds(1)), te(0, 0, inIds(2)), te(0, 0, inIds(3)), te(0, 0)),
        testEdges(te(0, 0, inIds(1, 2, 3))),
        options);
    assertEquals(EdgeType.DIRECTED, options.edgeType());
  }

  /**
   * Test that when SiblingPairs.DISCARD or SiblingPairs.DISCARD_EXCESS is specified, the edge
   * labels of degenerate edges are merged together (for consistency, since these options merge the
   * labels of all non-degenerate edges as well).
   */
  @Test
  public void testProcessEdgesSiblingPairsDiscardMergesDegenerateEdgeLabels() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.KEEP,
        DuplicateEdges.KEEP,
        SiblingPairs.DISCARD);
    testProcessEdges(
        testEdges(te(0, 0, inIds(1)), te(0, 0, inIds(2)), te(0, 0, inIds(3))),
        testEdges(te(0, 0, inIds(1, 2, 3)), te(0, 0, inIds(1, 2, 3)), te(0, 0, inIds(1, 2, 3))),
        options);
    options.setSiblingPairs(SiblingPairs.DISCARD_EXCESS);
    testProcessEdges(
        testEdges(te(0, 0, inIds(1)), te(0, 0, inIds(2)), te(0, 0, inIds(3))),
        testEdges(te(0, 0, inIds(1, 2, 3)), te(0, 0, inIds(1, 2, 3)), te(0, 0, inIds(1, 2, 3))),
        options);
  }

  /** Test that siblings are not discarded in a directed graph. */
  @Test
  public void testProcessEdgesKeepSiblingPairs() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP,
        SiblingPairs.KEEP);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
  }

  /** Test that duplicate edges are merged, even if provided in unsorted order. */
  @Test
  public void testProcessEdgesMergeUnsortedDuplicates() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.MERGE,
        SiblingPairs.KEEP);
    testProcessEdges(
        testEdges(te(0, 1), te(2, 0), te(0, 1), te(2, 1), te(2, 0), te(1, 0)),
        testEdges(te(0, 1), te(1, 0), te(2, 0), te(2, 1)),
        options);
  }

  /** Test that duplicate siblings are discarded, but one sibling is retained. */
  @Test
  public void testProcessEdgesMergeDuplicateSiblingPairs() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.MERGE,
        SiblingPairs.KEEP);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
  }

  /** Check that matched pairs are discarded, leaving behind any excess edges. */
  @Test
  public void testProcessEdgesDiscardSiblingPairs() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP,
        SiblingPairs.DISCARD);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0)),
        testEdges(),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(1, 0), te(1, 0)),
        testEdges(),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(0, 1), te(1, 0)),
        testEdges(te(0, 1), te(0, 1)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0), te(1, 0), te(1, 0)),
        testEdges(te(1, 0), te(1, 0)),
        options);
  }

  /** Check that matched pairs are discarded, and then any remaining edges are merged. */
  @Test
  public void testProcessEdgesDiscardSiblingPairsMergeDuplicates() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.MERGE,
        SiblingPairs.DISCARD);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(1, 0), te(1, 0)),
        testEdges(),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(0, 1), te(1, 0)),
        testEdges(te(0, 1)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0), te(1, 0), te(1, 0)),
        testEdges(te(1, 0)),
        options);
  }

  /**
   * An undirected sibling pair consists of four edges, two in each direction: see the {@link
   * S2Builder} class javadoc. Since undirected edges always come in pairs, this means that the
   * result always consists of either 0 or 2 edges.
   */
  @Test
  public void testProcessEdgesDiscardUndirectedSiblingPairs() {
    GraphOptions options = new GraphOptions(
        EdgeType.UNDIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP,
        SiblingPairs.DISCARD);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(1, 0), te(1, 0)),
        testEdges(),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(0, 1), te(1, 0), te(1, 0), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
  }

  /**
   * Like SiblingPairs.DISCARD, except that one sibling pair is kept if the result would otherwise
   * be empty.
   */
  @Test
  public void testProcessEdgesDiscardExcessSiblingPairs() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP,
        SiblingPairs.DISCARD_EXCESS);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(1, 0), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(0, 1), te(1, 0)),
        testEdges(te(0, 1), te(0, 1)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0), te(1, 0), te(1, 0)),
        testEdges(te(1, 0), te(1, 0)),
        options);
  }

  /**
   * Like SiblingPairs.DISCARD, except that one sibling pair is kept if the result would otherwise
   * be empty.
   */
  @Test
  public void testProcessEdgesDiscardExcessSiblingPairsMergeDuplicates() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.MERGE,
        SiblingPairs.DISCARD_EXCESS);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(1, 0), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(0, 1), te(1, 0)),
        testEdges(te(0, 1)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0), te(1, 0), te(1, 0)),
        testEdges(te(1, 0)),
        options);
  }

  /**
   * Like SiblingPairs.DISCARD, except that one undirected sibling pair (4 edges) is kept if the
   * result would otherwise be empty.
   */
  @Test
  public void testProcessEdgesDiscardExcessUndirectedSiblingPairs() {
    GraphOptions options = new GraphOptions(
        EdgeType.UNDIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP,
        SiblingPairs.DISCARD_EXCESS);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(1, 0), te(1, 0)),
        testEdges(te(0, 1), te(0, 1), te(1, 0), te(1, 0)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(0, 1), te(1, 0), te(1, 0), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
  }

  /** Verify that siblings are created when needed, including for duplicates. */
  @Test
  public void testProcessEdgesCreateSiblingPairs() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP,
        SiblingPairs.CREATE);
    testProcessEdges(
        testEdges(te(0, 1)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1)),
        testEdges(te(0, 1), te(0, 1), te(1, 0), te(1, 0)),
        options);
  }

  /** Like SiblingPairs.CREATE, but generates an error. */
  @Test
  public void testProcessEdgesRequireSiblingPairs() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP,
        SiblingPairs.REQUIRE);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    testProcessEdges(
        testEdges(te(0, 1)),
        testEdges(te(0, 1), te(1, 0)),
        options,
        S2Error.Code.BUILDER_MISSING_EXPECTED_SIBLING_EDGES);
  }

  /**
   * An undirected sibling pair consists of 4 edges, but SiblingPairs.CREATE also converts the graph
   * to EdgeType.DIRECTED and cuts the number of edges in half.
   */
  @Test
  public void testProcessEdgesCreateUndirectedSiblingPairs() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP,
        SiblingPairs.CREATE);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    assertEquals(EdgeType.DIRECTED, options.edgeType());

    options.setEdgeType(EdgeType.UNDIRECTED);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(1, 0), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    assertEquals(EdgeType.DIRECTED, options.edgeType());

    options.setEdgeType(EdgeType.UNDIRECTED);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(0, 1), te(1, 0), te(1, 0), te(1, 0)),
        testEdges(te(0, 1), te(0, 1), te(1, 0), te(1, 0)),
        options);
    assertEquals(EdgeType.DIRECTED, options.edgeType());
  }

  /** Verify that siblings are created when needed, while duplicate inputs are discarded. */
  @Test
  public void testProcessEdgesCreateSiblingPairsMergeDuplicates() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.MERGE,
        SiblingPairs.CREATE);
    testProcessEdges(
        testEdges(te(0, 1)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1)),
        testEdges(te(0, 1), te(1, 0)),
        options);
  }

  /** Tests that duplicates are merged but one siblings in each direction is retained. */
  @Test
  public void testProcessEdgesCreateUndirectedSiblingPairsMergeDuplicates() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.MERGE,
        SiblingPairs.CREATE);
    testProcessEdges(
        testEdges(te(0, 1), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    assertEquals(EdgeType.DIRECTED, options.edgeType());

    options.setEdgeType(EdgeType.UNDIRECTED);
    testProcessEdges(
        testEdges(te(0, 1), te(0, 1), te(0, 1), te(1, 0), te(1, 0), te(1, 0)),
        testEdges(te(0, 1), te(1, 0)),
        options);
    assertEquals(EdgeType.DIRECTED, options.edgeType());
  }

  private static Edge e(int src, int dst) {
    return Edge.of(src, dst);
  }

  private void checkMakeSubgraph(
      S2BuilderGraph graph,
      IdSetLexicon newInputEdgeIdSetLexicon,
      GraphOptions newOptions,
      EdgeList newEdges,
      IntVector newInputEdgeIdSetIds,
      GraphOptions expectedOptions,
      EdgeList expectedEdges,
      IntVector expectedInputEdgeIdSetIds) {
    S2Error error = new S2Error();
    S2BuilderGraph newGraph = graph.makeSubgraph(
        newOptions, newEdges, newInputEdgeIdSetIds,
        newInputEdgeIdSetLexicon, null, error);

    // The newGraph and original graph use exactly the same vertices, labelSetIds, and
    // labelSetLexicon.
    assertSame(newGraph.vertices(), graph.vertices());
    assertSame(newGraph.labelSetIds(), graph.labelSetIds());
    assertSame(newGraph.labelSetLexicon(), graph.labelSetLexicon());

    // The newGraph edges, inputEdgeIdSetIds, and inputEdgeIdSetLexicon should be stored in the
    // containers that were actually provided.
    assertSame(newGraph.edges(), newEdges);
    assertSame(newGraph.inputEdgeIdSetIds(), newInputEdgeIdSetIds);
    assertSame(newGraph.inputEdgeIdSetLexicon(), newInputEdgeIdSetLexicon);

    // The new graph should have the expected options.
    assertEquals(expectedOptions.edgeType(), newGraph.options().edgeType());
    assertEquals(expectedOptions.degenerateEdges(), newGraph.options().degenerateEdges());
    assertEquals(expectedOptions.duplicateEdges(), newGraph.options().duplicateEdges());
    assertEquals(expectedOptions.siblingPairs(), newGraph.options().siblingPairs());

    // The edges should be updated according to the requested options.
    assertTrue(expectedEdges.isEqualTo(newGraph.edges()));
    assertTrue(expectedInputEdgeIdSetIds.isEqualTo(newGraph.inputEdgeIdSetIds()));
    assertTrue(newInputEdgeIdSetLexicon.isEqualTo(newGraph.inputEdgeIdSetLexicon()));
  }

  /**
   * Test that makeSubgraph() doesn't transform edges into edge pairs when creating an undirected
   * subgraph of an undirected graph.
   */
  @Test
  public void testMakeSubgraphUndirectedToUndirected() {
    GraphOptions options = new GraphOptions(
        EdgeType.UNDIRECTED,
        DegenerateEdges.KEEP,
        DuplicateEdges.KEEP,
        SiblingPairs.KEEP);
    List<S2Point> vertices = S2TextFormat.parsePointsOrDie("0:0, 0:1, 1:1");
    EdgeList edges = EdgeList.of(e(0, 0), e(0, 0), e(1, 2), e(2, 1));
    IntVector inputEdgeIdSetIds = IntVector.of(0, 0, 1, 1);
    IntVector labelSetIds = new IntVector();
    IdSetLexicon inputEdgeIdSetLexicon = new IdSetLexicon();
    IdSetLexicon labelSetLexicon = new IdSetLexicon();

    S2BuilderGraph graph = new S2BuilderGraph(
      options, vertices, edges, inputEdgeIdSetIds,
      inputEdgeIdSetLexicon, labelSetIds, labelSetLexicon, null);

    // Now create a subgraph that also has undirected edges but discards degenerate edges.
    GraphOptions newOptions = new GraphOptions(
        EdgeType.UNDIRECTED,
        DegenerateEdges.DISCARD,
        DuplicateEdges.KEEP,
        SiblingPairs.KEEP);
    EdgeList expectedEdges = EdgeList.of(e(1, 2), e(2, 1));
    IntVector expectedInputEdgeIdSetIds = IntVector.of(1, 1);
    checkMakeSubgraph(
        graph, inputEdgeIdSetLexicon,
        newOptions, edges, inputEdgeIdSetIds,
        newOptions, expectedEdges, expectedInputEdgeIdSetIds);
  }

  /** Test transforming directed edges into undirected edges (which doubles the number of edges). */
  @Test
  public void testMakeSubgraphDirectedToUndirected() {
    GraphOptions options = new GraphOptions(
        EdgeType.DIRECTED,
        DegenerateEdges.KEEP,
        DuplicateEdges.KEEP,
        SiblingPairs.KEEP);
    List<S2Point> vertices = S2TextFormat.parsePointsOrDie("0:0, 0:1, 1:1");
    EdgeList edges = EdgeList.of(e(0, 0), e(0, 1), e(1, 2), e(1, 2), e(2, 1));
    IntVector inputEdgeIdSetIds = IntVector.of(1, 2, 3, 3, 3);
    IntVector labelSetIds = new IntVector();
    IdSetLexicon inputEdgeIdSetLexicon = new IdSetLexicon();
    IdSetLexicon labelSetLexicon = new IdSetLexicon();
    S2BuilderGraph graph = new S2BuilderGraph(
      options, vertices, edges, inputEdgeIdSetIds,
      inputEdgeIdSetLexicon, labelSetIds, labelSetLexicon, null);

    // Now create a subgraph with undirected edges and different options.
    GraphOptions newOptions = new GraphOptions(
        EdgeType.UNDIRECTED,
        DegenerateEdges.KEEP,
        DuplicateEdges.KEEP,
        SiblingPairs.DISCARD_EXCESS);
    EdgeList expectedEdges = EdgeList.of(
      e(0, 0), e(0, 0),  // Undirected degenerate edge.
      e(0, 1), e(1, 0),  // Undirected edge.
      e(1, 2), e(2, 1)   // Undirected edge after discarding sibling pair.
    );
    IntVector expectedInputEdgeIdSetIds = IntVector.of(
      1, 1, 2, IdSetLexicon.EMPTY_SET_ID, 3, 3);

    checkMakeSubgraph(
        graph, inputEdgeIdSetLexicon,
        newOptions, edges, inputEdgeIdSetIds,
        newOptions, expectedEdges, expectedInputEdgeIdSetIds);
  }
}
