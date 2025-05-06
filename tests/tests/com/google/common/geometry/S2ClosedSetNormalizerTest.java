/*
 * Copyright 2024 Google Inc.
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

import static java.util.stream.Collectors.groupingBy;
import static java.util.stream.Collectors.joining;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2BuilderGraph.PolylineType;
import com.google.common.geometry.S2BuilderUtil.IndexedLayer;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public final class S2ClosedSetNormalizerTest {

  private boolean suppressLowerDimensions;
  GraphOptions pointGraphOptions = new GraphOptions();
  GraphOptions lineGraphOptions = new GraphOptions();
  GraphOptions polygonGraphOptions = new GraphOptions();

  @Before
  public void setup() {
    suppressLowerDimensions = true;

    pointGraphOptions.setEdgeType(EdgeType.DIRECTED);
    pointGraphOptions.setDegenerateEdges(DegenerateEdges.KEEP);
    pointGraphOptions.setDuplicateEdges(DuplicateEdges.KEEP);
    pointGraphOptions.setSiblingPairs(SiblingPairs.KEEP);

    lineGraphOptions.setEdgeType(EdgeType.UNDIRECTED);
    lineGraphOptions.setDegenerateEdges(DegenerateEdges.KEEP);
    lineGraphOptions.setDuplicateEdges(DuplicateEdges.KEEP);
    lineGraphOptions.setSiblingPairs(SiblingPairs.KEEP);

    polygonGraphOptions.setEdgeType(EdgeType.DIRECTED);
    polygonGraphOptions.setDegenerateEdges(DegenerateEdges.KEEP);
    polygonGraphOptions.setDuplicateEdges(DuplicateEdges.KEEP);
    polygonGraphOptions.setSiblingPairs(SiblingPairs.KEEP);
  }

  private void run(String inputStr, String expectedStr) {
    S2Builder builder = new S2Builder.Builder().build();

    Map<Integer, List<S2Shape>> expectedShapesByDim =
        S2TextFormat.makeIndexOrDie(expectedStr).getShapes().stream()
            .collect(groupingBy(S2Shape::dimension));
    CopyGraphLayer expectedPointLayer =
        addCopyGraphLayer(expectedShapesByDim.get(0), pointGraphOptions, builder);
    CopyGraphLayer expectedLineLayer =
        addCopyGraphLayer(expectedShapesByDim.get(1), lineGraphOptions, builder);
    CopyGraphLayer expectedPolygonLayer =
        addCopyGraphLayer(expectedShapesByDim.get(2), polygonGraphOptions, builder);

    CopyGraphLayer actualPointLayer = new CopyGraphLayer(pointGraphOptions);
    CopyGraphLayer actualLineLayer = new CopyGraphLayer(lineGraphOptions);
    CopyGraphLayer actualPolygonLayer = new CopyGraphLayer(polygonGraphOptions);

    S2ClosedSetNormalizer.Options options = new S2ClosedSetNormalizer.Options();
    options.setSuppressLowerDimensions(suppressLowerDimensions);
    S2ClosedSetNormalizer normalizer =
        new S2ClosedSetNormalizer(options, actualPointLayer, actualLineLayer, actualPolygonLayer);

    Map<Integer, List<S2Shape>> inputShapesByDim =
        S2TextFormat.makeIndexOrDie(inputStr).getShapes().stream()
            .collect(groupingBy(S2Shape::dimension));
    addLayer(inputShapesByDim.get(0), normalizer.pointLayer(), builder);
    addLayer(inputShapesByDim.get(1), normalizer.lineLayer(), builder);
    addLayer(inputShapesByDim.get(2), normalizer.polygonLayer(), builder);

    S2Error error = new S2Error();
    assertTrue(error.toString(), builder.build(error));

    assertGraphEquals(expectedPointLayer.graph, actualPointLayer.graph);
    assertGraphEquals(expectedLineLayer.graph, actualLineLayer.graph);
    assertGraphEquals(expectedPolygonLayer.graph, actualPolygonLayer.graph);
  }

  private static CopyGraphLayer addCopyGraphLayer(
      List<S2Shape> shapes, GraphOptions graphOptions, S2Builder builder) {
    CopyGraphLayer layer = new CopyGraphLayer(graphOptions);
    addLayer(shapes, layer, builder);
    return layer;
  }

  private static void addLayer(List<S2Shape> shapes, S2BuilderLayer layer, S2Builder builder) {
    builder.startLayer(layer);
    if (shapes != null) {
      shapes.forEach(builder::addShape);
    }
  }

  private static void assertGraphEquals(S2BuilderGraph expected, S2BuilderGraph actual) {
    String expectedStr = graphToString(expected);
    String actualStr = graphToString(actual);
    assertEquals(expectedStr, actualStr);
  }

  private static String graphToString(S2BuilderGraph graph) {
    return IntStream.range(0, graph.numEdges())
        .mapToObj(i -> edgeToString(graph, i))
        .collect(joining("; "));
  }

  private static String edgeToString(S2BuilderGraph graph, int i) {
    return Stream.of(
            S2TextFormat.toString(graph.vertex(graph.edgeSrcId(i))),
            S2TextFormat.toString(graph.vertex(graph.edgeDstId(i))))
        .collect(joining("->"));
  }

  @Test
  public void emptyGraphs_returnsEmptyGraphs() {
    run("# #", "# #");
  }

  @Test
  public void nonDegenerateInputs_unchanged() {
    run("0:0 # 1:0, 1:1 | 1:2, 1:3 # 2:2, 2:3, 3:2", "0:0 # 1:0, 1:1 | 1:2, 1:3 # 2:2, 2:3, 3:2");
  }

  @Test
  public void pointShell_demotedToPoint() {
    run("# # 0:0", "0:0 # #");
  }

  @Test
  public void pointHole_disappears() {
    run("# # 0:0, 0:3, 3:0 | 1:1", "# # 0:0, 0:3, 3:0");
  }

  @Test
  public void pointPolyline_demotedToPoint() {
    // Verify that a single degenerate polyline edge is transformed into a single point. Note that
    // since the polyline layer is undirected while the point layer is not, this tests that the edge
    // count is halved when the edge is demoted.
    run("# 0:0, 0:0 #", "0:0 # #");
  }

  @Test
  public void siblingPairShell_demotedToLine() {
    run("# # 0:0, 1:0 ", "# 0:0, 1:0 #");
  }

  @Test
  public void siblingPairHole_disappears() {
    run("# # 0:0, 0:3, 3:0; 0:0, 1:1", "# # 0:0, 0:3, 3:0");
  }

  @Test
  public void pointAtPolygonVertex_suppressed() {
    run("0:0 | 0:1 | 1:0 # # 0:0, 0:1, 1:0", "# # 0:0, 0:1, 1:0");
  }

  @Test
  public void pointAtPolygonVertex_unsuppressed() {
    suppressLowerDimensions = false;
    run("0:0 | 0:1 | 1:0 # # 0:0, 0:1, 1:0", "0:0 | 0:1 | 1:0 # # 0:0, 0:1, 1:0");
  }

  @Test
  public void pointAtPolylineVertex_suppressed() {
    run("0:0 | 0:1 # 0:0, 0:1 #", "# 0:0, 0:1 #");
  }

  @Test
  public void pointAtPolylineVertex_unsuppressed() {
    suppressLowerDimensions = false;
    run("0:0 | 0:1 # 0:0, 0:1 #", "0:0 | 0:1 # 0:0, 0:1 #");
  }

  @Test
  public void pointShell_suppressedByPolylineEdge() {
    // This tests the case where a single-point shell is demoted to a point which is then suppressed
    // by a matching polyline vertex.
    run("# 0:0, 1:0 # 0:0; 1:0", "# 0:0, 1:0 #");
  }

  @Test
  public void pointShell_unsuppressedByPolylineEdge() {
    suppressLowerDimensions = false;
    run("# 0:0, 1:0 # 0:0; 1:0", "0:0 | 1:0 # 0:0, 1:0 #");
  }

  @Test
  public void polylineEdge_suppressedByPolygonEdge() {
    run("# 0:0, 0:1 # 0:0, 0:1, 1:0", "# # 0:0, 0:1, 1:0");
  }

  @Test
  public void polylineEdge_unsuppressedByPolygonEdge() {
    suppressLowerDimensions = false;
    run("# 0:0, 0:1 # 0:0, 0:1, 1:0", "# 0:0, 0:1 # 0:0, 0:1, 1:0");
  }

  @Test
  public void polylineEdge_suppressedByReversePolygonEdge() {
    lineGraphOptions.setEdgeType(EdgeType.DIRECTED);
    run("# 1:0, 0:0 # 0:0, 0:1, 1:0", "# # 0:0, 0:1, 1:0");
  }

  @Test
  public void polylineEdge_unsuppressedByReversePolygonEdge() {
    suppressLowerDimensions = false;
    lineGraphOptions.setEdgeType(EdgeType.DIRECTED);
    run("# 1:0, 0:0 # 0:0, 0:1, 1:0", "# 1:0, 0:0 # 0:0, 0:1, 1:0");
  }

  @Test
  public void duplicateEdgeMergingNotRequested_edgesNotMerged() {
    // Verify that when DuplicateEdges::KEEP is specified, demoted edges are added as new edges
    // rather than being merged with existing ones. (Note that NormalizeTest specifies
    // DuplicateEdges::KEEP by default.)
    run(
        "0:0 | 0:0 # 0:0, 0:0 | 0:1, 0:2 # 0:0; 0:1, 0:2",
        "0:0 | 0:0 | 0:0 | 0:0 # 0:1, 0:2 | 0:1, 0:2 #");
  }

  @Test
  public void duplicateEdgeMergingRequested_edgesMerged() {
    pointGraphOptions.setDuplicateEdges(DuplicateEdges.MERGE);
    lineGraphOptions.setDuplicateEdges(DuplicateEdges.MERGE);
    run("0:0 | 0:0 # 0:0, 0:0 | 0:1, 0:2 # 0:0; 0:1, 0:2", "0:0 # 0:1, 0:2 #");
  }

  /**
   * The example from the S2ClosedSetNormalizer Javadoc. If this code changes, please update the
   * Javadoc to match.
   */
  static boolean computeUnion(S2ShapeIndex a, S2ShapeIndex b, S2ShapeIndex output, S2Error error) {
    var polylineOptions =
        new S2PolylineVectorLayer.Options()
            .setEdgeType(EdgeType.UNDIRECTED)
            .setPolylineType(PolylineType.WALK)
            .setDuplicateEdges(DuplicateEdges.MERGE);

    S2ClosedSetNormalizer normalizer =
        new S2ClosedSetNormalizer(
            new IndexedLayer<>(output, new S2PointVectorLayer()),
            new IndexedLayer<>(output, new S2PolylineVectorLayer(polylineOptions)),
            new IndexedLayer<>(output, new S2PolygonLayer()));

    S2BooleanOperation op =
        new S2BooleanOperation.Builder()
            .build(
                S2BooleanOperation.OpType.UNION,
                normalizer.pointLayer(),
                normalizer.lineLayer(),
                normalizer.polygonLayer());
    return op.build(a, b, error);
  }

  @Test
  public void computeUnionWithMixedGeometry_computesNormalizedUnion() {
    // Verifies that the code above works.  Features tested include:
    //  - Points and polylines in the interior of the other polygon are removed
    //  - Degenerate polygon shells are converted to points/polylines
    //  - Degenerate polygon holes are removed
    //  - Points coincident with polyline or polygon edges are removed
    //  - Polyline edges coincident with polygon edges are removed
    S2ShapeIndex a =
        S2TextFormat.makeIndexOrDie(
            "0:0 | 10:10 | 20:20 # 0:0, 0:10 | 0:0, 10:0 | 15:15, 16:16 # 0:0, 0:10, 10:10, 10:0;"
                + " 0:0, 1:1; 2:2; 10:10, 11:11; 12:12");
    S2ShapeIndex b =
        S2TextFormat.makeIndexOrDie(
            "0:10 | 10:0 | 3:3 | 16:16 # 10:10, 0:10 | 10:10, 10:0 | 5:5, 6:6 # 19:19, 19:21,"
                + " 21:21, 21:19");
    S2ShapeIndex result = new S2ShapeIndex();
    S2Error error = new S2Error();
    assertTrue(computeUnion(a, b, result, error));
    assertEquals(
        "12:12 # 15:15, 16:16 | 10:10, 11:11 # 19:19, 19:21, 21:21, 21:19; 0:0, 0:10, 10:10, 10:0",
        S2TextFormat.toString(result));
  }

  private static class CopyGraphLayer implements S2BuilderLayer {
    private S2BuilderGraph graph;

    private final GraphOptions graphOptions;

    CopyGraphLayer(GraphOptions graphOptions) {
      this.graphOptions = graphOptions;
    }

    @Override
    public GraphOptions graphOptions() {
      return graphOptions;
    }

    @Override
    public boolean build(S2BuilderGraph graph, S2Error error) {
      this.graph = S2ClosedSetNormalizer.copyGraphWithOptions(graph, graph.options());
      return true;
    }
  }
}
