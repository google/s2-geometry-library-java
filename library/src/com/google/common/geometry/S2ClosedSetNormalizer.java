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

import com.google.common.annotations.VisibleForTesting;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2BuilderGraph.Edge;
import com.google.common.geometry.S2BuilderGraph.EdgeList;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * The purpose of this class is to allow S2Builder::Layer implementations to remove polygon and
 * polyline degeneracies by converting them to polylines or points.
 *
 * <p>A polyline degeneracy is a polyline consisting of a single degenerate edge. A polygon
 * degeneracy is either a single-vertex loop (a degenerate edge from a vertex to itself) or a
 * sibling edge pair (consisting of an edge and its corresponding reverse edge). Polygon
 * degeneracies are further classified as shells or holes depending on whether they are located in
 * the exterior or interior of the polygon respectively. For example, a single-vertex loop contained
 * within a polygon shell would be classified as a hole.
 *
 * <p>All objects are modeled as closed, i.e. polygons contain their boundaries and polylines
 * contain their endpoints. Note that under this model, degenerate polygon shells and holes need to
 * be handled differently. Degenerate shells are converted to polylines or points, whereas
 * degenerate holes do not affect the set of points contained by the polygon and are simply
 * discarded.
 *
 * <p>Specifically, given three S2Builder::Graphs (corresponding to points, polylines, and
 * polygons), this class makes the following transformations:
 *
 * <ul>
 *   <li>Polygon sibling edge pairs are either discarded (for holes) or converted to a pair of
 *       polyline edges (for shells).
 *   <li>Degenerate polygon edges are either discarded (for holes) or converted to points (for
 *       shells).
 *   <li>Degenerate polyline edges are converted to points.
 * </ul>
 *
 * <p>Optionally, this class further will normalize the graphs by suppressing edges that are
 * duplicates of higher-dimensional edges. In other words:
 *
 * <ul>
 *   <li>Polyline edges that coincide with polygon edges are discarded.
 *   <li>Points that coincide with polyline or polygon vertices are discarded.
 * </ul>
 *
 * <p>(When edges are discarded, any labels attached to those edges are discarded as well.)
 *
 * <p>This class takes three graphs as input and yields three graphs as output. However note that
 * the output graphs are *not* independent objects; they may point to data in the input graphs or
 * data owned by the ClosedSetNormalizer itself. For this reason the input graphs and
 * ClosedSetNormalizer must persist until the output graphs are no longer needed.
 *
 * <p>Finally, note that although this class may be necessary in some situations (e.g., to implement
 * the OGC Simple Features Access spec), in general the recommended approach to degeneracies is
 * simply to keep them (by using a representation such as S2LaxPolygonShape or S2LaxPolylineShape).
 * Keeping degeneracies has many advantages, such as not needing to deal with geometry of multiple
 * dimensions, and being able to preserve polygon boundaries accurately (including degenerate
 * holes).
 *
 * <p>Below is an example showing how to compute the union of two {@link S2ShapeIndex}es containing
 * points, lines and/or polygons, and save the result as a collection of {@link S2Point}s, {@link
 * S2Polyline}s and {@link S2Polygon}s in another S2ShapeIndex (where degeneracies have been
 * normalized to objects of lower dimension). Note that in this example, IndexedLayer is used so
 * that the layers passed to S2ClosedSetNormalizer automatically add their shapes to the output
 * S2ShapeIndex. For other types of output, one may need to save reference to the layers passed to
 * S2ClosedSetNormalizer to retrieve their output.
 *
 * {@snippet :
 * static boolean computeUnion(S2ShapeIndex a, S2ShapeIndex b, S2ShapeIndex output, S2Error error) {
 *   var polylineOptions =
 *       new S2PolylineVectorLayer.Options()
 *           .setEdgeType(EdgeType.UNDIRECTED)
 *           .setPolylineType(PolylineType.WALK)
 *           .setDuplicateEdges(DuplicateEdges.MERGE);
 *
 *   S2ClosedSetNormalizer normalizer =
 *       new S2ClosedSetNormalizer(
 *           new IndexedLayer<>(output, new S2PointVectorLayer()),
 *           new IndexedLayer<>(output, new S2PolylineVectorLayer(polylineOptions)),
 *           new IndexedLayer<>(output, new S2PolygonLayer()));
 *
 *   S2BooleanOperation op =
 *       new S2BooleanOperation.Builder()
 *           .build(
 *               S2BooleanOperation.OpType.UNION,
 *               normalizer.pointLayer(),
 *               normalizer.lineLayer(),
 *               normalizer.polygonLayer());
 *   return op.build(a, b, error);
 * }
 * }
 */
public final class S2ClosedSetNormalizer {

  /**
   * Options for the S2ClosedSetNormalizer:
   *
   * <ul>
   *   <li>suppressLowerDimensions: If true, the graphs are further normalized by discarding
   *       lower-dimensional edges that coincide with higher-dimensional edges. True by default.
   * </ul>
   */
  public static class Options {
    private boolean suppressLowerDimensions;

    public Options(boolean suppressLowerDimensions) {
      this.suppressLowerDimensions = suppressLowerDimensions;
    }

    public Options() {
      this(true);
    }

    @CanIgnoreReturnValue
    public Options setSuppressLowerDimensions(boolean suppressLowerDimensions) {
      this.suppressLowerDimensions = suppressLowerDimensions;
      return this;
    }

    public boolean suppressLowerDimensions() {
      return suppressLowerDimensions;
    }
  }

  private final Options options;

  /** Intermediate data needed to normalize one dimension of the input. */
  private class NormalizedLayer {
    public GraphOptions inputGraphOptions = null;
    public GraphOptions outputGraphOptions = null;
    public final EdgeList newEdges = new EdgeList();
    public final IntVector newInputEdgeIds = new IntVector();
    public IdSetLexicon newInputEdgeIdSetLexicon = null;
    public S2BuilderGraph inputGraph = null;
    public final S2BuilderLayer inputLayer = new DimensionLayer(this);
    public S2BuilderGraph outputGraph = null;
    public S2BuilderLayer outputLayer = null;
  }

  private final NormalizedLayer normalizedPointLayer = new NormalizedLayer();
  private final NormalizedLayer normalizedLineLayer = new NormalizedLayer();
  private final NormalizedLayer normalizedPolygonLayer = new NormalizedLayer();

  /**
   * A vector of incoming polygon edges sorted in lexicographic order. This is used to suppress
   * directed polyline edges that match a polygon edge in the reverse direction.
   */
  private IntVector polygonInputEdgeIds = null;

  /**
   * excludeFromPointOutput[i] is true if vertex[i] belongs to a non-degenerate edge, and therefore
   * should be suppressed from the output graph for points.
   */
  private boolean[] excludeFromPointOutput = null;

  /**
   * Constructs a ClosedSetNormalizer that will pass a normalized S2BuilderGraph to each of the
   * three given output layers.
   *
   * <ul>
   *   <li>REQUIRES: pointLayer.graphOptions().edgeType() == DIRECTED
   *   <li>REQUIRES: lineLayer.graphOptions().siblingPairs() != {CREATE, REQUIRE}
   *   <li>REQUIRES: polygonLayer.graphOptions().edgeType() == DIRECTED
   * </ul>
   */
  public S2ClosedSetNormalizer(
      Options options,
      S2BuilderLayer pointLayer,
      S2BuilderLayer lineLayer,
      S2BuilderLayer polygonLayer) {
    normalizedPointLayer.outputLayer = pointLayer;
    normalizedLineLayer.outputLayer = lineLayer;
    normalizedPolygonLayer.outputLayer = polygonLayer;

    GraphOptions pointGraphOptions = pointLayer.graphOptions();
    GraphOptions lineGraphOptions = lineLayer.graphOptions();
    GraphOptions polygonGraphOptions = polygonLayer.graphOptions();

    assert pointGraphOptions.edgeType() == EdgeType.DIRECTED;
    assert lineGraphOptions.siblingPairs() != SiblingPairs.CREATE;
    assert lineGraphOptions.siblingPairs() != SiblingPairs.REQUIRE;
    assert polygonGraphOptions.edgeType() == EdgeType.DIRECTED;

    this.options = options;

    normalizedPointLayer.outputGraphOptions = pointGraphOptions;
    normalizedLineLayer.outputGraphOptions = lineGraphOptions;
    normalizedPolygonLayer.outputGraphOptions = polygonGraphOptions;

    normalizedPointLayer.inputGraphOptions = new GraphOptions(pointGraphOptions);
    normalizedLineLayer.inputGraphOptions = new GraphOptions(lineGraphOptions);
    normalizedPolygonLayer.inputGraphOptions = new GraphOptions(polygonGraphOptions);

    /*
     * Set the GraphOptions for the input graphs to ensure that (1) they share a common set of
     * vertices, (2) degenerate edges are kept only if they are isolated, and (3) multiple copies of
     * siblings pairs are discarded. (Note that there may be multiple copies of isolated degenerate
     * edges; clients can eliminate them if desired using DuplicateEdges::MERGE.)
     */
    normalizedPointLayer.inputGraphOptions.setAllowVertexFiltering(false);
    normalizedLineLayer.inputGraphOptions.setAllowVertexFiltering(false);
    normalizedLineLayer.inputGraphOptions.setDegenerateEdges(DegenerateEdges.DISCARD_EXCESS);
    normalizedPolygonLayer.inputGraphOptions.setAllowVertexFiltering(false);
    normalizedPolygonLayer.inputGraphOptions.setDegenerateEdges(DegenerateEdges.DISCARD_EXCESS);
    normalizedPolygonLayer.inputGraphOptions.setSiblingPairs(SiblingPairs.DISCARD_EXCESS);
  }

  public S2ClosedSetNormalizer(
      S2BuilderLayer pointLayer, S2BuilderLayer lineLayer, S2BuilderLayer polygonLayer) {
    this(new Options(), pointLayer, lineLayer, polygonLayer);
  }

  /** Returns the layer that should receive points to be normalized. */
  public S2BuilderLayer pointLayer() {
    return normalizedPointLayer.inputLayer;
  }

  /** Returns the layer that should receive lines to be normalized. */
  public S2BuilderLayer lineLayer() {
    return normalizedLineLayer.inputLayer;
  }

  /** Returns the layer that should receive polygons to be normalized. */
  public S2BuilderLayer polygonLayer() {
    return normalizedPolygonLayer.inputLayer;
  }

  private void maybeRun(S2Error error) {
    if (normalizedPointLayer.inputGraph == null
        || normalizedLineLayer.inputGraph == null
        || normalizedPolygonLayer.inputGraph == null) {
      return;
    }
    run(error);
  }

  private void run(S2Error error) {
    S2BuilderGraph pointGraph = normalizedPointLayer.inputGraph;
    S2BuilderGraph lineGraph = normalizedLineLayer.inputGraph;
    S2BuilderGraph polygonGraph = normalizedPolygonLayer.inputGraph;

    assert pointGraph.options().equals(this.normalizedPointLayer.inputGraphOptions);
    assert lineGraph.options().equals(this.normalizedLineLayer.inputGraphOptions);
    assert polygonGraph.options().equals(this.normalizedPolygonLayer.inputGraphOptions);
    if (options.suppressLowerDimensions()) {
      // Build the auxiliary data needed to suppress lower dimensions.
      polygonInputEdgeIds = polygonGraph.getInEdgeIds();
      excludeFromPointOutput = new boolean[pointGraph.numVertices()];
      Arrays.fill(excludeFromPointOutput, false);
      initExcludeFromPointOutput(lineGraph);
      initExcludeFromPointOutput(polygonGraph);
    }

    // Compute the edges that belong in the output graphs.
    normalizeEdges();

    /*
     * If any edges were added or removed, we need to run S2BuilderGraph.processEdges to ensure that
     * the edges satisfy the requested GraphOptions. Note that since edges are never added to
     * dimension 2, we can use the edge count to test whether any edges were removed. If no edges
     * were removed from dimension 2, then no edges were added to dimension 1, and so we can again
     * use the edge count to test whether any edges were removed, etc.
     */
    boolean polygonsModified = normalizedPolygonLayer.newEdges.size() != polygonGraph.numEdges();
    boolean linesModified =
        polygonsModified || normalizedLineLayer.newEdges.size() != lineGraph.numEdges();
    boolean pointsModified =
        linesModified || normalizedPointLayer.newEdges.size() != pointGraph.numEdges();
    if (!pointsModified && !linesModified && !polygonsModified) {
      // No edges added or removed. Copy the graphs to ensure that they have the GraphOptions that
      // were originally requested.
      normalizedPointLayer.outputGraph =
          copyGraphWithOptions(pointGraph, normalizedPointLayer.outputGraphOptions);
      normalizedLineLayer.outputGraph =
          copyGraphWithOptions(lineGraph, normalizedLineLayer.outputGraphOptions);
      normalizedPolygonLayer.outputGraph =
          copyGraphWithOptions(polygonGraph, normalizedPolygonLayer.outputGraphOptions);
    } else {
      // Make a copy of inputEdgeIdSetLexicon so that ProcessEdges can merge edges if necessary.
      var inputEdgeIdSetLexicon = new IdSetLexicon(pointGraph.inputEdgeIdSetLexicon());
      normalizedPointLayer.newInputEdgeIdSetLexicon = inputEdgeIdSetLexicon;
      normalizedLineLayer.newInputEdgeIdSetLexicon = inputEdgeIdSetLexicon;
      normalizedPolygonLayer.newInputEdgeIdSetLexicon = inputEdgeIdSetLexicon;

      buildOutputGraph(pointGraph, normalizedPointLayer, pointsModified, error);
      buildOutputGraph(lineGraph, normalizedLineLayer, linesModified, error);
      buildOutputGraph(polygonGraph, normalizedPolygonLayer, polygonsModified, error);
    }

    boolean unused =
        normalizedPointLayer.outputLayer.build(normalizedPointLayer.outputGraph, error)
            && normalizedLineLayer.outputLayer.build(normalizedLineLayer.outputGraph, error)
            && normalizedPolygonLayer.outputLayer.build(normalizedPolygonLayer.outputGraph, error);
  }

  /**
   * For each non-degenerate edge in the given graph, set excludeFromPointOutput[v] = true for both
   * the edge's vertices.
   */
  private void initExcludeFromPointOutput(S2BuilderGraph graph) {
    for (int e = 0; e < graph.numEdges(); ++e) {
      int first = graph.edges().getFirst(e);
      int second = graph.edges().getSecond(e);
      if (first != second) {
        excludeFromPointOutput[first] = true;
        excludeFromPointOutput[second] = true;
      }
    }
  }

  private boolean isExcludedFromPointOutput(int vertex) {
    return options.suppressLowerDimensions() && excludeFromPointOutput[vertex];
  }

  @VisibleForTesting
  static S2BuilderGraph copyGraphWithOptions(S2BuilderGraph graph, GraphOptions options) {
    IntVector inputEdgeIdSetIds = new IntVector();
    inputEdgeIdSetIds.addAll(graph.inputEdgeIdSetIds());
    IntVector labelSetIds = new IntVector();
    labelSetIds.addAll(graph.labelSetIds());
    return new S2BuilderGraph(
        new GraphOptions(options),
        new ArrayList<>(graph.vertices()),
        new EdgeList(graph.edges()),
        inputEdgeIdSetIds,
        new IdSetLexicon(graph.inputEdgeIdSetLexicon()),
        labelSetIds,
        new IdSetLexicon(graph.labelSetLexicon()),
        graph.isFullPolygonPredicate());
  }

  private static void buildOutputGraph(
      S2BuilderGraph inputGraph, NormalizedLayer outputLayer, boolean modified, S2Error error) {
    if (modified) {
      S2BuilderGraph.processEdges(
          outputLayer.outputGraphOptions,
          outputLayer.newEdges,
          outputLayer.newInputEdgeIds,
          outputLayer.newInputEdgeIdSetLexicon,
          error);
    }
    IntVector labelSetIds = new IntVector();
    labelSetIds.addAll(inputGraph.labelSetIds());
    outputLayer.outputGraph =
        new S2BuilderGraph(
            outputLayer.outputGraphOptions,
            new ArrayList<>(inputGraph.vertices()),
            outputLayer.newEdges,
            outputLayer.newInputEdgeIds,
            outputLayer.newInputEdgeIdSetLexicon,
            labelSetIds,
            new IdSetLexicon(inputGraph.labelSetLexicon()),
            inputGraph.isFullPolygonPredicate());
  }

  private static int advance(S2BuilderGraph graph, int index, Edge out) {
    if (++index >= graph.numEdges()) {
      out.setSentinel();
      return Integer.MAX_VALUE;
    }
    graph.edges().get(index, out);
    return index;
  }

  private int advanceIncoming(S2BuilderGraph graph, int index, Edge out) {
    if (polygonInputEdgeIds == null || ++index >= polygonInputEdgeIds.size()) {
      out.setSentinel();
      return Integer.MAX_VALUE;
    }
    graph.edges().getReverse(polygonInputEdgeIds.get(index), out);
    return index;
  }

  private void normalizeEdges() {
    S2BuilderGraph pointGraph = normalizedPointLayer.inputGraph;
    S2BuilderGraph lineGraph = normalizedLineLayer.inputGraph;
    S2BuilderGraph polygonGraph = normalizedPolygonLayer.inputGraph;

    // Find the degenerate polygon edges and sibling pairs, and classify each edge as belonging to
    // either a shell or a hole.
    var degeneracies = S2PolygonDegeneracyFinder.findPolygonDegeneracies(polygonGraph);
    int iDegeneracy = 0;

    // Walk through the three edge vectors performing a merge join. We also maintain positions in
    // two other auxiliary vectors: the vector of sorted polygon degeneracies (degeneracies), and
    // the vector of incoming polygon edges (if we are suppressing lower-dimensional duplicate
    // edges).
    Edge ePoint = new Edge();
    int iPoint = advance(pointGraph, -1, ePoint);
    Edge eLine = new Edge();
    int iLine = advance(lineGraph, -1, eLine);
    Edge ePolygon = new Edge();
    int iPolygon = advance(polygonGraph, -1, ePolygon);

    Edge eIncoming = new Edge();
    int iIncoming = advanceIncoming(polygonGraph, -1, eIncoming);

    for (; ; ) {
      if (!eLine.isLessThan(ePolygon) && !ePoint.isLessThan(ePolygon)) {
        if (ePolygon.isSentinel()) {
          break;
        }
        if (iDegeneracy >= degeneracies.size() || degeneracies.edgeId(iDegeneracy) != iPolygon) {
          // Normal polygon edge (not part of a degeneracy).
          addEdge(polygonGraph, iPolygon, normalizedPolygonLayer);
          while (options.suppressLowerDimensions() && eLine.isEqualTo(ePolygon)) {
            iLine = advance(lineGraph, iLine, eLine);
          }
        } else if (!degeneracies.isHole(iDegeneracy++)) {
          // Edge belongs to a degenerate shell.
          if (ePolygon.srcId != ePolygon.dstId) {
            addEdge(polygonGraph, iPolygon, normalizedLineLayer);
            // Since this edge was demoted, make sure that it does not suppress any coincident
            // polyline edge(s).
            while (eLine.isEqualTo(ePolygon)) {
              addEdge(lineGraph, iLine, normalizedLineLayer);
              iLine = advance(lineGraph, iLine, eLine);
            }
          } else {
            // The test below is necessary because a single-vertex polygon shell can be discarded by
            // a polyline edge incident to that vertex.
            if (!isExcludedFromPointOutput(ePolygon.srcId)) {
              addEdge(polygonGraph, iPolygon, normalizedPointLayer);
            }
          }
        }
        iPolygon = advance(polygonGraph, iPolygon, ePolygon);
      } else if (!ePoint.isLessThan(eLine)) {
        if (eLine.srcId != eLine.dstId) {
          // Non-degenerate polyline edge. (Note that polygonInputEdgeIds is empty whenever
          // "suppressLowerDimensions" is false.)
          while (eIncoming.isLessThan(eLine)) {
            iIncoming = advanceIncoming(polygonGraph, iIncoming, eIncoming);
          }
          if (!eLine.isEqualTo(eIncoming)) {
            addEdge(lineGraph, iLine, normalizedLineLayer);
          }
        } else {
          // Degenerate polyline edge.
          if (!isExcludedFromPointOutput(eLine.srcId)) {
            addEdge(lineGraph, iLine, normalizedPointLayer);
          }
          if (lineGraph.options().edgeType() == EdgeType.UNDIRECTED) {
            ++iLine;
          }
        }
        iLine = advance(lineGraph, iLine, eLine);
      } else {
        // Input point.
        if (!isExcludedFromPointOutput(ePoint.srcId)) {
          addEdge(pointGraph, iPoint, normalizedPointLayer);
        }
        iPoint = advance(pointGraph, iPoint, ePoint);
      }
    }
  }

  private static void addEdge(S2BuilderGraph in, int edgeId, NormalizedLayer out) {
    out.newEdges.add(in.edges().getSrcId(edgeId), in.edges().getDstId(edgeId));
    out.newInputEdgeIds.add(in.inputEdgeIdSetId(edgeId));
  }

  /**
   * A simple layer implementation that gathers graphs from each dimension's input layer and
   * processes them together when graphs for all dimensions have been gathered.
   */
  private class DimensionLayer implements S2BuilderLayer {
    private final NormalizedLayer layerData;

    public DimensionLayer(NormalizedLayer layerData) {
      this.layerData = layerData;
    }

    @Override
    public GraphOptions graphOptions() {
      return layerData.inputGraphOptions;
    }

    @Override
    public boolean build(S2BuilderGraph graph, S2Error error) {
      layerData.inputGraph = graph;
      maybeRun(error);
      return error.ok();
    }
  }
}
