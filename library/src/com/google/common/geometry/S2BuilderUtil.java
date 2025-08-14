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

import static com.google.common.geometry.S2.DBL_EPSILON;

import com.google.common.base.Preconditions;
import com.google.common.collect.Multiset;
import com.google.common.collect.Multisets;
import com.google.common.collect.TreeMultiset;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.IsFullPolygonPredicate;
import com.google.common.geometry.S2BuilderGraph.EdgeList;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Comparator;
import java.util.List;

/**
 * A collection of utility classes and methods used by S2Builder, S2BooleanOperation, and their unit
 * tests.
 */
public class S2BuilderUtil {
  private S2BuilderUtil() {}

  /**
   * Given a List of some Comparable 'T', removes adjacent duplicates and truncates the list.
   * Adjacent objects are determined to be duplicates (i.e. equal) using compareTo().
   */
  static <T extends Comparable<T>> void deduplicateSortedList(List<T> input) {
    if (input.size() < 2) {
      return;
    }

    // Copy downward from src to dst, skipping src locations that compare equal to dst.
    int dst = 0;
    for (int src = 1; src < input.size(); src++) {
      if (input.get(src).compareTo(input.get(dst)) != 0) {
        input.set(++dst, input.get(src));
      }
    }
    // Truncate the list to dst + 1 elements by deleting elements from position dst+1 to the end.
    input.subList(dst + 1, input.size()).clear();
  }

  /**
   * GraphShape is an S2Shape implementation used to represent the entire collection of S2Builder
   * input edges, or all the edges in an S2BuilderGraph. Edges are represented as pairs of indices
   * into a vertex list to save space.
   *
   * <p>The GraphShape has a single-edge chain for every edge.
   */
  static class GraphShape implements S2Shape {
    /** A reference to the provided list of vertices. */
    private final List<S2Point> vertices;

    /** A reference to the provided EdgeList, which indexes into the provided vertices. */
    private final EdgeList edges;

    /** Constructs a new GraphShape wrapping the vertices and edges of the given S2BuilderGraph. */
    public GraphShape(S2BuilderGraph graph) {
      this.edges = graph.edges();
      this.vertices = graph.vertices();
    }

    /**
     * Constructs a new GraphShape holding references to the provided EdgeList and List of vertices,
     * which should not be modified for the life of this object.
     */
    public GraphShape(EdgeList edges, List<S2Point> vertices) {
      this.edges = edges;
      this.vertices = vertices;
    }

    @Override
    public int dimension() {
      return 1;
    }

    @Override
    public int numEdges() {
      return edges.size();
    }

    @Override
    public void getEdge(int edgeId, MutableEdge result) {
      result.set(vertices.get(edges.getSrcId(edgeId)), vertices.get(edges.getDstId(edgeId)));
    }

    @Override
    public boolean hasInterior() {
      return false;
    }

    @Override
    public boolean containsOrigin() {
      return false;
    }

    @Override
    public int numChains() {
      // Every edge is in its own chain.
      return edges.size();
    }

    @Override
    public int getChainStart(int chainId) {
      Preconditions.checkElementIndex(chainId, numChains());
      return chainId;
    }

    @Override
    public int getChainLength(int chainId) {
      Preconditions.checkElementIndex(chainId, numChains());
      return 1;
    }

    @Override
    public void getChainEdge(int chainId, int offset, MutableEdge result) {
      // getChainLength validates chainId.
      Preconditions.checkElementIndex(offset, getChainLength(chainId));
      getEdge(chainId, result);
    }

    @Override
    public S2Point getChainVertex(int chainId, int edgeOffset) {
      Preconditions.checkElementIndex(chainId, numChains());
      Preconditions.checkElementIndex(edgeOffset, 2); // Each chain is one edge with two vertices.
      return edgeOffset == 0
          ? vertices.get(edges.getSrcId(chainId))
          : vertices.get(edges.getDstId(chainId));
    }

    @Override
    public void getChainPosition(int edgeId, ChainPosition result) {
      // Each edge is in its own chain.
      Preconditions.checkElementIndex(edgeId, numEdges());
      result.set(edgeId, 0);
    }
  }

  /**
   * A GraphClone copies an S2BuilderGraph and owns the underlying data, unlike S2BuilderGraph,
   * which is just a view of data.
   */
  static class GraphClone {
    private GraphOptions options;
    private ArrayList<S2Point> vertices;
    private EdgeList edges;
    private IntArrayList inputEdgeIdSetIds;
    private IdSetLexicon inputEdgeIdSetLexicon;
    private IntVector labelSetIds;
    private IdSetLexicon labelSetLexicon;
    private IsFullPolygonPredicate isFullPolygonPredicate;

    private S2BuilderGraph graph;

    /** Construct an empty GraphClone. */
    public GraphClone() {}

    /**
     * Construct a GraphClone as a deep copy of the given S2BuilderGraph and its underlying data.
     */
    public GraphClone(S2BuilderGraph graph) {
      init(graph);
    }

    /**
     * Discard existing data and make this GraphClone be a deep copy of the given S2BuilderGraph and
     * its underlying data.
     */
    public void init(S2BuilderGraph graph) {
      // Copy the data
      this.options = new GraphOptions(graph.options());
      this.vertices = new ArrayList<>(graph.vertices());
      this.edges = new EdgeList(graph.edges());
      this.inputEdgeIdSetIds = new IntArrayList(graph.inputEdgeIdSetIds());
      this.inputEdgeIdSetLexicon = new IdSetLexicon(graph.inputEdgeIdSetLexicon());
      this.labelSetIds = IntVector.copyOf(graph.labelSetIds());
      this.labelSetLexicon = new IdSetLexicon(graph.labelSetLexicon());
      this.isFullPolygonPredicate = graph.isFullPolygonPredicate();

      // Wrap the data in a S2BuilderGraph.
      this.graph =
          new S2BuilderGraph(
              options,
              vertices,
              edges,
              inputEdgeIdSetIds,
              inputEdgeIdSetLexicon,
              labelSetIds,
              labelSetLexicon,
              isFullPolygonPredicate);
    }

    /** Returns the S2BuilderGraph of this GraphClone. */
    public S2BuilderGraph graph() {
      return graph;
    }
  }

  /**
   * An S2BuilderLayer implementation that wraps an S2BuilderShapesLayer, and when build() completes
   * successfully, adds the built shapes (if not empty) to an index.
   */
  public static class IndexedLayer<L extends S2BuilderShapesLayer> implements S2BuilderShapesLayer {
    private final S2ShapeIndex index;
    private final L layer;

    public IndexedLayer(S2ShapeIndex index, L wrappedLayer) {
      this.index = index;
      this.layer = wrappedLayer;
    }

    @Override
    public GraphOptions graphOptions() {
      return layer.graphOptions();
    }

    @Override
    public boolean build(S2BuilderGraph graph, S2Error error) {
      boolean success = layer.build(graph, error);
      if (success) {
        for (S2Shape shape : layer.shapes()) {
          if (!shape.isEmpty()) {
            index.add(shape);
          }
        }
      }
      return success;
    }

    @Override
    public Iterable<? extends S2Shape> shapes() {
      return layer.shapes();
    }

    @Override
    public String toString() {
      return "IndexedLayer wrapping a " + layer;
    }
  }

  /**
   * An S2BuilderLayer implementation that copies an S2BuilderGraph into a GraphClone object, which
   * owns the underlying data, unlike S2BuilderGraph itself.
   */
  public static class GraphCloningLayer implements S2BuilderLayer {
    private final GraphOptions graphOptions;
    private final GraphClone graphClone;

    /** Constructs a GraphCloningLayer. */
    public GraphCloningLayer(GraphOptions graphOptions, GraphClone graphClone) {
      this.graphOptions = graphOptions;
      this.graphClone = graphClone;
    }

    @Override
    public GraphOptions graphOptions() {
      return graphOptions;
    }

    @Override
    public boolean build(S2BuilderGraph graph, S2Error error) {
      graphClone.init(graph);
      return true;
    }

    @Override
    public String toString() {
      return "GraphCloningLayer with GraphOptions " + graphOptions;
    }
  }

  /**
   * An S2BuilderLayer implementation that expects that the edges in the S2BuilderGraph passed to
   * its build() method should match the edges in the given S2ShapeIndex (including multiplicities).
   * This allows testing whether an algorithm produces a given multiset of edges without needing to
   * specify a particular ordering of those edges.
   *
   * <p>Edges A and B are considered matching if their corresponding endpoints have a 3D Cartesian
   * distance less than a tolerance which may be specified, or otherwise defaults to a very small
   * but non-zero value. It may be set to zero to force exact matching.
   */
  static class IndexMatchingLayer implements S2BuilderLayer {
    private final S2Builder.GraphOptions graphOptions;
    private final S2ShapeIndex expectedIndex;
    private final int dimension;
    private final double tolerance;
    private final Multiset<S2Edge> expectedEdges;
    private final Multiset<S2Edge> actualEdges;

    /**
     * Constructs an IndexMatchingLayer that will check whether the edges passed to its build()
     * method match the edges in the provided S2ShapeIndex (including multiplicities). If any
     * differences are found, sets "error" to a descriptive error message.
     *
     * <p>If "dimension" is non-negative then only shapes of the given dimension are used. (This
     * allows for use with classes such as S2BooleanOperation that output one S2BuilderGraph for
     * each dimension.)
     *
     * <p>Edges are considered matching if their corresponding endpoints are separated by a 3D
     * Cartesian distance less than or equal to the given 'tolerance'.
     */
    public IndexMatchingLayer(
        S2Builder.GraphOptions graphOptions, S2ShapeIndex index, int dimension, double tolerance) {
      this.graphOptions = new S2Builder.GraphOptions(graphOptions);
      this.expectedIndex = index;
      this.dimension = dimension;
      this.tolerance = tolerance;
      this.expectedEdges = TreeMultiset.create(tolerantEdgeComparator(tolerance));
      this.actualEdges = TreeMultiset.create(tolerantEdgeComparator(tolerance));

      computeExpectedEdges();
    }

    /** As above but with the default small tolerance for matching edges. */
    public IndexMatchingLayer(
        S2Builder.GraphOptions graphOptions, S2ShapeIndex index, int dimension) {
      this(graphOptions, index, dimension, DBL_EPSILON);
    }

    /**
     * As above but with a default dimension of -1 to use shapes of all dimensions, and the default
     * small tolerance for matching edges.
     */
    public IndexMatchingLayer(S2Builder.GraphOptions graphOptions, S2ShapeIndex index) {
      this(graphOptions, index, -1);
    }

    /**
     * Orders S2Edges lexicographically, but considers them to be equal if their endpoints are at
     * most 'tolerance' apart. Note that this sense of equality is reflexive but NOT transitive.
     */
    public Comparator<S2Edge> tolerantEdgeComparator(double tolerance) {
      return (S2Edge left, S2Edge right) -> {
        if (left.getStart().getDistance(right.getStart()) > tolerance) {
          // Edges are not equal, order lexicographically by start point.
          return left.getStart().compareTo(right.getStart());
        }
        // Start points match. Consider end points.
        if (left.getEnd().getDistance(right.getEnd()) > tolerance) {
          // Edges are not equal, order lexicographically by end point.
          return left.getEnd().compareTo(right.getEnd());
        }
        // Start and end points match.
        return 0;
      };
    }

    /**
     * Computes expected edges. This is done as soon as the layer is constructed so that they can be
     * used in toString(), which is helpful when debugging.
     */
    private void computeExpectedEdges() {
      for (S2Shape expectedShape : expectedIndex.getShapes()) {
        if (dimension >= 0 && expectedShape.dimension() != dimension) {
          continue;
        }
        S2Shape.MutableEdge mutableEdge = new S2Shape.MutableEdge();
        for (int e = expectedShape.numEdges(); --e >= 0; ) {
          expectedShape.getEdge(e, mutableEdge);
          expectedEdges.add(new S2Edge(mutableEdge.getStart(), mutableEdge.getEnd()));
        }
      }
    }

    private String toString(Collection<S2Edge> edges) {
      StringBuilder msg = new StringBuilder();
      for (S2Edge edge : edges) {
        if (msg.length() > 0) {
          msg.append("; ");
        }
        msg.append("[")
            .append(S2TextFormat.toString(edge.getStart()))
            .append(", ")
            .append(S2TextFormat.toString(edge.getEnd()))
            .append("]");
        // The following may be useful in debugging. In some cases, lat/lng to 15 decimal places
        // does not have sufficient precision to show differences. To print the points as XYZ:
        //    .append(" {")
        //    .append(edge.getStart().toString())
        //    .append(", ")
        //    .append(edge.getEnd().toString())
        //    .append("}");
      }
      return msg.toString();
    }

    @Override
    public String toString() {
      return "IndexMatchingLayer for dimension "
          + dimension
          + " with tolerance="
          + tolerance + ", and "
          + expectedEdges.size()
          + " expected edges: "
          + toString(expectedEdges);
    }

    @Override
    public GraphOptions graphOptions() {
      return graphOptions;
    }

    @Override
    public boolean build(S2BuilderGraph g, S2Error error) {
      actualEdges.clear();
      for (int e = 0; e < g.numEdges(); ++e) {
        actualEdges.add(
            new S2Edge(g.vertex(g.edges().getSrcId(e)), g.vertex(g.edges().getDstId(e))));
      }

      Multiset<S2Edge> missingEdges = Multisets.difference(expectedEdges, actualEdges);
      Multiset<S2Edge> extraEdges = Multisets.difference(actualEdges, expectedEdges);

      if (!missingEdges.isEmpty() || !extraEdges.isEmpty()) {
        // There may be errors in more than one dimension, so we append to the existing error text.
        String label = "";
        if (dimension >= 0) {
          label = Platform.formatString("Dimension %d: ", dimension);
        }
        error.init(
            S2Error.Code.FAILED_PRECONDITION,
            "%s%s\n   Missing edges: %s\n   Extra edges: %s\n",
            error.text(), // previous layer's error message, if any, ends in newline
            label,
            toString(missingEdges),
            toString(extraEdges));
      }
      return error.ok();
    }
  }
}
