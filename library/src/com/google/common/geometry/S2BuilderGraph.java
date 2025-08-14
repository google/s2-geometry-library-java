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

import static java.lang.Math.max;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2Builder.IsFullPolygonPredicate;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IdSetLexicon.IdSet;
import com.google.common.geometry.primitives.IntPairVector;
import com.google.common.geometry.primitives.IntVector;
import com.google.common.geometry.primitives.Ints.IntBiConsumer;
import com.google.common.geometry.primitives.Ints.IntConsumer;
import com.google.common.geometry.primitives.Ints.IntList;
import com.google.common.geometry.primitives.Ints.IntSequence;
import com.google.common.geometry.primitives.Ints.OfInt;
import com.google.common.geometry.primitives.Sorter;
import com.google.common.primitives.Ints;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntComparator;
import it.unimi.dsi.fastutil.ints.IntIterator;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

/**
 * A collection of snapped edges that is passed to a {@link S2BuilderLayer} for assembly. (Example
 * Layers include polygons, polylines, and polygon meshes.) You will only need this interface if you
 * want to implement a new Layer subtype.
 *
 * <p>The graph consists of vertices and directed edges. Vertices are S2Points, provided to the
 * S2BuilderGraph in a list indexed sequentially starting from zero. The index of a vertex in this
 * list is the vertex id.
 *
 * <p>An Edge is represented as a pair of vertex ids. The edges are stored in an EdgeList and sorted
 * in lexicographic order, so that the outgoing edges from a particular vertex form a contiguous
 * range. The index of an edge in this list is the edge id.
 */
@SuppressWarnings("Assertion")
public class S2BuilderGraph {
  // An edge from SENTINEL_ID to SENTINEL_ID is non-existent and compares as larger than any real
  // edge.
  private static final int SENTINEL_ID = Integer.MAX_VALUE;

  // TODO(user): MAX_INPUT_EDGE_ID is kMaxInputEdgeId in C++. The name implies this is a
  // valid edge id, but the comment says it is not. Also, NO_INPUT_EDGE_ID (kNoInputEdgeId) is
  // less than max, so numerically appears to be a valid input edge id. But it is not.
  // (It isn't actually used anywhere?)

  /** Defines a value larger than any valid InputEdgeId. */
  public static final int MAX_INPUT_EDGE_ID = Integer.MAX_VALUE;

  /**
   * The following value of InputEdgeId means that an edge does not corresponds to any input edge.
   */
  public static final int NO_INPUT_EDGE_ID = Integer.MAX_VALUE - 1;

  // TODO(torrey): Move the definition of GraphOptions into S2BuilderGraph.java.
  /** The options for this S2BuilderGraph. */
  private final S2Builder.GraphOptions options;

  /** The vertices in this S2BuilderGraph. The index of a Vertex in this list is the VertexId. */
  private final List<S2Point> vertices;

  /** The edges in this S2BuilderGraph. The index of an Edge in this list is the EdgeId. */
  private final EdgeList edges;

  /**
   * A vector mapping edge id to IdSet ids. The keys are edge ids, and values are InputEdgeIdSetIds,
   * i.e. the ids of IdSets in the inputEdgeIdSetLexicon. Each of those IdSets is the set of input
   * edge ids that were mapped to the edge.
   */
  private final IntArrayList inputEdgeIdSetIds;

  /** Mapping from ints which are InputEdgeIdSetIds to sets of input edge ids. */
  private final IdSetLexicon inputEdgeIdSetLexicon;

  /**
   * A vector mapping input edge ids to LabelSet ids. The keys are input edge ids, and values are
   * LabelSetIds, i.e. the ids of IdSets in the labelSetLexicon. Each of those IdSets is the set of
   * labels were attached to the given input edge. This vector may be empty to indicate that no
   * labels are present.
   */
  private final IntList labelSetIds;

  /**
   * Mapping from ints which are LabelSetIds to a sets of Labels, which are also ints. This lexicon
   * will exist even if no labels are present.
   */
  private final IdSetLexicon labelSetLexicon;

  /**
   * If the geometry in this graph has no edges, this predicate is used to determine if it is a full
   * or empty polygon. See {@link S2Builder.IsFullPolygonPredicate}.
   */
  private final IsFullPolygonPredicate isFullPolygonPredicate;

  /**
   * The constructor for S2BuilderGraph.
   *
   * @param options the GraphOptions that will be used to build the S2BuilderGraph. In some cases
   *     these can be different than the options provided by the Layer.
   * @param vertices a list of S2Points. The index of a point in the list is its VertexId.
   * @param edges a list of Edges, which are VertexId pairs. They must be provided presorted in
   *     lexicographic order. The index of an edge in this list is its EdgeId.
   * @param inputEdgeIdSetIds a vector indexed by int EdgeIds that allows access to the set of
   *     InputEdgeIds (ints) that were mapped to the given edge id, by looking up the returned value
   *     (an int that is a InputEdgeIdSetId) in "inputEdgeIdSetLexicon".
   * @param inputEdgeIdSetLexicon a class that maps an int InputEdgeIdSetId to a set of ints that
   *     are InputEdgeIds.
   * @param labelSetIds a vector indexed by InputEdgeId that allows access to the set of labels that
   *     were attached to the given input edge, by looking up the returned value (an int LabelSetId)
   *     in the "labelSetLexicon". This vector may be empty to indicate that no labels are present.
   * @param labelSetLexicon a class that maps an int LabelSetId to a set of S2Builder labels (ints).
   *     Must be provided even if no labels are present.
   * @param isFullPolygonPredicate a predicate called to determine whether a graph consisting only
   *     of polygon degeneracies represents the empty polygon or the full polygon (see {@link
   *     S2Builder.IsFullPolygonPredicate} for details).
   */
  public S2BuilderGraph(
      GraphOptions options,
      List<S2Point> vertices,
      EdgeList edges,
      IntArrayList inputEdgeIdSetIds,
      IdSetLexicon inputEdgeIdSetLexicon,
      IntList labelSetIds,
      IdSetLexicon labelSetLexicon,
      IsFullPolygonPredicate isFullPolygonPredicate) {
    this.options = options;
    this.vertices = vertices;
    this.edges = edges;
    this.inputEdgeIdSetIds = inputEdgeIdSetIds;
    this.inputEdgeIdSetLexicon = inputEdgeIdSetLexicon;
    this.labelSetIds = labelSetIds;
    this.labelSetLexicon = labelSetLexicon;
    this.isFullPolygonPredicate = isFullPolygonPredicate;
    // TODO(torrey): Add a debug-mode check that the given edges are sorted in lexicographic order.
  }

  /** Returns a copy of the GraphOptions used by this S2BuilderGraph. */
  public GraphOptions options() {
    return new GraphOptions(options);
  }

  // TODO(torrey): Consider removing these redundant getters; clients can use vertices().size() etc.
  /** Returns the number of vertices in the graph. */
  public int numVertices() {
    return vertices.size();
  }

  /** Returns the vertex with the given VertexId (its index in vertices). */
  public S2Point vertex(int vertexId) {
    return vertices.get(vertexId);
  }

  /** Returns the entire list of vertices. */
  public List<S2Point> vertices() {
    return vertices;
  }

  // TODO(torrey): Consider removing these redundant getters, clients can use edges().
  /** Returns the total number of Edges in the graph. */
  public int numEdges() {
    return edges.size();
  }

  /** Returns the source vertex id of the Edge with the given edgeId (its index in edges). */
  public int edgeSrcId(int edgeId) {
    return edges.getSrcId(edgeId);
  }

  /** Returns the destination vertex id of the Edge with the given edgeId (its index in edges). */
  public int edgeDstId(int edgeId) {
    return edges.getDstId(edgeId);
  }

  /** Returns the source vertex of the Edge with the given edgeId (its index in edges). */
  public S2Point edgeSrc(int edgeId) {
    return vertices.get(edges.getSrcId(edgeId));
  }

  /** Returns the destination vertex of the Edge with the given edgeId (its index in edges). */
  public S2Point edgeDst(int edgeId) {
    return vertices.get(edges.getDstId(edgeId));
  }

  /** Returns the entire ordered list of Edges. */
  public EdgeList edges() {
    return edges;
  }

  /**
   * Returns an IntComparator that compares EdgeIds in the given EdgeList by destination vertex id,
   * breaking ties with source vertex ids, and breaking ties for equal edges with edge ids.
   */
  @VisibleForTesting
  static IntComparator byDstVertexEdge(EdgeList edges) {
    return (int edgeIdA, int edgeIdB) -> {
      int dstVertexIdA = edges.getDstId(edgeIdA);
      int dstVertexIdB = edges.getDstId(edgeIdB);
      if (dstVertexIdA != dstVertexIdB) {
        return Integer.compare(dstVertexIdA, dstVertexIdB);
      }
      int srcVertexIdA = edges.getSrcId(edgeIdA);
      int srcVertexIdB = edges.getSrcId(edgeIdB);
      if (srcVertexIdA != srcVertexIdB) {
        return Integer.compare(srcVertexIdA, srcVertexIdB);
      }
      return Integer.compare(edgeIdA, edgeIdB);
    };
  }

  private static void fillConsecutive(IntArrayList list, int n) {
    list.clear();
    for (int i = 0; i < n; i++) {
      list.add(i);
    }
  }

  /**
   * Returns EdgeIds sorted in lexicographic order by (destination, origin) vertex ids. All of the
   * incoming edges to each vertex form a contiguous subrange of this ordering.
   */
  IntArrayList getInEdgeIds() {
    IntArrayList inEdgeIds = new IntArrayList(numEdges());
    fillConsecutive(inEdgeIds, numEdges());
    inEdgeIds.sort(byDstVertexEdge(edges));
    return inEdgeIds;
  }

  /**
   * For a graph such that every directed edge has a sibling, returns a map from each edge id to the
   * sibling edge id. This method is identical to getInEdgeIds() except that (1) it requires edges
   * to have siblings, and (2) undirected degenerate edges are grouped together in pairs such that
   * one edge is the sibling of the other. (The sibling of a directed degenerate edge is itself.)
   * Handles duplicate edges correctly and is also consistent with fillLeftTurnMap().
   *
   * <p>REQUIRES: An option is chosen that guarantees sibling pairs: (options.siblingPairs() == {
   * REQUIRE, CREATE } || options.edgeType() == UNDIRECTED)
   */
  IntArrayList getSiblingMap() {
    IntArrayList inEdgeIds = getInEdgeIds();
    makeSiblingMap(inEdgeIds);
    return inEdgeIds;
  }

  /** Verifies that the sibling of a sibling is an identity function. */
  private static boolean validateSiblingMap(IntArrayList siblingMap, int numEdges) {
    for (int edgeId = 0; edgeId < numEdges; ++edgeId) {
      if (siblingMap.getInt(siblingMap.getInt(edgeId)) != edgeId) {
        return false;
      }
    }
    return true;
  }

  /** Verifies that edges(edgeId) == reverse(edges(inEdgeIds(edgeId))). */
  private boolean validateInEdgeIds(IntArrayList inEdgeIds) {
    for (int edgeId = 0; edgeId < numEdges(); ++edgeId) {
      if (edges.getSrcId(edgeId) != edges.getDstId(inEdgeIds.getInt(edgeId))
          || edges.getDstId(edgeId) != edges.getSrcId(inEdgeIds.getInt(edgeId))) {
        return false;
      }
    }
    return true;
  }

  /**
   * Like getSiblingMap(), but constructs the map from edge id to sibling edge id starting from the
   * provided list of incoming edge ids, which is returned by getInEdgeIds().
   *
   * <p>(This operation is a no-op unless undirected degenerate edges are present, in which case
   * such edges are grouped together in pairs to satisfy the requirement that every edge must have a
   * sibling edge.)
   */
  void makeSiblingMap(IntArrayList inEdgeIds) {
    Preconditions.checkState(
        options.siblingPairs() == SiblingPairs.REQUIRE
            || options.siblingPairs() == SiblingPairs.CREATE
            || options.edgeType() == EdgeType.UNDIRECTED);
    assert validateInEdgeIds(inEdgeIds);

    if (options.edgeType() == EdgeType.DIRECTED
        || options.degenerateEdges() == DegenerateEdges.DISCARD) {
      assert validateSiblingMap(inEdgeIds, numEdges());
      return;
    }

    for (int edgeId = 0; edgeId < numEdges(); ++edgeId) {
      int vertexId = edges.getSrcId(edgeId);
      if (edges.getDstId(edgeId) == vertexId) {
        // Edge is degenerate
        assert edgeId + 1 < numEdges();
        assert edges.getSrcId(edgeId + 1) == vertexId;
        assert edges.getDstId(edgeId + 1) == vertexId;
        assert inEdgeIds.getInt(edgeId) == edgeId;
        assert inEdgeIds.getInt(edgeId + 1) == edgeId + 1;
        inEdgeIds.set(edgeId, edgeId + 1);
        inEdgeIds.set(edgeId + 1, edgeId);
        edgeId++;
      }
    }
    assert validateSiblingMap(inEdgeIds, numEdges());
  }

  /**
   * Returns the set of input edge ids that were snapped to the given edge. ("Input edge ids" are
   * assigned to input edges sequentially in the order they are added to the builder.) For example,
   * if input edges 2 and 17 were snapped to edge 12, then inputEdgeIds(12) returns a set containing
   * the numbers 2 and 17. Example usage:
   *
   * <pre>
   * graph.inputEdgeIds(e).forEach(inputEdgeId -> {
   *   ...
   * });
   * </pre>
   *
   * <p>Please note the following:
   *
   * <ul>
   *   <li>When edge chains are simplified, the simplified edge is assigned all the input edge ids
   *       associated with edges of the chain.
   *   <li>Edges can also have multiple input edge ids due to edge merging (if DuplicateEdges.MERGE
   *       is specified).
   *   <li>Sibling edges automatically created by EdgeType.UNDIRECTED or SiblingPairs.CREATE have an
   *       empty set of input edge ids. (However you can use a LabelFetcher to retrieve the set of
   *       labels associated with both edges of a given sibling pair.)
   * </ul>
   */
  IdSetLexicon.IdSet inputEdgeIds(int edgeId) {
    return inputEdgeIdSetLexicon().idSet(inputEdgeIdSetIds.getInt(edgeId));
  }

  /**
   * Low-level method that returns an InputEdgeIdSetId, which is an integer that is the id of an
   * IdSet containing the set of input edge ids that were snapped to the output edge with the given
   * 'edgeId'. The elements of the IdSet can be accessed using {@link inputEdgeIdSetLexicon()}.
   */
  int inputEdgeIdSetId(int edgeId) {
    return inputEdgeIdSetIds.getInt(edgeId);
  }

  /**
   * Low-level method that returns an IntArrayList of InputEdgeIdSetIds. For each index i in the
   * vector, the element at 'i' is an InputEdgeIdSetId, an integer that is the id of the IdSet
   * containing the input edge ids that were snapped to output edge 'i'.
   */
  IntArrayList inputEdgeIdSetIds() {
    return inputEdgeIdSetIds;
  }

  /** Returns a mapping from an InputEdgeIdSetId to a set of input edge ids. */
  IdSetLexicon inputEdgeIdSetLexicon() {
    return inputEdgeIdSetLexicon;
  }

  /**
   * Returns the minimum input edge id that was snapped to the edge with given id 'edgeId', or -1 if
   * no input edges were snapped. (See {@link SiblingPairs#CREATE}). This is useful for layers that
   * wish to preserve the input edge ordering as much as possible (e.g., to ensure idempotency).
   */
  public int minInputEdgeId(int edgeId) {
    IdSetLexicon.IdSet idSet = inputEdgeIds(edgeId);
    return idSet.isEmpty() ? NO_INPUT_EDGE_ID : idSet.first();
  }

  /**
   * Returns a vector containing the minimum input edge id for every edge, by edge id. If an edge
   * has no input ids, NO_INPUT_EDGE_ID is used.
   */
  IntArrayList getMinInputEdgeIds() {
    IntArrayList minInputEdgeIds = new IntArrayList(numEdges());
    fillMinInputEdgeIds(minInputEdgeIds);
    return minInputEdgeIds;
  }

  /**
   * Fills the given IntArrayList with the minimum input edge id for every edge, by edge id. If an
   * edge has no input ids, such as automatically created edges for EdgeType.UNDIRECTED,
   * NO_INPUT_EDGE_ID is used.
   */
  void fillMinInputEdgeIds(IntArrayList minInputEdgeIds) {
    minInputEdgeIds.clear();
    minInputEdgeIds.size(numEdges());
    for (int edgeId = 0; edgeId < numEdges(); ++edgeId) {
      minInputEdgeIds.set(edgeId, minInputEdgeId(edgeId));
    }
  }

  /**
   * Returns edge ids sorted by minimum input edge id. This is an approximation of the input edge
   * ordering.
   */
  IntArrayList getInputEdgeOrder(IntArrayList minInputEdgeIds) {
    IntArrayList orderedEdgeIds = new IntArrayList(numEdges());
    fillConsecutive(orderedEdgeIds, numEdges());
    orderedEdgeIds.sort(byMinInputEdgeId(minInputEdgeIds));
    return orderedEdgeIds;
  }

  /**
   * Returns the set of labels associated with a given input edge. Example:
   *
   * {@snippet :
   * for (Label label : g.labels(inputEdgeId)) {
   *   ...
   * }
   * }
   *
   * See also LabelFetcher, which returns the labels for a given graph edge.
   */
  public IdSet labels(int inputEdgeId) {
    return labelSetLexicon().idSet(labelSetId(inputEdgeId));
  }

  /**
   * Low-level method that returns an integer representing the set of labels associated with a given
   * input edge. The elements of the IdSet can be accessed using labelSetLexicon().
   */
  public int labelSetId(int inputEdgeId) {
    return labelSetIds().isEmpty() ? IdSetLexicon.EMPTY_SET_ID : labelSetIds().get(inputEdgeId);
  }

  /**
   * Low-level method that returns a vector where each element represents the set of labels
   * associated with a particular input edge. Note that this vector may be empty, which indicates
   * that no labels are present.
   */
  public IntList labelSetIds() {
    return labelSetIds;
  }

  /** Returns a mapping from an int LabelSetId to a set of labels. */
  public IdSetLexicon labelSetLexicon() {
    return labelSetLexicon;
  }

  /**
   * Returns a method that determines whether a graph that consists only of polygon degeneracies
   * represents the empty polygon or the full polygon. See {@link S2Builder.IsFullPolygonPredicate}
   * for details.
   */
  IsFullPolygonPredicate isFullPolygonPredicate() {
    // TODO(user): IsFullPolygonPredicate should be a property of layers.
    return isFullPolygonPredicate;
  }

  /** Returns true if this graph is a full polygon, according to its IsFullPolygonPredicate. */
  boolean isFullPolygon() {
    // TODO(user): IsFullPolygonPredicate should be a property of layers.
    return isFullPolygonPredicate().test(this);
  }

  /**
   * Fills in the provided 'leftTurnMap', which must have length {@code numEdges()}, so it maps each
   * edge e=(v0, v1) id to the following outgoing edge around "v1" in clockwise order. (This
   * corresponds to making a "left turn" at the vertex.) By starting at a given edge and making only
   * left turns, you can construct a loop whose interior does not contain any edges in the same
   * connected component.
   *
   * <p>If the incoming and outgoing edges around a vertex do not alternate perfectly (e.g., there
   * are two incoming edges in a row), then adjacent (incoming, outgoing) pairs are repeatedly
   * matched and removed. This is similar to finding matching parentheses in a string such as
   * "(()())()".
   *
   * <p>For sibling edge pairs, the incoming edge is assumed to immediately follow the outgoing edge
   * in clockwise order. Thus a left turn is made from an edge to its sibling only if there are no
   * other outgoing edges. With respect to the parentheses analogy, a sibling pair is ")(".
   * Similarly, if there are multiple copies of a sibling edge pair then the duplicate incoming and
   * outgoing edges are sorted in alternating order, e.g., ")()(".
   *
   * <p>Degenerate edges (edges from a vertex to itself) are treated as loops consisting of a single
   * edge. This avoids the problem of deciding the connectivity and ordering of such edges when they
   * share a vertex with other edges (possibly including other degenerate edges).
   *
   * <p>If it is not possible to make a left turn from every input edge, this method returns false
   * and sets "error" appropriately. In this situation the left turn map is still valid except that
   * any incoming edge where it is not possible to make a left turn will have its entry set to -1.
   * Otherwise returns true and leaves error unchanged.
   *
   * @param "inEdgeIds" should be equal to getInEdgeIds() or getSiblingMap().
   */
  @CanIgnoreReturnValue
  boolean fillLeftTurnMap(IntArrayList inEdgeIds, int[] leftTurnMap, S2Error error) {
    Preconditions.checkArgument(leftTurnMap.length == numEdges());
    Arrays.fill(leftTurnMap, -1);
    if (numEdges() == 0) {
      return true;
    }

    // The collection of all incoming and outgoing edges to/from the current vertex 'v0'.
    VertexEdgeList v0Edges = new VertexEdgeList(vertices);
    // unmatchedIncomingEdges is a stack of unmatched incoming edge ids to a vertex.
    IntArrayList unmatchedIncomingEdges = new IntArrayList();
    // unmatchedOutgoingEdges is a stack of outgoing edge ids with no previous incoming edge.
    IntArrayList unmatchedOutgoingEdges = new IntArrayList();

    // Walk through the two sorted arrays of edges (outgoing and incoming) and gather all the edges
    // incident to each vertex. Then we sort those edges and add an entry to the left turn map from
    // each incoming edge to the immediately following outgoing edge in clockwise order.
    int out = 0;
    int in = 0;
    final Edge outEdge = new Edge();
    edges.get(out, outEdge);
    final Edge reverseInEdge = new Edge();
    edges.getReverse(inEdgeIds.getInt(in), reverseInEdge);
    final Edge minEdge = new Edge();

    if (EdgeList.compareEdges(outEdge, reverseInEdge) < 0) {
      minEdge.set(outEdge);
    } else {
      minEdge.set(reverseInEdge);
    }

    while (!minEdge.isSentinel()) {
      // Gather all incoming and outgoing edges around the source vertex 'v0' of minEdge into the
      // VertexEdgeList 'v0Edges'.
      int v0 = minEdge.srcId;
      while (minEdge.srcId == v0) {
        int v1 = minEdge.dstId;
        // Count the number of copies of "minEdge" in each direction.
        int outBegin = out;
        int inBegin = in;

        // Advance 'out' and 'outEdge' over copies of minEdge.
        while (outEdge.isEqualTo(minEdge)) {
          // outEdge = (++out == numEdges()) ? sentinel : edge(out);
          if (++out == numEdges()) {
            outEdge.setSentinel();
          } else {
            edges.get(out, outEdge);
          }
        }

        // Advance 'in' and 'reverseInEdge' over copies of minEdge.
        while (reverseInEdge.isEqualTo(minEdge)) {
          if (++in == numEdges()) {
            reverseInEdge.setSentinel();
          } else {
            edges.getReverse(inEdgeIds.getInt(in), reverseInEdge);
          }
        }

        // If minEdge is not degenerate, add those incoming and outgoing edges to v0Edges.
        if (v0 != v1) {
          addVertexEdges(outBegin, out, inBegin, in, v1, v0Edges);
        } else {
          // Otherwise, each degenerate edge simply becomes its own loop.
          for (; inBegin < in; ++inBegin) {
            leftTurnMap[inBegin] = inBegin;
          }
        }

        // Advance to the next minEdge.
        if (EdgeList.compareEdges(outEdge, reverseInEdge) < 0) {
          minEdge.set(outEdge);
        } else {
          minEdge.set(reverseInEdge);
        }
      }

      if (v0Edges.isEmpty()) {
        continue;
      }

      // Sort the edges in clockwise order around "v0".
      v0Edges.sortClockwiseAround(v0);

      // Match incoming with outgoing edges. We do this by keeping a stack of unmatched incoming
      // edges. We also keep a stack of outgoing edges with no previous incoming edge, and match
      // these at the end by wrapping around circularly to the start of the edge ordering.

      // Iterate through the VertexEdgeList "VertexEdges" in clockwise order.
      for (int veIndex = 0; veIndex < v0Edges.size(); ++veIndex) {
        if (v0Edges.incoming(veIndex)) {
          // Add an unmatched incoming edge id.
          unmatchedIncomingEdges.add(inEdgeIds.getInt(v0Edges.indexEdgeId(veIndex)));
        } else if (!unmatchedIncomingEdges.isEmpty()) {
          // e is outgoing, match it with a previous incoming edge as a left turn from the incoming
          // edge to the outgoing edge 'e'.
          leftTurnMap[unmatchedIncomingEdges.popInt()] = v0Edges.indexEdgeId(veIndex);
        } else {
          // Add an unmatched outgoing edge id.
          unmatchedOutgoingEdges.add(v0Edges.indexEdgeId(veIndex));
        }
      }

      // Pair up additional edges using the fact that the ordering is circular.
      reverse(unmatchedOutgoingEdges);
      while (!unmatchedOutgoingEdges.isEmpty() && !unmatchedIncomingEdges.isEmpty()) {
        leftTurnMap[unmatchedIncomingEdges.popInt()] = unmatchedOutgoingEdges.popInt();
      }

      // We only need to process unmatched incoming edges, since we are only responsible for
      // creating left turn map entries for those edges.
      if (!unmatchedIncomingEdges.isEmpty() && error.ok()) {
        // TODO(user): Could return immediately upon finding an error, here and in C++.
        error.init(
            S2Error.Code.BUILDER_EDGES_DO_NOT_FORM_LOOPS,
            "Given edges do not form loops (indegree != outdegree)");
      }
      unmatchedIncomingEdges.clear();
      unmatchedOutgoingEdges.clear();
      v0Edges.clear();
    }
    return error.ok();
  }

  private static void reverse(IntArrayList ial) {
    for (int i = 0, j = ial.size() - 1; i < j; i++, j--) {
      int t = ial.getInt(i);
      ial.set(i, ial.getInt(j));
      ial.set(j, t);
    }
  }

  /**
   * Rotates the edges of "loop" if necessary so that the edge(s) with the largest input edge ids
   * are last. This ensures that when an output loop is equivalent to an input loop, their cyclic
   * edge orders are the same.
   *
   * <p>"minInputIds" is the output of {@link #getMinInputEdgeIds()}, i.e. a list containing the
   * minimum input edge id for every edge id.
   */
  static void canonicalizeLoopOrder(IntArrayList minInputIds, int[] loop) {
    if (loop.length < 2) {
      return;
    }

    // Find the position of the element with the highest input edge id. If there are multiple such
    // elements together (i.e., the edge was split into several pieces by snapping it to several
    // vertices), then we choose the last such position in cyclic order (this attempts to preserve
    // the original loop order even when new vertices are added). For example, if the input edge id
    // sequence is (7, 7, 4, 5, 6, 7) then we would rotate it to obtain (4, 5, 6, 7, 7, 7).

    // The reason that we put the highest-numbered edge last, rather than the lowest-numbered edge
    // first, is that S2Loop.invert() reverses the loop edge order *except* for the last edge. For
    // example, the loop ABCD (with edges AB, BC, CD, DA) becomes DCBA (with edges DC, CB, BA, AD).
    // Note that the last edge is the same except for its direction (DA vs. AD). This has the
    // advantage that if an undirected loop is assembled with the wrong orientation and later
    // inverted (e.g. by S2Polygon.initOriented), we still end up preserving the original cyclic
    // vertex order.
    int pos = 0;
    boolean sawGap = false;
    for (int i = 1; i < loop.length; ++i) {
      int cmp = minInputIds.getInt(loop[i]) - minInputIds.getInt(loop[pos]);
      if (cmp < 0) {
        sawGap = true;
      } else if (cmp > 0 || !sawGap) {
        pos = i;
        sawGap = false;
      }
    }
    if (++pos == loop.length) {
      pos = 0; // Convert loop end to loop start.
    }
    // Rotate all the elements of loop left, so that the element at 'pos' becomes the first element.
    Ints.rotate(loop, -pos);
  }

  /**
   * Sorts the given directed components (i.e. lists of edge chains, which may be loops or
   * polylines) by the minimum input edge id of each component's first chain's first edge.
   */
  static void canonicalizeDirectedComponentOrder(
      IntArrayList minInputIds, List<DirectedComponent> components) {
    Collections.sort(
        components,
        (DirectedComponent componentA, DirectedComponent componentB) -> {
          int minInputEdgeIdFirstLoopEdgeA = minInputIds.getInt(componentA.get(0)[0]);
          int minInputEdgeIdFirstLoopEdgeB = minInputIds.getInt(componentB.get(0)[0]);
          // First sort by minimum input edge id of the components's first loop's first edge.
          if (minInputEdgeIdFirstLoopEdgeA != minInputEdgeIdFirstLoopEdgeB) {
            return Integer.compare(minInputEdgeIdFirstLoopEdgeA, minInputEdgeIdFirstLoopEdgeB);
          }
          // Tie-break with the components's first loops's first edge ids.
          return Integer.compare(componentA.get(0)[0], componentB.get(0)[0]);
        });
  }

  /**
   * Sorts the given edge chains (i.e., loops or polylines) by the minimum input edge id of each
   * chains first edge. This ensures that when the output consists of multiple loops or polylines,
   * they are sorted in the same order as they were provided in the input.
   */
  static void canonicalizeEdgeChainOrder(IntArrayList minInputIds, ArrayList<int[]> chains) {
    Collections.sort(
        chains,
        (int[] a, int[] b) -> {
          int minInputEdgeIdFirstEdgeA = minInputIds.getInt(a[0]);
          int minInputEdgeIdFirstEdgeB = minInputIds.getInt(b[0]);
          // First sort by minimum input edge id of the loop's first edge.
          if (minInputEdgeIdFirstEdgeA != minInputEdgeIdFirstEdgeB) {
            return Integer.compare(minInputEdgeIdFirstEdgeA, minInputEdgeIdFirstEdgeB);
          }
          // Tie-break with the loops first edge id.
          return Integer.compare(a[0], b[0]);
        });
  }

  /**
   * Builds loops from a set of directed edges, turning left at each vertex until either a repeated
   * vertex (for LoopType.SIMPLE) or a repeated edge (for LoopType.CIRCUIT) is found. (Use
   * LoopType.SIMPLE if you intend to construct an S2Loop.)
   *
   * <p>Each loop is represented as a sequence of edges. The edge ordering and loop ordering are
   * automatically canonicalized in order to preserve the input ordering as much as possible. Loops
   * are non-crossing provided that the graph contains no crossing edges. If some edges cannot be
   * turned into loops, returns false and sets "error" appropriately.
   *
   * <p>If any degenerate edges are present, then each such edge is treated as a separate loop. This
   * is mainly useful in conjunction with options.degenerateEdges() == DISCARD_EXCESS, in order to
   * build polygons that preserve degenerate geometry.
   *
   * <p>REQUIRES: options.degenerateEdges() == {DISCARD, DISCARD_EXCESS}
   *
   * <p>REQUIRES: options.edgeType() == DIRECTED
   */
  public boolean getDirectedLoops(LoopType loopType, ArrayList<int[]> loops, S2Error error) {
    Preconditions.checkState(
        options.degenerateEdges() == DegenerateEdges.DISCARD
            || options.degenerateEdges() == DegenerateEdges.DISCARD_EXCESS);
    Preconditions.checkState(options.edgeType() == EdgeType.DIRECTED);

    int[] leftTurnMap = new int[numEdges()];
    if (!fillLeftTurnMap(getInEdgeIds(), leftTurnMap, error)) {
      return false;
    }
    IntArrayList minInputEdgeIds = getMinInputEdgeIds();

    // If LoopType is SIMPLE then we will break loops at repeated vertices. To do this we maintain a
    // map from vertexId to the index of an edge id in "pathEdgeIds" with that vertex as the source
    // vertex.
    int[] pathIndex = new int[numVertices()];
    if (loopType == LoopType.SIMPLE) {
      Arrays.fill(pathIndex, -1);
    }

    // The edge ids of a path through the graph.
    IntArrayList pathEdgeIds = new IntArrayList();

    // Visit edges in arbitrary order, and try to build a loop from each edge.
    for (int startEdgeId = 0; startEdgeId < numEdges(); ++startEdgeId) {
      if (leftTurnMap[startEdgeId] < 0) {
        continue; // We cannot make a left turn from this edge.
      }

      // Build a loop by making left turns at each vertex until we return to "start". We use
      // "leftTurnMap" to keep track of which edges have already been visited by setting its entries
      // to -1 as we go along. If we are building vertex cycles, then whenever we encounter a vertex
      // that is already part of the path, we "peel off" a loop by removing those edges from the
      // path so far.
      int nextEdgeId;

      // Keep making left turns and adding the edge to the path until we cannot make a left turn.
      for (int edgeId = startEdgeId; leftTurnMap[edgeId] >= 0; edgeId = nextEdgeId) {
        pathEdgeIds.add(edgeId);
        nextEdgeId = leftTurnMap[edgeId];
        leftTurnMap[edgeId] = -1; // Mark this edge as visited.

        if (loopType == LoopType.SIMPLE) {
          // Save the location of this edge in the pathIndex for the source vertex.
          pathIndex[edges.getSrcId(edgeId)] = pathEdgeIds.size() - 1;

          // Is this edge's destination vertex already in the pathIndex?
          int loopStart = pathIndex[edges.getDstId(edgeId)];
          if (loopStart < 0) {
            // We haven't seen this destination vertex, keep going.
            continue;
          }

          // This edge ends at the vertex where a previously seen edge starts, forming a loop.
          // Peel off that loop from the path and add it to the output loops.

          // Copy the edge ids that form a loop from the end of the path.
          int[] loop = pathEdgeIds.subList(loopStart, pathEdgeIds.size()).toIntArray();
          // Remove these edges from the path.
          pathEdgeIds.size(loopStart);

          for (int loopEdgeId : loop) {
            // Reset the pathIndex for the vertices in this loop.
            pathIndex[edges.getSrcId(loopEdgeId)] = -1;
          }
          canonicalizeLoopOrder(minInputEdgeIds, loop);
          loops.add(loop);
        }
      }

      if (loopType == LoopType.SIMPLE) {
        // When no more left turns are available and loops have been removed, there must be no more
        // edges in the path.
        Preconditions.checkState(pathEdgeIds.isEmpty());
      } else {
        int[] loop = pathEdgeIds.toIntArray();
        canonicalizeLoopOrder(minInputEdgeIds, loop);
        loops.add(loop);
        pathEdgeIds.clear();
      }
    }
    canonicalizeEdgeChainOrder(minInputEdgeIds, loops);
    return true;
  }

  /**
   * Builds loops from a set of directed edges, turning left at each vertex until a repeated edge is
   * found (i.e., LoopType.CIRCUIT). The loops are further grouped into connected components, where
   * each component consists of one or more loops connected by shared vertices.
   *
   * <p>This method is used to build polygon meshes from directed or undirected input edges. To
   * convert the output of this method into a mesh, the client must determine how the loops in
   * different components are related to each other: for example, several loops from different
   * components may bound the same region on the sphere, in which case all of those loops are
   * combined into a single polygon. See {@link S2LaxPolygonLayer} for details.)
   *
   * <p>Note that loops may include both edges of a sibling pair. When several such edges are
   * connected in a chain or a spanning tree, they form a zero-area "filament". The entire loop may
   * be a filament (i.e., a degenerate loop with an empty interior), or the loop may have have
   * non-empty interior with several filaments that extend inside it, or the loop may consist of
   * several "holes" connected by filaments. These filaments do not change the interior of any loop,
   * so if you are only interested in point containment then they can safely be removed by setting
   * the "degenerateBoundaries" parameter to DISCARD. (They can't be removed by setting
   * (options.siblingPairs() == DISCARD) because the two siblings might belong to different polygons
   * of the mesh.) Note that you can prevent multiple copies of sibling pairs by specifying
   * options.duplicateEdges() == MERGE.
   *
   * <p>Each loop is represented as a sequence of edges. The edge ordering and loop ordering are
   * automatically canonicalized in order to preserve the input ordering as much as possible. Loops
   * are non-crossing provided that the graph contains no crossing edges. If some edges cannot be
   * turned into loops, returns false and sets "error" appropriately.
   *
   * <p>REQUIRES: options.degenerateEdges() == { DISCARD, DISCARD_EXCESS } (but requires DISCARD if
   * degenerateBoundaries == DISCARD) REQUIRES: options.siblingPairs() == { REQUIRE, CREATE } [i.e.,
   * every edge must have a sibling edge]
   */
  public boolean getDirectedComponents(
      DegenerateBoundaries degenerateBoundaries,
      List<DirectedComponent> outputComponents,
      S2Error error) {
    // TODO(torrey): Use the new new index validation system's component detection logic instead of
    // the current algorithm for finding DirectedComponents and UndirectedComponents.
    Preconditions.checkState(
        options.degenerateEdges() == DegenerateEdges.DISCARD
            || (options.degenerateEdges() == DegenerateEdges.DISCARD_EXCESS
                && degenerateBoundaries == DegenerateBoundaries.KEEP));
    Preconditions.checkState(
        options.siblingPairs() == SiblingPairs.REQUIRE
            || options.siblingPairs() == SiblingPairs.CREATE);
    Preconditions.checkState(options.edgeType() == EdgeType.DIRECTED); // Implied by above.

    IntArrayList siblingMap = getSiblingMap();
    int[] leftTurnMap = new int[numEdges()];
    if (!fillLeftTurnMap(siblingMap, leftTurnMap, error)) {
      return false;
    }

    IntArrayList minInputEdgeIds = getMinInputEdgeIds();
    IntArrayList frontier = new IntArrayList(); // Unexplored sibling edge ids.

    // A map from edge id to the position of that edge in "path". Only needed if degenerate
    // boundaries are being discarded.
    int[] pathIndex = new int[numEdges()];
    if (degenerateBoundaries == DegenerateBoundaries.DISCARD) {
      Arrays.fill(pathIndex, -1);
    }

    for (int startEdgeId = 0; startEdgeId < numEdges(); ++startEdgeId) {
      if (leftTurnMap[startEdgeId] < 0) {
        continue; // Already used.
      }

      // Build a connected component by keeping a stack of unexplored siblings of the edges used so
      // far.
      DirectedComponent component = new DirectedComponent();
      frontier.add(startEdgeId);
      while (!frontier.isEmpty()) {
        int edgeId = frontier.popInt();
        if (leftTurnMap[edgeId] < 0) {
          continue; // Already used this edge.
        }

        // Build a path by making left turns at each vertex until we complete a loop. Whenever we
        // encounter an edge that is a sibling of an edge that is already on the path, we
        // "peel off" a loop consisting of any edges that were between these two edges.
        IntArrayList path = new IntArrayList();
        for (int nextEdgeId = 0; leftTurnMap[edgeId] >= 0; edgeId = nextEdgeId) {
          path.add(edgeId);
          nextEdgeId = leftTurnMap[edgeId];
          leftTurnMap[edgeId] = -1;

          // If the sibling hasn't been visited yet, add it to the frontier.
          int siblingEdgeId = siblingMap.getInt(edgeId);
          if (leftTurnMap[siblingEdgeId] >= 0) {
            frontier.add(siblingEdgeId);
          }

          if (degenerateBoundaries == DegenerateBoundaries.DISCARD) {
            // Record the position of this edgeId in the path.
            pathIndex[edgeId] = path.size() - 1;
            int siblingIndex = pathIndex[siblingEdgeId];
            if (siblingIndex < 0) {
              // The sibling edge hasn't been added to the path yet, keep going.
              continue;
            }

            // Common special case: the edge and its sibling are adjacent, in which case we can
            // simply remove them from the path and continue.
            if (siblingIndex == path.size() - 2) {
              path.size(siblingIndex);
              // We don't need to update "pathIndex" for these two edges because both edges of the
              // sibling pair have now been used.
              continue;
            }

            // Peel off a loop from the path.
            int[] loop = path.subList(siblingIndex + 1, path.size() - 2).toIntArray();
            path.size(siblingIndex);

            // Mark the edges that are no longer part of the path.
            for (int loopEdgeId : loop) {
              pathIndex[loopEdgeId] = -1;
            }
            canonicalizeLoopOrder(minInputEdgeIds, loop);
            component.add(loop);
          }
          edgeId = nextEdgeId;
        }

        // Mark the edges that are no longer part of the path.
        if (degenerateBoundaries == DegenerateBoundaries.DISCARD) {
          path.forEach(pathEdgeId -> pathIndex[pathEdgeId] = -1);
        }
        int[] loop = path.toIntArray();
        canonicalizeLoopOrder(minInputEdgeIds, loop);
        component.add(loop);
      }
      canonicalizeEdgeChainOrder(minInputEdgeIds, component);
      outputComponents.add(component);
    }
    canonicalizeDirectedComponentOrder(minInputEdgeIds, outputComponents);
    return true;
  }

  /**
   * There are two complements of each component (a.k.a. the "slot", either 0 or 1). Maps a slot to
   * the other slot.
   */
  private static int markEdgeUsed(int slot) {
    return -1 - slot;
  }

  /**
   * Builds loops from a set of undirected edges, turning left at each vertex until either a
   * repeated vertex (for LoopType.SIMPLE) or a repeated edge (for LoopType.CIRCUIT) is found. The
   * loops are further grouped into "components" such that all the loops in a component are
   * connected by shared vertices. Finally, the loops in each component are divided into two
   * "complements" such that every edge in one complement is the sibling of an edge in the other
   * complement. This corresponds to the fact that given any set of non-crossing undirected loops,
   * there are exactly two possible interpretations of the region that those loops represent (where
   * one possibility is the complement of the other). This method does not attempt to resolve this
   * ambiguity, but instead returns both possibilities for each connected component and lets the
   * client choose among them.
   *
   * <p>This method is used to build single polygons. (Use getDirectedComponents to build polygon
   * meshes, even when the input edges are undirected.) To convert the output of this method into a
   * polygon, the client must choose one complement from each component such that the entire set of
   * loops is oriented consistently (i.e., they define a region such that the interior of the region
   * is always on the left). The non-chosen complements form another set of loops that are also
   * oriented consistently but represent the complementary region on the sphere. Finally, the client
   * needs to choose one of these two sets of loops based on heuristics (e.g., the area of each
   * region), since both sets of loops are equally valid interpretations of the input.
   *
   * <p>Each loop is represented as a sequence of edges. The edge ordering and loop ordering are
   * automatically canonicalized in order to preserve the input ordering as much as possible. Loops
   * are non-crossing provided that the graph contains no crossing edges. If some edges cannot be
   * turned into loops, returns false and sets "error" appropriately.
   *
   * <p>REQUIRES: options.degenerateEdges() == { DISCARD, DISCARD_EXCESS }
   *
   * <p>REQUIRES: options.edgeType() == UNDIRECTED
   *
   * <p>REQUIRES: options.siblingsPairs() == { DISCARD, DISCARD_EXCESS, KEEP } since REQUIRE and
   * CREATE convert the edgeType() to DIRECTED.
   */
  public boolean getUndirectedComponents(
      LoopType loopType, List<UndirectedComponent> outputComponents, S2Error error) {
    Preconditions.checkState(
        options.degenerateEdges() == DegenerateEdges.DISCARD
            || options.degenerateEdges() == DegenerateEdges.DISCARD_EXCESS);
    Preconditions.checkState(options.edgeType() == EdgeType.UNDIRECTED);

    // This is not actually a siblingMap until makeSiblingMap() is called on it. The output of
    // getInEdgeIds is first used to build the left turn map and then converted to a sibling map.
    IntArrayList siblingMap = getInEdgeIds();
    int[] leftTurnMap = new int[numEdges()];
    if (!fillLeftTurnMap(siblingMap, leftTurnMap, error)) {
      return false;
    }

    makeSiblingMap(siblingMap);
    IntArrayList minInputEdgeIds = getMinInputEdgeIds();

    // A stack of unexplored sibling edge ids (the first element of the pair). Each sibling edge id
    // has a "slot" (0 or 1) in the second element of the pair that indicates which of the two
    // complements it belongs to.
    IntPairVector frontier = new IntPairVector();

    // If we are breaking loops at repeated vertices, we maintain a map from VertexId to its
    // position in "path".
    int[] pathIndex = new int[numVertices()];
    if (loopType == LoopType.SIMPLE) {
      Arrays.fill(pathIndex, -1);
    }

    for (int minStartEdgeId = 0; minStartEdgeId < numEdges(); ++minStartEdgeId) {
      if (leftTurnMap[minStartEdgeId] < 0) {
        continue; // Already used.
      }

      // Build a connected component by keeping a stack of unexplored siblings of the edges used
      // so far.
      UndirectedComponent component = new UndirectedComponent();
      frontier.pushPair(minStartEdgeId, 0);
      while (!frontier.isEmpty()) {
        // Take the startEdgeId and its slot from the top of the frontier stack.
        int startEdgeId = frontier.peekFirst();
        int slot = frontier.peekSecond();
        frontier.pop();
        if (leftTurnMap[startEdgeId] < 0) {
          continue; // Already used.
        }

        // Build a path by making left turns at each vertex until we return to "startEdgeId".
        // We use "leftTurnMap" to keep track of which edges have already been visited, and which
        // complement they were assigned to, by setting its entries to negative values as we go
        // along.
        // TODO(torrey): Eric points out that much of the copying here can be avoided. In both the
        // C++ implementation and here, filling in 'path' can be replaced with keeping track of a
        // subset of the leftTurnMap with a pair of integers, since the set of edges to follow is
        // bounded, and the algorithm to read out the leftTurnMap in canonical order is to start at
        // offset0 (which is the position computed in the canonicalize method), iterate through the
        // leftTurnMap up to the position we just found a duplicate vertex in, and then cycle around
        // to the original vertex that was duplicated.
        IntArrayList path = new IntArrayList();
        int nextEdgeId;
        for (int edgeId = startEdgeId; leftTurnMap[edgeId] >= 0; edgeId = nextEdgeId) {
          path.add(edgeId);
          nextEdgeId = leftTurnMap[edgeId];
          leftTurnMap[edgeId] = markEdgeUsed(slot);

          // If the sibling hasn't been visited yet, add it to the frontier.
          int siblingEdgeId = siblingMap.getInt(edgeId);
          if (leftTurnMap[siblingEdgeId] >= 0) {
            frontier.pushPair(siblingEdgeId, 1 - slot);
          } else if (leftTurnMap[siblingEdgeId] != markEdgeUsed(1 - slot)) {
            // Two siblings edges can only belong to the same complement if the given undirected
            // edges do not form loops.
            error.init(
                S2Error.Code.BUILDER_EDGES_DO_NOT_FORM_LOOPS,
                "Given undirected edges do not form loops");
            return false;
          }

          if (loopType == LoopType.SIMPLE) {
            // Whenever we encounter a vertex that is already part of the path, we "peel off" a
            // loop by removing those edges from the path.
            pathIndex[edgeSrcId(edgeId)] = path.size() - 1;
            int loopStart = pathIndex[edgeDstId(edgeId)];
            if (loopStart < 0) {
              continue;
            }
            int[] loop = path.subList(loopStart, path.size()).toIntArray();
            path.size(loopStart);

            // Mark the vertices that are no longer part of the path.
            for (int loopEdgeId : loop) {
              pathIndex[edgeSrcId(loopEdgeId)] = -1;
            }
            canonicalizeLoopOrder(minInputEdgeIds, loop);
            component.addEdgeLoop(slot, loop);
          }
        }
        if (loopType == LoopType.SIMPLE) {
          Preconditions.checkState(path.isEmpty()); // Invariant.
        } else {
          int[] loop = path.toIntArray();
          canonicalizeLoopOrder(minInputEdgeIds, loop);
          component.addEdgeLoop(slot, loop);
        }
      }
      canonicalizeEdgeChainOrder(minInputEdgeIds, component.getComplement(0));
      canonicalizeEdgeChainOrder(minInputEdgeIds, component.getComplement(1));

      // To save some work in S2PolygonLayer, we swap the two complements or loop sets of the
      // component so that the complement whose first loop most closely follows the input edge
      // ordering is first. (If the input was a valid S2Polygon, then this component will contain
      // normalized loops.)
      if (minInputEdgeIds.getInt(component.getComplement(0).get(0)[0])
          > minInputEdgeIds.getInt(component.getComplement(1).get(0)[0])) {
        component.swapComplements();
      }
      outputComponents.add(component);
    }

    // Sort the components to correspond to the input edge ordering.
    Collections.sort(outputComponents, new UndirectedComponentComparator(minInputEdgeIds));
    return true;
  }

  /**
   * A comparator of UndirectedComponents that uses the mapping of edge id to minimum input edge id
   * provided by {@link #getMinInputEdgeIds()} to order UndirectedComponents by the minimum input
   * edge id corresponding to the first edge of the first loop of the complement in slot 0.
   */
  private static final class UndirectedComponentComparator
      implements Comparator<UndirectedComponent> {
    private final IntArrayList minInputEdgeIds;

    public UndirectedComponentComparator(IntArrayList minInputEdgeIds) {
      this.minInputEdgeIds = minInputEdgeIds;
    }

    @Override
    public int compare(UndirectedComponent a, UndirectedComponent b) {
      return Integer.compare(
          minInputEdgeIds.getInt(a.getComplement(0).get(0)[0]),
          minInputEdgeIds.getInt(b.getComplement(0).get(0)[0]));
    }
  }

  /**
   * Builds polylines (either walks or paths) from the edges in an S2BuilderGraph. May be reused by
   * calling {@link #init(S2BuilderGraph)}.
   */
  public static final class PolylineBuilder {
    private S2BuilderGraph graph;
    // TODO(torrey): Make sure PolylineBuilder is actually reused when building layers. Make
    // VertexInMap and VertexOutMap reusable and final. Make the arrays reusable, perhaps by
    // switching to IntArrayLists.
    private VertexInMap inMap;
    private VertexOutMap outMap;

    /** Map from edgeId to sibling edgeId. Only present for graphs with undirected edges. */
    private IntArrayList siblingMap = null;

    /**
     * Maps edge id to the smallest input edge id snapped to that edge. Edges with no input edges
     * snapped to them have NO_INPUT_EDGE_ID.
     */
    private final IntArrayList minInputEdgeIds = new IntArrayList();

    private boolean directed;
    private int edgesLeft;
    // TODO(torrey): Benchmark using a Bitset here.
    private boolean[] used;
    // Map from vertexId to (outDegree(v) - inDegree(v)) considering used edges only.
    private int[] excessUsed;

    /** Constructs a new PolylineBuilder. Before use, init() must be called. */
    public PolylineBuilder() {}

    /** Prepares this PolylineBuilder for building walks or paths on the given graph. */
    public void init(S2BuilderGraph graph) {
      Preconditions.checkNotNull(graph);
      Preconditions.checkState(
          graph.options.siblingPairs() == SiblingPairs.DISCARD
              || graph.options.siblingPairs() == SiblingPairs.DISCARD_EXCESS
              || graph.options.siblingPairs() == SiblingPairs.KEEP);
      this.graph = graph;
      inMap = new VertexInMap(graph);
      outMap = new VertexOutMap(graph);
      graph.fillMinInputEdgeIds(minInputEdgeIds);
      directed = graph.options().edgeType() == EdgeType.DIRECTED;
      edgesLeft = graph.numEdges() / (directed ? 1 : 2);
      used = new boolean[graph.numEdges()];
      Arrays.fill(used, false);
      excessUsed = new int[graph.numVertices()];

      if (!directed) {
        if (siblingMap == null) {
          siblingMap = new IntArrayList();
        } else {
          siblingMap.clear();
        }
        inMap.inEdgeIds.forEach(i -> siblingMap.add(i));
        graph.makeSiblingMap(siblingMap);
      } else {
        // siblingMap is only used for the undirected case.
        if (siblingMap != null) {
          siblingMap.clear();
        }
      }
    }

    /**
     * Builds paths from the edges in the graph. A path is an edge chain with no repeated edges or
     * vertices. Only vertices of indegree and outdegree 1 (or degree 2 in the case of undirected
     * edges) will appear in the interior of paths. This essentially generates one polyline for each
     * edge chain in the graph.
     *
     * <p>This method attempts to preserve the input edge ordering in order to implement
     * idempotency, even when there are repeated edges or loops. This is true whether directed or
     * undirected edges are used. Degenerate edges are also handled appropriately.
     */
    public List<int[]> buildPaths() {
      // First build polylines starting at all the vertices that cannot be in the polyline interior
      // (i.e., indegree != 1 or outdegree != 1 for directed edges, or degree != 2 for undirected
      // edges). We consider the possible starting edges in input edge id order so that we preserve
      // the input path direction even when undirected edges are used. (Undirected edges are
      // represented by sibling pairs where only the edge in the input direction is labeled with an
      // input edge id.)
      // TODO(torrey): This creates an int[] for every path, and there may be a very large number of
      // paths. Also, this does not reuse memory even when the builder and output layer are reused.
      // So, (a) revise this structure, and (b) use instance fields for the memory.
      ArrayList<int[]> polylines = new ArrayList<>();
      IntArrayList edges = graph.getInputEdgeOrder(minInputEdgeIds);
      IntArrayList polyline = new IntArrayList();
      // TODO(torrey): Benchmark visitor vs. IntIterator approach. This is more compact but may be
      // slower. Revise the code in all the relevant places if so.
      edges.forEach(
          edgeId -> {
            if (!used[edgeId] && !isInterior(graph.edgeSrcId(edgeId))) {
              fillPath(edgeId, polyline);
              polylines.add(polyline.toIntArray());
            }
          });
      // If there are any edges left, they form non-intersecting loops. We build each loop and then
      // canonicalize its edge order. We consider candidate starting edges in input edge id order in
      // order to preserve the input direction of undirected loops. Even so, we still need to
      // canonicalize the edge order to ensure that when an input edge is split into an edge chain,
      // the loop does not start in the middle of such a chain.
      for (IntIterator edgeIter = edges.intIterator(); edgesLeft > 0 && edgeIter.hasNext(); ) {
        int edgeId = edgeIter.nextInt();
        if (used[edgeId]) {
          continue;
        }
        fillPath(edgeId, polyline);
        int[] path = polyline.toIntArray();
        canonicalizeLoopOrder(minInputEdgeIds, path);
        polylines.add(path);
      }
      Preconditions.checkState(edgesLeft == 0);

      // Sort the paths to correspond to the input order (if possible).
      canonicalizeEdgeChainOrder(minInputEdgeIds, polylines);
      return polylines;
    }

    /**
     * Builds walks from the edges in the graph. A walk is a chain of edge ids which may pass
     * through the same vertex or even the same edge multiple times (if duplicate edges are
     * present). Each polyline will be as long as possible. This is useful for reconstructing a
     * polyline that has been snapped to a lower resolution, since snapping can cause edges to
     * become identical.
     *
     * <p>This method attempts to preserve the input edge ordering in order to implement
     * idempotency, even when there are repeated edges or loops. This is true whether directed or
     * undirected edges are used. Degenerate edges are also handled appropriately.
     */
    public List<int[]> buildWalks() {
      // Note that some of this code is worst-case quadratic in the maximum vertex degree. This
      // could be fixed with a few extra arrays, but it should not be a problem in practice.

      // First, build polylines from all vertices where outdegree > indegree (or for undirected
      // edges, vertices whose degree is odd). We consider the possible starting edges in input
      // edge id order, for idempotency in the case where multiple input polylines share vertices
      // or edges.
      ArrayList<IntArrayList> polylines = new ArrayList<>();
      IntArrayList edges = graph.getInputEdgeOrder(minInputEdgeIds);
      for (IntIterator edgeIter = edges.intIterator(); edgeIter.hasNext(); ) {
        int edgeId = edgeIter.nextInt();
        if (used[edgeId]) {
          continue;
        }

        int vertexId = graph.edgeSrcId(edgeId);
        int excess = excessDegree(vertexId); // Out-edges minus in-edges.
        if (excess <= 0) {
          continue;
        }
        excess -= excessUsed[vertexId];
        if (directed ? (excess <= 0) : (excess % 2 == 0)) {
          continue;
        }

        // There is at least one unused out-edge at vertexId. Build a polyline starting there.
        IntArrayList polyline = new IntArrayList();
        // At the start vertex, only an out-edge will be used by the new polyline.
        ++excessUsed[vertexId]; // outdegree - indegree of used edges
        // For all the interior vertices, equal number of in and out edges are used.
        fillWalk(vertexId, polyline);
        polylines.add(polyline);
        // At the last vertex of the polyline, only an in-edge is used.
        --excessUsed[graph.edgeDstId(polyline.getInt(polyline.size() - 1))];
      }

      // Now all vertices have outdegree == indegree (or even degree if undirected edges are being
      // used). Therefore all remaining edges can be assembled into loops. We first try to expand
      // the existing polylines if possible by adding loops to them.
      if (edgesLeft > 0) {
        for (IntArrayList polyline : polylines) {
          maximizeWalk(polyline);
        }
      }

      // Convert the IntArrayLists produced so far to arrays of ints.
      ArrayList<int[]> walks = new ArrayList<>();
      for (IntArrayList polyline : polylines) {
        walks.add(polyline.toIntArray());
      }
      polylines.clear();

      // Finally, if there are still unused edges then we build loops. If the input is a polyline
      // that forms a loop, then for idempotency we need to start from the edge with minimum input
      // edge id. If the minimal input edge was split into several edges, then we start from the
      // first edge of the chain.
      IntArrayList polyline = new IntArrayList(); // Reused for each loop

      // Find every unused edge in the graph, in input edge order.
      for (int i = 0; i < edges.size() && edgesLeft > 0; ++i) {
        int edgeId = edges.getInt(i);
        if (used[edgeId]) {
          continue;
        }

        // Determine whether the origin of this edge is the start of an edge chain. To do this, we
        // test whether (outdegree - indegree == 1) for the origin, considering only unused edges
        // with the same minimum input edge id. (Undirected edges have input edge ids in one
        // direction only.)
        int vertexId = graph.edgeSrcId(edgeId);
        int inputEdgeId = minInputEdgeIds.getInt(edgeId);
        int excess = 0; // outdegree - indegree at the origin
        for (int j = i;
            j < edges.size() && minInputEdgeIds.getInt(edges.getInt(j)) == inputEdgeId;
            ++j) {
          int edgeId2 = edges.getInt(j);
          if (used[edgeId2]) {
            continue;
          }
          if (graph.edgeSrcId(edgeId2) == vertexId) {
            ++excess;
          }
          if (graph.edgeDstId(edgeId2) == vertexId) {
            --excess;
          }
        }

        // If there is one excess out-edge, build a loop from this vertex.
        // It is also acceptable to start a polyline from any degenerate edge.
        if (excess == 1 || graph.edgeDstId(edgeId) == vertexId) {
          fillWalk(vertexId, polyline);
          maximizeWalk(polyline);
          walks.add(polyline.toIntArray());
        }
      }
      Preconditions.checkState(edgesLeft == 0);

      // Sort the walks to correspond to the input order (if possible).
      canonicalizeEdgeChainOrder(minInputEdgeIds, walks);
      return walks;
    }

    /**
     * True if vertexId is part of the interior of an edge chain. For directed graphs, this means
     * the outDegree and indegree are both 1. For undirected graphs, the degree is 2.
     */
    private boolean isInterior(int vertexId) {
      return directed
          ? inMap.inDegree(vertexId) == 1 && outMap.outDegree(vertexId) == 1
          : outMap.outDegree(vertexId) == 2;
    }

    /**
     * For directed graphs, the outdegree - indegree at vertexId. For undirected graphs, 1 if the
     * outdegree is odd, 0 if the outdegree is even.
     */
    private int excessDegree(int vertexId) {
      return directed
          ? outMap.outDegree(vertexId) - inMap.inDegree(vertexId)
          : outMap.outDegree(vertexId) % 2;
    }

    /** Clears the provided polyline and fills it with a path starting with the provided edgeId. */
    private void fillPath(int edgeId, IntArrayList polyline) {
      polyline.clear();
      // We simply follow edges until either we reach a vertex where there is a choice about which
      // way to go (where isInterior(v) is false), or we return to the starting vertex (if the
      // polyline is actually a loop).
      int startVertexId = graph.edgeSrcId(edgeId);
      for (; ; ) {
        polyline.add(edgeId);
        Preconditions.checkState(!used[edgeId]);
        used[edgeId] = true;
        if (!directed) {
          used[siblingMap.getInt(edgeId)] = true;
        }
        --edgesLeft;
        int vertexId = graph.edgeDstId(edgeId);
        if (!isInterior(vertexId) || vertexId == startVertexId) {
          break;
        }
        if (directed) {
          Preconditions.checkState(outMap.outDegree(vertexId) == 1);
          edgeId = outMap.edgeIds(vertexId).begin();
        } else {
          Preconditions.checkState(outMap.outDegree(vertexId) == 2);
          for (OfInt iter = outMap.edgeIds(vertexId).intIterator(); iter.hasNext(); ) {
            int edgeId2 = iter.nextInt();
            if (!used[edgeId2]) {
              edgeId = edgeId2;
            }
          }
        }
      }
    }

    /** Clears the provided polyline and fills it with a walk starting at the provided vertexId. */
    private void fillWalk(int vertexId, IntArrayList polyline) {
      polyline.clear();
      for (; ; ) {
        // Follow the best unused out-edge from vertexId with the smallest input edge id. For
        // undirected graphs, automatically created edges have no input edges snapped to them, so
        // have minInputEdgeId == NO_INPUT_EDGE_ID and thus will be used last.
        int bestEdgeId = -1;
        int bestOutId = Integer.MAX_VALUE; // Higher than NO_INPUT_EDGE_ID.
        for (OfInt outEdgeIter = outMap.edgeIds(vertexId).intIterator(); outEdgeIter.hasNext(); ) {
          int edgeId = outEdgeIter.nextInt();
          if (used[edgeId]) {
            continue;
          }
          int minInputEdgeId = minInputEdgeIds.getInt(edgeId);
          if (minInputEdgeId >= bestOutId) {
            continue;
          }
          bestOutId = minInputEdgeId;
          bestEdgeId = edgeId;
        }
        if (bestEdgeId < 0) {
          return; // No unused out-edge was found.
        }

        // For idempotency when there are multiple input polylines, we stop the walk early if
        // "bestEdgeId" might be a continuation of a different incoming edge.
        int excessUnused = excessDegree(vertexId) - excessUsed[vertexId];
        if (directed ? (excessUnused < 0) : (excessUnused % 2) == 1) {
          for (OfInt inEdgeIter = inMap.edgeIds(vertexId).intIterator(); inEdgeIter.hasNext(); ) {
            int inEdgeId = inEdgeIter.nextInt();
            if (!used[inEdgeId] && minInputEdgeIds.getInt(inEdgeId) <= bestOutId) {
              return;
            }
          }
        }

        // Extend the walk.
        polyline.add(bestEdgeId);
        used[bestEdgeId] = true;
        if (!directed) {
          used[siblingMap.getInt(bestEdgeId)] = true;
        }
        --edgesLeft;
        vertexId = graph.edgeDstId(bestEdgeId);
      }
    }

    /**
     * Examine all vertices of the polyline and check whether there are any unused outgoing edges.
     * If so, then build a loop starting at that vertex and insert it into the polyline. (The walk
     * is guaranteed to be a loop because this method is only called when all vertices have equal
     * numbers of unused incoming and outgoing edges.)
     */
    private void maximizeWalk(IntArrayList polyline) {
      IntArrayList loop = new IntArrayList(); // Reused for each inserted loop
      for (int i = 0; i <= polyline.size(); ++i) {
        int vertexId =
            i == 0 ? graph.edgeSrcId(polyline.getInt(i)) : graph.edgeDstId(polyline.getInt(i - 1));
        for (OfInt outEdgeIter = outMap.edgeIds(vertexId).intIterator(); outEdgeIter.hasNext(); ) {
          int edgeId = outEdgeIter.nextInt();
          if (!used[edgeId]) {
            fillWalk(vertexId, loop);
            Preconditions.checkState(vertexId == graph.edgeDstId(loop.getInt(loop.size() - 1)));
            polyline.addAll(i, loop);
            Preconditions.checkState(used[edgeId]); // All outgoing edges from "v" are now used.
            break;
          }
        }
      }
    }
  }

  /**
   * Constructs a new graph with the given GraphOptions and containing the given edges. Each edge is
   * associated with a (possibly empty) set of input edges as specified by newInputEdgeIdSetIds
   * (which must be the same length as "newEdges") and the given IdSetLexicon (which allows looking
   * up the set of input edge ids associated with a graph edge). Finally, the subgraph may also
   * specify a new IsFullPolygonPredicate (which is used to distinguish an empty polygon possibly
   * with degenerate shells from a full polygon possibly with degenerate holes).
   *
   * <p>The output graph has the same set of vertices and edge labels as the input graph (noting
   * that edge labels are associated with *input* edges, not graph edges).
   *
   * <p>If newOptions.edgeType() is UNDIRECTED then note the following:
   *
   * <ul>
   *   <li>If this.options().edgeType() is DIRECTED then each input edge will be transformed into a
   *       pair of directed edges before calling processEdges() above.
   *   <li>If newOptions.siblingPairs() is CREATE or REQUIRE then processEdges() will discard half
   *       of the edges in each direction and change edgeType() to DIRECTED. See {@link
   *       S2Builder.GraphOptions.SiblingPairs}.
   * </ul>
   */
  @SuppressWarnings("ParameterName") // Deliberately switching src and dst to reverse the edge.
  public S2BuilderGraph makeSubgraph(
      GraphOptions newOptions,
      EdgeList newEdges,
      IntArrayList newInputEdgeIdSetIds,
      IdSetLexicon newInputEdgeIdSetLexicon,
      IsFullPolygonPredicate isFullPolygonPredicate,
      S2Error error) {
    // If the subgraph is undirected but the graph is directed, create a reversed edge for every
    // edge in newEdges.
    if (options().edgeType() == EdgeType.DIRECTED && newOptions.edgeType() == EdgeType.UNDIRECTED) {
      int n = newEdges.size();
      newEdges.ensureCapacity(2 * n);
      newInputEdgeIdSetIds.ensureCapacity(2 * n);

      for (int i = 0; i < n; ++i) {
        // Add a reversed copy of the edge.
        newEdges.add(/* srcId= */ newEdges.getDstId(i), /* dstId= */ newEdges.getSrcId(i));
        newInputEdgeIdSetIds.add(IdSetLexicon.EMPTY_SET_ID);
      }
    }
    processEdges(newOptions, newEdges, newInputEdgeIdSetIds, newInputEdgeIdSetLexicon, error);

    return new S2BuilderGraph(
        newOptions,
        vertices,
        newEdges,
        newInputEdgeIdSetIds,
        newInputEdgeIdSetLexicon,
        labelSetIds,
        labelSetLexicon,
        isFullPolygonPredicate);
  }

  /**
   * Given an unsorted collection of edges, transform them according to the given set of
   * GraphOptions. This includes actions such as discarding degenerate edges; merging duplicate
   * edges; and canonicalizing sibling edge pairs in several possible ways (e.g. discarding or
   * creating them). The output is suitable for passing to the Graph constructor.
   *
   * <p>All of the provided parameters are modified. The 'edges', 'inputIds', and 'idSetLexicon' are
   * modified by the graph transformation. The provided 'error' may be set if the graph is invalid,
   * such as missing sibling pairs when they are required by the options.
   *
   * <p>If options.edgeType() == EdgeType.UNDIRECTED, then all input edges should already have been
   * transformed into a pair of directed edges.
   *
   * <p>"inputIds" is a vector of the same length as "edges" that indicates which input edges were
   * snapped to each edge, by mapping edgeId to a set of input edge ids in the idSetLexicon. The
   * "inputIds" vector and "idSetLexicon" are updated appropriately as edges are discarded, merged,
   * etc.
   *
   * <p>Note that "options" may be modified by this method: in particular, if edgeType() is
   * UNDIRECTED and siblingPairs() is CREATE or REQUIRE, then half of the edges in each direction
   * will be discarded and edgeType() will be changed to DIRECTED. See {@link
   * S2Builder.GraphOptions.SiblingPairs}.
   */
  public static void processEdges(
      GraphOptions options,
      EdgeList edges,
      IntArrayList inputIds,
      IdSetLexicon idSetLexicon,
      S2Error error) {
    EdgeProcessor processor = new EdgeProcessor(options, edges, inputIds, idSetLexicon);
    processor.run(error);

    // Certain values of siblingPairs() discard half of the edges and change the edgeType() to
    // DIRECTED. See the description in {@link S2Builder.GraphOptions.SiblingPairs}.
    if (options.siblingPairs() == SiblingPairs.REQUIRE
        || options.siblingPairs() == SiblingPairs.CREATE) {
      options.setEdgeType(EdgeType.DIRECTED);
    }
  }

  /**
   * Given a list of vertices and edges, removes all vertices that do not have any edges and returns
   * the new, minimal set of vertices. Also updates each edge in "edges" to correspond to the new
   * vertex numbering. (Note that this method does *not* merge duplicate vertices, it simply removes
   * vertices of degree zero.)
   *
   * <p>The new vertex ordering is a subsequence of the original ordering, therefore if the edges
   * were lexicographically sorted before calling this method then they will still be sorted after
   * calling this method.
   *
   * <p>The "vmap" and "usedVertexIds" are temporary storage used by this method. All calls to this
   * method from a single thread can reuse the same temporary storage. They should initially point
   * to empty vectors. This can make a big difference to efficiency when this method is called many
   * times (e.g. to extract the vertices for different layers), since the incremental running time
   * for each layer becomes O(edges.size()) rather than O(vertices.size() + edges.size()).
   */
  public static ArrayList<S2Point> filterVertices(
      List<S2Point> vertices, EdgeList edges, IntVector vmap, IntVector usedVertexIds) {
    // Gather the vertices that are actually used by edges.
    usedVertexIds.clear();
    usedVertexIds.ensureCapacity(2 * edges.size());
    int lastDstId = -1;
    for (int e = 0; e < edges.size(); e++) {
      int edgeSrcId = edges.getSrcId(e);
      if (edgeSrcId != lastDstId) {
        usedVertexIds.add(edgeSrcId); // Cheap but effective partial deduplication.
      }
      lastDstId = edges.getDstId(e);
      usedVertexIds.add(lastDstId);
    }
    // Sort the vertex ids and find the distinct ones.
    usedVertexIds.sort();
    usedVertexIds.unique();

    // Build the list of new vertices, and generate a map from old vertex id to new vertex id.
    vmap.resize(vertices.size());
    // TODO(torrey): Consider providing a view of the filtered vertices, instead of copying them.
    // Benchmark first.
    ArrayList<S2Point> newVertices = new ArrayList<>(usedVertexIds.size());
    for (int i = 0; i < usedVertexIds.size(); ++i) {
      newVertices.add(vertices.get(usedVertexIds.get(i)));
      vmap.set(usedVertexIds.get(i), i);
    }

    // Update the edges with the new vertex ids.
    for (int e = 0; e < edges.size(); e++) {
      edges.set(e, vmap.get(edges.getSrcId(e)), vmap.get(edges.getDstId(e)));
    }
    return newVertices;
  }

  /**
   * Returns an IntComparator that compares EdgeIds in the given EdgeList by source vertex id,
   * breaking ties with destination vertex ids, and breaking ties for equal edges with edge ids.
   */
  private static IntComparator bySrcVertex(EdgeList edges) {
    return new IntComparator() {
      @Override
      public int compare(int edgeIdA, int edgeIdB) {
        int srcVertexIdA = edges.getSrcId(edgeIdA);
        int srcVertexIdB = edges.getSrcId(edgeIdB);
        if (srcVertexIdA != srcVertexIdB) {
          return Integer.compare(srcVertexIdA, srcVertexIdB);
        }
        int dstVertexIdA = edges.getDstId(edgeIdA);
        int dstVertexIdB = edges.getDstId(edgeIdB);
        if (dstVertexIdA != dstVertexIdB) {
          return Integer.compare(dstVertexIdA, dstVertexIdB);
        }
        return Integer.compare(edgeIdA, edgeIdB);
      }
    };
  }

  /**
   * Returns an IntComparator that compares EdgeIds in the given EdgeList by destination vertex id,
   * breaking ties with source vertex ids, and breaking ties for equal edges with edge ids.
   */
  private static IntComparator byDstVertex(EdgeList edges) {
    return (int edgeIdA, int edgeIdB) -> {
      int dstVertexIdA = edges.getDstId(edgeIdA);
      int dstVertexIdB = edges.getDstId(edgeIdB);
      if (dstVertexIdA != dstVertexIdB) {
        return Integer.compare(dstVertexIdA, dstVertexIdB);
      }
      int srcVertexIdA = edges.getSrcId(edgeIdA);
      int srcVertexIdB = edges.getSrcId(edgeIdB);
      if (srcVertexIdA != srcVertexIdB) {
        return Integer.compare(srcVertexIdA, srcVertexIdB);
      }
      return Integer.compare(edgeIdA, edgeIdB);
    };
  }

  /** A simple container for source and destination vertex ids defining an edge. */
  public static class Edge {
    int srcId;
    int dstId;

    /** Constructs and returns a new Edge with the given srcId and dstId. */
    public static Edge of(int srcId, int dstId) {
      Edge e = new Edge();
      e.srcId = srcId;
      e.dstId = dstId;
      return e;
    }

    /** Sets srcId and dstId to the given values. */
    public void set(int src, int dst) {
      this.srcId = src;
      this.dstId = dst;
    }

    /** Sets srcId and dstId to the values on the other Edge. */
    public void set(Edge other) {
      this.srcId = other.srcId;
      this.dstId = other.dstId;
    }

    /** Sets srcId and dstId to SENTINEL_ID. */
    public void setSentinel() {
      this.srcId = SENTINEL_ID;
      this.dstId = SENTINEL_ID;
    }

    /** Returns true if both srcId and dstId == SENTINEL_ID. */
    public boolean isSentinel() {
      return srcId == SENTINEL_ID && dstId == SENTINEL_ID;
    }

    /** Returns true if this Edge is lexicographically less than 'other'. */
    public boolean isLessThan(Edge other) {
      if (srcId == other.srcId) {
        return dstId < other.dstId;
      }
      return srcId < other.srcId;
    }

    /** Returns true if this.srcId == other.srcId and this.dstId == other.dstId. */
    public boolean isEqualTo(Edge other) {
      return srcId == other.srcId && dstId == other.dstId;
    }

    @Override
    public boolean equals(Object o) {
      return o instanceof Edge && isEqualTo((Edge) o);
    }

    @Override
    public int hashCode() {
      return srcId * 39119 + dstId;
    }

    @Override
    public String toString() {
      return Platform.formatString("%d-%d", srcId, dstId);
    }
  }

  /**
   * An EdgeList is like an ArrayList of pairs of source and destination vertex ids, but is more
   * efficient by avoiding boxing. EdgeList also defines the lexicographic ordering of edges,
   * sorting first by source vertex id and then destination vertex id.
   */
  public static class EdgeList extends IntPairVector {
    /** Constructs an empty EdgeList. */
    public EdgeList() {
      super();
    }

    /** Copy constructor. */
    public EdgeList(EdgeList other) {
      super(other);
    }

    /** Factory method for creating an EdgeList containing the provided Edges. For testing. */
    public static EdgeList of(Edge... edges) {
      EdgeList list = new EdgeList();
      for (Edge e : edges) {
        list.add(e);
      }
      return list;
    }

    /**
     * Compares two edges 'left' and 'right' defined by their source and destination vertex ids.
     * Returns 0 if the two edges are equal, -1 if 'left' is ordered before 'right', and 1 if 'left'
     * is ordered after 'right'.
     */
    public static int compareEdges(Edge left, Edge right) {
      if (left.srcId == right.srcId) {
        return Integer.compare(left.dstId, right.dstId);
      }
      return Integer.compare(left.srcId, right.srcId);
    }

    /**
     * Compares two edges 'left' and 'right' defined by their source and destination vertex ids.
     * Returns 0 if the two edges are equal, -1 if 'left' is ordered before 'right', and 1 if 'left'
     * is ordered after 'right'.
     */
    public static int compareEdges(int leftSrcId, int leftDstId, int rightSrcId, int rightDstId) {
      if (leftSrcId == rightSrcId) {
        return Integer.compare(leftDstId, rightDstId);
      }
      return Integer.compare(leftSrcId, rightSrcId);
    }

    /** Adds the given Edge to this EdgeList. */
    public void add(Edge e) {
      add(e.srcId, e.dstId);
    }

    /** Copies the edge from the given EdgeList at index {@code other[index]} to this EdgeList. */
    // TODO(torrey): Rename this 'addEdgeFrom'
    public void copyEdge(EdgeList other, int index) {
      add(other.getSrcId(index), other.getDstId(index));
    }

    /** Returns the source vertex id of the edge at the given index. */
    public int getSrcId(int index) {
      return getFirst(index);
    }

    /** Returns the destination vertex id of the edge at the given index. */
    public int getDstId(int index) {
      return getSecond(index);
    }

    /**
     * Fills the given Edge with the source and destination vertex ids of the edge at the given
     * index.
     */
    public void get(int index, Edge edge) {
      edge.set(getFirst(index), getSecond(index));
    }

    /**
     * Fills the given Edge with the reversed source and destination vertex ids of the edge at the
     * given index.
     */
    public void getReverse(int index, Edge edge) {
      edge.set(getSecond(index), getFirst(index));
    }

    /**
     * Replaces the edge at the specified 'index' with the given source and destination vertex ids.
     */
    public void set(int index, int srcId, int dstId) {
      setPair(index, srcId, dstId);
    }

    /**
     * Returns true if this EdgeList contains the sibling of the edge at id 'index'. Requires that
     * this EdgeList is sorted as defined by compareEdges() above. A degenerate edge is its own
     * sibling.
     */
    public boolean containsSibling(int index) {
      return contains(getDstId(index), getSrcId(index));
    }
  }

  /**
   * VertexEdgeList is like an ArrayList of "VertexEdges", where a "VertexEdge" represents an
   * incoming or outgoing edge to or from some common vertex "v0". These can be sorted in clockwise
   * order. However, VertexEdge isn't actually an interface, and VertexEdges are never instantiated.
   */
  public static class VertexEdgeList implements Sorter.SortableCollection {
    /**
     * For each edge 'e' at v0, incomings.get(e) is 1 if the edge is incoming to "v0", 0 if the edge
     * is outgoing.
     */
    private final IntVector incomings = new IntVector();

    /**
     * For each edge 'e' at v0, indexEdgeIds.get(e) is the index of the edge in the graph edges or
     * inEdgeIds.
     */
    private final IntVector indexEdgeIds = new IntVector();

    /**
     * For each edge 'e' at v0, endpointVertexIds.get(e) is the vertex id of the other end of the
     * edge.
     */
    private final IntVector endpointVertexIds = new IntVector();

    /**
     * For each edge 'e' at v0, ranks.get(e) is a rank assigned by {@link #addVertexEdges(int, int,
     * int, int, int, VertexEdgeList)} to provide a secondary key for ordering edges that have the
     * same endpoints.
     */
    private final IntVector ranks = new IntVector();

    /**
     * While sorting the VertexEdgeList in clockwise order, 'v0' is the vertex id of the common
     * endpoint.
     */
    private int v0;

    /**
     * While sorting the VertexEdgeList in clockwise order, 'minEndpointVertexId' is the vertex id
     * of the first edge endpoint.
     */
    private int minEndpointVertexId;

    /** A reference to the list of vertices in the graph that the edges are part of. */
    private final List<S2Point> vertices;

    /** Constructs a new VertexEdgeList using the given graph vertices. */
    public VertexEdgeList(List<S2Point> vertices) {
      this.vertices = vertices;
    }

    /** The number of VertexEdges in this VertexEdgeList. */
    @Override
    public int size() {
      return incomings.size();
    }

    @Override
    public void truncate(int start) {
      incomings.truncate(start);
      indexEdgeIds.truncate(start);
      endpointVertexIds.truncate(start);
      ranks.truncate(start);
    }

    /** True if this VertexEdgeList contains no VertexEdges. */
    public boolean isEmpty() {
      return incomings.isEmpty();
    }

    /** Clears the content of this VertexEdgeList, making it empty. */
    public void clear() {
      incomings.clear();
      indexEdgeIds.clear();
      endpointVertexIds.clear();
      ranks.clear();
    }

    /** Gets the 'incoming' value of the VertexEdge at the provided 'index'. */
    public boolean incoming(int index) {
      return incomings.get(index) == 1;
    }

    /** Gets the 'indexEdgeId' value of the VertexEdge at the provided 'index'. */
    public int indexEdgeId(int index) {
      return indexEdgeIds.get(index);
    }

    /** Gets the 'endpointVertexId' value of the VertexEdge at the provided 'index'. */
    public int endpointVertexId(int index) {
      return endpointVertexIds.get(index);
    }

    /** Gets the 'rank' of the VertexEdge at the provided 'index'. */
    public int rank(int index) {
      return ranks.get(index);
    }

    /** Adds a VertexEdge to the VertexEdgeList with the given VertexEdge parameters. */
    public void add(boolean incoming, int indexEdgeId, int endpointVertexId, int rank) {
      incomings.add(incoming ? 1 : 0);
      indexEdgeIds.add(indexEdgeId);
      endpointVertexIds.add(endpointVertexId);
      ranks.add(rank);
    }

    /** Sorts the VertexEdges around the vertex with id 'v0' in clockwise order. */
    public void sortClockwiseAround(int v0) {
      this.v0 = v0;
      this.minEndpointVertexId = endpointVertexIds.get(0);
      // Edge 0 defines the minimum, so only edges from 1 need to be sorted. The Sorter orders Edges
      // using the less() method below.
      Sorter.sort(this, 1, incomings.size() - 1);
    }

    /** Exchanges the VertexEdges at the given indices 'a' and 'b'. */
    @Override
    public void swap(int a, int b) {
      incomings.swap(a, b);
      indexEdgeIds.swap(a, b);
      endpointVertexIds.swap(a, b);
      ranks.swap(a, b);
    }

    /** True if the VertexEdge at 'a' comes before the VertexEdge at 'b' in a clockwise order. */
    @Override
    public boolean less(int a, int b) {
      if (endpointVertexIds.get(a) == endpointVertexIds.get(b)) {
        // The edges are duplicates, so order by rank.
        return ranks.get(a) < ranks.get(b);
      }
      if (endpointVertexIds.get(a) == minEndpointVertexId) {
        // Edge 'a' has the minimum endpoint so must precede 'b'.
        return true;
      }
      if (endpointVertexIds.get(b) == minEndpointVertexId) {
        // Edge 'b' has the minimum endpoint so must precede 'a'.
        return false;
      }

      // Otherwise, use orderedCCW() to determine if the edges (v0, a.endpoint), (v0, b.endpoint),
      // and (v0, minEndpoint) are encountered in that order when sweeping counter-clockwise
      // around v0. If not, 'a' is less than 'b'.
      return !S2Predicates.orderedCCW(
          vertices.get(endpointVertexIds.get(a)),
          vertices.get(endpointVertexIds.get(b)),
          vertices.get(minEndpointVertexId),
          vertices.get(v0));
    }
  }

  /**
   * Given a set of duplicate outgoing edges (v0, v1) and a set of duplicate incoming edges (v1,
   * v0), this method assigns each edge an integer "rank" so that the edges are sorted in a
   * consistent order with respect to their orderings around "v0" and "v1". Usually there is just
   * one edge, in which case this is easy. Sometimes there is one edge in each direction, in which
   * case the outgoing edge is always ordered before the incoming edge.
   *
   * <p>In general, we allow any number of duplicate edges in each direction, in which case outgoing
   * edges are interleaved with incoming edges so as to create as many degenerate (two-edge) loops
   * as possible. In order to get a consistent ordering around "v0" and "v1", we move forwards
   * through the list of outgoing edges and backwards through the list of incoming edges. If there
   * are more incoming edges, they go at the beginning of the ordering, while if there are more
   * outgoing edges then they go at the end.
   *
   * <p>For example, suppose there are 2 edges "a,b" from "v0" to "v1", and 4 edges "w,x,y,z" from
   * "v1" to "v0". Using lower/upper case letters to represent incoming/outgoing edges, the
   * clockwise ordering around v0 would be zyAxBw, and the clockwise ordering around v1 would be
   * WbXaYZ. (Try making a diagram with each edge as a separate arc.)
   */
  private static void addVertexEdges(
      int outBeginEdgeId,
      int outEndEdgeId,
      int inBeginEdgeId,
      int inEndEdgeId,
      int v1,
      VertexEdgeList v0Edges) {
    int rank = 0;
    // Any extra incoming edges go at the beginning of the ordering.
    while (inEndEdgeId - inBeginEdgeId > outEndEdgeId - outBeginEdgeId) {
      v0Edges.add(true, --inEndEdgeId, v1, rank++);
    }
    // Next we interleave as many outgoing and incoming edges as possible.
    while (inEndEdgeId > inBeginEdgeId) {
      v0Edges.add(false, outBeginEdgeId++, v1, rank++);
      v0Edges.add(true, --inEndEdgeId, v1, rank++);
    }
    // Any extra outgoing edges to at the end of the ordering.
    while (outEndEdgeId > outBeginEdgeId) {
      v0Edges.add(false, outBeginEdgeId++, v1, rank++);
    }
  }

  /** A helper class for VertexOutMap that represents outgoing edges from a given vertex. */
  public static class VertexOutEdges {
    // A shared reference to the S2BuilderGraph's EdgeList, which is sorted in lexicographic order
    // by source vertex, then destination vertex, so that all the outgoing edges from a given vertex
    // are in consecutive order.
    private final EdgeList edges;

    // The index of the first edge in 'edges' from the given vertex (inclusive).
    private final int begin;

    // The index of the last edge in 'edges' from the given vertex (exclusive). Must be >= begin.
    // If end == begin there are no outgoing edges.
    private final int end;

    /** Constructs a new VertexOutEdges as a begin-end range over the given EdgeList. */
    VertexOutEdges(EdgeList edges, int begin, int end) {
      Preconditions.checkArgument(begin >= 0, "begin: %s", begin);
      Preconditions.checkArgument(end >= begin, "end: %s, begin: %s", end, begin);
      Preconditions.checkArgument(
          end <= edges.size(), "end: %s, edges.size(): %s", end, edges.size());
      this.edges = edges;
      this.begin = begin;
      this.end = end;
    }

    /** The number of out-edges from the given vertex. */
    public int size() {
      return end - begin;
    }

    /**
     * Returns the source vertex id of the outgoing edge from the given vertex at 'index', which
     * must be between 0 and {@link #size()}. (This will always be the given vertex.)
     */
    public int getSrcId(int index) {
      Preconditions.checkArgument(index >= 0 && index < size());
      return edges.getSrcId(begin + index);
    }

    /** Returns the destination vertex id of the edge at the given index. */
    public int getDstId(int index) {
      Preconditions.checkArgument(index >= 0 && index < size());
      return edges.getDstId(begin + index);
    }

    /**
     * Fills the given Edge with the source and destination vertex ids of the outgoing edge from the
     * given vertex at the given index.
     */
    public void get(int index, Edge edge) {
      Preconditions.checkArgument(index >= 0 && index < size());
      edge.set(edges.getSrcId(begin + index), edges.getDstId(begin + index));
    }

    /**
     * Fills the given Edge with the reversed source and destination vertex ids of the outgoing edge
     * from the given vertex at the given index.
     */
    public void getReverse(int index, Edge edge) {
      Preconditions.checkArgument(index >= 0 && index < size());
      edge.set(edges.getDstId(begin + index), edges.getSrcId(begin + index));
    }
  }

  // TODO(torrey): VertexOutEdgeIds might be generally useful. Consider renaming it RangeSequence or
  // Range and moving it into the primitives package? Or, it could just be IntStream.range().
  /**
   * A helper class for VertexOutMap that represents the set of outgoing edge *ids* from a given
   * vertex, VertexOutEdgeIds is simply a sequence of ints from 'begin' to 'end', with additional
   * accessors begin() and end().
   *
   * <p>An instance of VertexOutEdgeIds represents the ids of the first (inclusive) and last
   * (exclusive) edges in the Edges list that originate at some vertex. As edges are stored in
   * lexicographic order, all the outgoing Edges from a given source vertex are stored consecutively
   * in the Edges list.
   */
  public static class VertexOutEdgeIds implements IntSequence {
    private final int beginEdgeId;
    private final int endEdgeId;

    /** Creates a new VertexOutEdgeIds from beginEdgeId (inclusive) to endEdgeId (exclusive). */
    VertexOutEdgeIds(int beginEdgeId, int endEdgeId) {
      this.beginEdgeId = beginEdgeId;
      this.endEdgeId = endEdgeId;
    }

    @Override
    public int size() {
      return endEdgeId - beginEdgeId;
    }

    public int begin() {
      return beginEdgeId;
    }

    public int end() {
      return endEdgeId;
    }

    @Override
    public OfInt intIterator() {
      return new OfInt() {
        int next = beginEdgeId;

        @Override
        public boolean hasNext() {
          return next < endEdgeId;
        }

        @Override
        public int nextInt() {
          return next++;
        }
      };
    }

    @Override
    public void forEach(IntConsumer action) {
      for (int i = beginEdgeId; i < endEdgeId; i++) {
        action.accept(i);
      }
    }

    @Override
    public void forEach(IntBiConsumer action) {
      for (int i = beginEdgeId; i < endEdgeId; i++) {
        action.accept(i - beginEdgeId, i);
      }
    }
  }

  /**
   * A class that maps vertices to their outgoing edge ids. Example usage:
   *
   * <pre>
   *   VertexOutMap out = new VertexOutMap(graph);
   *   for (int edgeId : out.edgeIds(vertexId)) {
   *     ...
   *   }
   * </pre>
   */
  public static class VertexOutMap {
    // A shared reference to the S2BuilderGraph's EdgeList.
    private EdgeList edges;
    // Maps each vertex id to the first edge id that starts at that vertex.
    private final IntVector edgeBegins = new IntVector();

    /** Constructs a new VertexOutMap and calls init(graph). */
    public VertexOutMap(S2BuilderGraph graph) {
      init(graph);
    }

    /** Initializes this VertexOutMap with the given S2BuilderGraph. */
    public void init(S2BuilderGraph graph) {
      edges = graph.edges();
      // Initialize the edgeBegins list with an entry for each vertex, containing the first edge id
      // with a source vertex id equal to the vertex. If the vertex has no out-edges, edgeBegins
      // will refer to the first edge with a higher source vertex id.
      edgeBegins.resize(graph.numVertices() + 1);
      edgeBegins.fill(-1);
      int edgeId = 0;
      for (int vertexId = 0; vertexId <= graph.numVertices(); vertexId++) {
        // Skip over the edges with source vertex ids below the current vertex id.
        while (edgeId < graph.numEdges() && graph.edges().getSrcId(edgeId) < vertexId) {
          edgeId++;
        }
        edgeBegins.set(vertexId, edgeId);
      }
    }

    /** Returns the number of edges that originate at the given vertex id. */
    public int outDegree(int vertexId) {
      int beginEdgeId = edgeBegins.get(vertexId);
      int endEdgeId = edgeBegins.get(vertexId + 1);
      return endEdgeId - beginEdgeId;
    }

    /** Returns the outgoing edges from the vertex with the given vertexId. */
    public VertexOutEdges edges(int vertexId) {
      return new VertexOutEdges(edges, edgeBegins.get(vertexId), edgeBegins.get(vertexId + 1));
    }

    /** Returns the outgoing edge ids from the vertex with the given vertexId. */
    public VertexOutEdgeIds edgeIds(int vertexId) {
      return new VertexOutEdgeIds(edgeBegins.get(vertexId), edgeBegins.get(vertexId + 1));
    }

    /** Returns the edge ids between the specified source and destination pair of vertices. */
    public VertexOutEdgeIds edgeIds(int srcVertexId, int dstVertexId) {
      // Start with [begin - end) as the range of all outgoing edges from srcVertexId.
      // Find the sub-range of edges (if any) that end at dstVertexId.
      int begin = edgeBegins.get(srcVertexId);
      int end = edgeBegins.get(srcVertexId + 1);

      // Advance 'begin' until the first edge ending at dstVertexId is found, or 'end' is reached.
      while (begin < end && edges.getDstId(begin) < dstVertexId) {
        begin++;
      }

      // If there is no edge from srcVertexId to dstVertexId, return an empty VertexOutEdgeIds.
      if (begin == end || edges.getDstId(begin) > dstVertexId) {
        return new VertexOutEdgeIds(begin, begin);
      }

      // Otherwise, advance 'last' until the last edge to dstVertexId is passed.
      int last = begin;
      while (last < end && edges.getDstId(last) == dstVertexId) {
        last++;
      }

      return new VertexOutEdgeIds(begin, last);
    }
  }

  /** A helper class for VertexInMap that represents the incoming edge *ids* to a given vertex. */
  public static class VertexInEdgeIds implements IntSequence {
    // The S2BuilderGraph's inEdgeIds, sorted in lexicographic order by destination vertex, then
    // source vertex, so that all the incoming edges to a given vertex are in consecutive order.
    private final IntArrayList inEdgeIds;
    // The index in inEdgeIds where edges to the given vertex begin.
    private final int beginIndex;
    // The index in inEdgeIds where edges to the given vertex end (exclusive).
    private final int endIndex;

    private VertexInEdgeIds(IntArrayList inEdgeIds, int begin, int end) {
      this.inEdgeIds = inEdgeIds;
      beginIndex = begin;
      endIndex = end;
    }

    @Override
    public int size() {
      return endIndex - beginIndex;
    }

    @Override
    public OfInt intIterator() {
      return new OfInt() {
        int position = beginIndex;

        @Override
        public boolean hasNext() {
          return position < endIndex;
        }

        @Override
        public int nextInt() {
          return inEdgeIds.getInt(position++);
        }
      };
    }

    @Override
    public void forEach(IntConsumer action) {
      for (int i = beginIndex; i < endIndex; i++) {
        action.accept(inEdgeIds.getInt(i));
      }
    }

    @Override
    public void forEach(IntBiConsumer action) {
      for (int i = beginIndex; i < endIndex; i++) {
        action.accept(i - beginIndex, inEdgeIds.getInt(i));
      }
    }
  }

  /**
   * A class that maps vertices to their incoming edge ids. Example usage:
   *
   * <pre>
   *   VertexInMap in = new VertexInMap(g);
   *   for (int edgeId : in.edgeIds(v)) {
   *     ...
   *   }
   * </pre>
   */
  public static class VertexInMap {
    // The S2BuilderGraph's inEdgeIds, which are edge ids sorted in lexicographic order by
    // destination vertex, then source vertex. The incoming edges to each vertex form a contiguous
    // subrange, and these subranges are in order by vertex id.
    private final IntArrayList inEdgeIds;

    // Constructed by this VertexInMap, maps each vertex id to the index in inEdgeIds where the
    // incoming edges to that vertex start. The incoming edges to this vertex end immediately before
    // the subrange for the next vertex. An entry is provided for a "next vertex" after the last
    // vertex.
    private final IntVector inEdgeBegins = new IntVector();

    /** Constructs a new VertexInMap for the given S2BuilderGraph. */
    public VertexInMap(S2BuilderGraph graph) {
      inEdgeIds = graph.getInEdgeIds();
      inEdgeBegins.resize(graph.numVertices() + 1);
      int edgeId = 0;
      for (int vertexId = 0; vertexId <= graph.numVertices(); vertexId++) {
        // Skip over the edges with destination vertex ids below the current vertex.
        while (edgeId < graph.numEdges() && graph.edgeDstId(inEdgeIds.getInt(edgeId)) < vertexId) {
          edgeId++;
        }
        inEdgeBegins.set(vertexId, edgeId);
      }
    }

    /** Returns the number of incoming edges to the given vertexId. */
    public int inDegree(int vertexId) {
      return edgeIds(vertexId).size();
    }

    /** Returns the edge ids of incoming edges to the given vertexId. */
    public VertexInEdgeIds edgeIds(int vertexId) {
      // Return a VertexInEdgeIds, which provides an iterator over inEdgeIds, from inEdgeBegins for
      // the given vertexId to inEdgeBegins for the next vertex (exclusive).
      return new VertexInEdgeIds(
          inEdgeIds, inEdgeBegins.get(vertexId), inEdgeBegins.get(vertexId + 1));
    }

    // Returns a sorted list of all incoming edges. See {@link S2BuilderGraph#getInEdgeIds()}.
    public IntArrayList inEdgeIds() {
      return inEdgeIds;
    }
  }

  /**
   * Returns an IntComparator of edge ids that orders edge ids by the smallest input edge ids mapped
   * to them, breaking ties with the edge ids themselves.
   */
  private static IntComparator byMinInputEdgeId(IntArrayList minInputEdgeIdsPerEdgeId) {
    return (int leftEdgeId, int rightEdgeId) -> {
      int minRight = minInputEdgeIdsPerEdgeId.getInt(rightEdgeId);
      int minLeft = minInputEdgeIdsPerEdgeId.getInt(leftEdgeId);
      if (minLeft < minRight) {
        return -1;
      }
      if (minLeft > minRight) {
        return 1;
      }
      // Break ties using the edge ids.
      return Integer.compare(leftEdgeId, rightEdgeId);
    };
  }

  /**
   * Indicates whether loops should be simple cycles (no repeated vertices) or circuits (which allow
   * repeated vertices but not repeated edges). In terms of how the loops are built, this
   * corresponds to closing off a loop at the first repeated vertex vs. the first repeated edge.
   */
  public enum LoopType {
    SIMPLE,
    CIRCUIT
  }

  /** TODO(torrey): Write javadoc. This isn't really documented in C++ either. */
  public enum DegenerateBoundaries {
    DISCARD,
    KEEP
  }

  /**
   * A DirectedComponent is an ArrayList of edge loops, which are arrays of edge ids (ints).
   *
   * <p>Edge loops are formed by turning left at each vertex until a repeated edge is found. Loops
   * with shared vertices form a connected component, and are grouped together into a
   * DirectedComponent. See {@link #getDirectedComponents(DegenerateBoundaries, List, S2Error)} for
   * details.
   */
  public static class DirectedComponent extends ArrayList<int[]> {}

  /**
   * An UndirectedComponent is two ArrayLists of edge loops. Edge loops are arrays of edge ids
   * (ints).
   *
   * <p>Edge loops are formed by turning left at each vertex until either a repeated vertex (for
   * LoopType.SIMPLE) or repeated edge (for LoopType.CIRCUIT) is found.
   *
   * <p>Loops with shared vertices form a connected component, and are grouped together into an
   * UndirectedComponent. Finally, the loops in each component are divided into two "complements",
   * corresponding to the two possible interpretations of the loops as either enclosing the area on
   * one side or the other of the edges. The complements are referred to as being in 'slot 0' and
   * 'slot 1'. See {@link #getUndirectedComponents(LoopType, List, S2Error)} for details.
   */
  public static class UndirectedComponent {
    // A two-element array of ArrayList<int[]> would be neater, but Java can't deal with arrays of
    // generic types.
    private ArrayList<int[]> comp0;
    private ArrayList<int[]> comp1;

    public UndirectedComponent() {
      comp0 = new ArrayList<>();
      comp1 = new ArrayList<>();
    }

    /** Adds the given edge loop to the complement in the given slot (0 or 1) */
    public void addEdgeLoop(int slot, int[] edgeLoop) {
      if (slot == 0) {
        comp0.add(edgeLoop);
      } else {
        comp1.add(edgeLoop);
      }
    }

    /** Gets the complement in the given slot. */
    public ArrayList<int[]> getComplement(int slot) {
      return (slot == 0) ? comp0 : comp1;
    }

    /** Exchanges the two complements. */
    public void swapComplements() {
      ArrayList<int[]> tmp = comp0;
      comp0 = comp1;
      comp1 = tmp;
    }
  }

  /**
   * Indicates whether polylines should be "paths" (which don't allow duplicate vertices, except
   * possibly the first and last vertex) or "walks" (which allow duplicate vertices and edges).
   */
  public enum PolylineType {
    PATH,
    WALK,
  }

  /**
   * Convenience class to return the set of labels associated with a given graph edge. Note that due
   * to snapping, one graph edge may correspond to several different input edges and will have all
   * of their labels. This class is the preferred way to retrieve edge labels.
   *
   * <p>The reason this is a class rather than a graph method is because for undirected edges, we
   * need to fetch the labels associated with both siblings. This is because only the original edge
   * of the sibling pair has labels; the automatically generated sibling edge does not.
   */
  public static class LabelFetcher {
    private S2BuilderGraph graph;
    private EdgeType edgeType;
    private IntArrayList siblingMap;

    /** Default constructor. Init must be called before use. */
    public LabelFetcher() {}

    /** Constructor which calls {@link #init(S2BuilderGraph, EdgeType)}. */
    public LabelFetcher(S2BuilderGraph graph, EdgeType edgeType) {
      init(graph, edgeType);
    }

    /**
     * Prepares to fetch labels associated with the given edge type. For EdgeType.UNDIRECTED, labels
     * associated with both edges of the sibling pair will be returned. "edgeType" is a parameter
     * (rather than using graph.options().edgeType()) so that clients can explicitly control whether
     * labels from one or both siblings are returned.
     */
    public void init(S2BuilderGraph graph, EdgeType edgeType) {
      this.graph = graph;
      this.edgeType = edgeType;
      if (edgeType == EdgeType.UNDIRECTED) {
        siblingMap = graph.getSiblingMap();
      }
    }

    /**
     * Add the labels of all the input edge ids that were snapped to the given edge to 'outLabels'.
     */
    private void addInputEdgeLabels(int edgeId, IntVector outLabels) {
      graph
          .inputEdgeIds(edgeId)
          .forEach(inputEdgeId -> outLabels.addAll(graph.labels(inputEdgeId)));
    }

    /**
     * Replaces the provided "outLabels" with the set of labels associated with edge "edgeId" (and
     * also the labels associated with the sibling of "edgeId" if edgeType() is UNDIRECTED). Labels
     * are sorted and duplicate labels are automatically removed.
     *
     * <p>This method uses an output parameter rather than returning by value in order to avoid
     * allocating a new collection for every call to this method.
     */
    public void fetch(int edgeId, IntVector outLabels) {
      outLabels.clear();
      addInputEdgeLabels(edgeId, outLabels);

      // For an undirected graph, add the labels of input edge ids snapped to the sibling.
      if (edgeType == EdgeType.UNDIRECTED) {
        addInputEdgeLabels(siblingMap.getInt(edgeId), outLabels);
      }

      // Sort and deduplicate labels.
      if (outLabels.size() > 1) {
        outLabels.sort();
        outLabels.unique();
      }
    }
  }

  /** Implements {@link S2BuilderGraph#processEdges()}. */
  private static class EdgeProcessor {

    // Specifies how the graph should be transformed.
    private final GraphOptions options;

    // The input list of Edges to process. The index of an Edge in the EdgeList is its edge id.
    private final EdgeList edges;

    // This 'inputIds' vector is a map from edge id to InputEdgeIdSetIds. An InputEdgeIdSetId
    // is a key to the idSetLexicon below.
    private final IntArrayList inputIds;

    // This IdSetLexicon maps an InputEdgeIdSetId to a set of "input edge ids".
    private final IdSetLexicon idSetLexicon;

    // 'outEdges' is an ordering of all the edge ids in 'edges' by source vertex id, so outgoing
    // edges from each vertex are consecutive. Ties are broken by destination vertex id, so
    // duplicate edges are consecutive.
    private final IntArrayList outEdges;

    // "inEdges" is an ordering of all the edge ids in 'edges' by destination vertex id, so
    // incoming edges from each vertex are consecutive. Ties are broken by source vertex id, so
    // duplicate edges are consecutive.
    private final IntArrayList inEdges;

    // As the EdgeProcessor runs, it builds a new EdgeList in 'newEdges', and a mapping from those
    // new edge ids to InputEdgeIdSetIds in 'newInputIds'. These are keys for 'idSetLexicon'. New
    // IdSets are added to the existing idSetLexicon (if needed) for the edges.
    private final EdgeList newEdges;
    private final IntArrayList newInputIds;

    /** The EdgeProcessor constructor. */
    public EdgeProcessor(
        GraphOptions options, EdgeList edges, IntArrayList inputIds, IdSetLexicon idSetLexicon) {
      this.options = new GraphOptions(options);
      this.edges = edges;
      this.inputIds = inputIds;
      this.idSetLexicon = idSetLexicon;

      // Initialize outEdges as all the edge ids sorted in lexicographic order by source vertex, so
      // all outgoing edges from a vertex are consecutive in the ordering.
      outEdges = new IntArrayList(edges.size());
      fillConsecutive(outEdges, edges.size());
      outEdges.sort(bySrcVertex(edges));

      // Initialize inEdges as all the edge ids sorted in lexicographic order by destination vertex,
      // so all incoming edges to a vertex are consecutive in the ordering.
      inEdges = new IntArrayList(edges.size());
      fillConsecutive(inEdges, edges.size());
      inEdges.sort(byDstVertex(edges));

      // Initialize the EdgeList that will contain the processed output edges, along with the new
      // vector mapping processed edge ids to input edge id sets.
      newEdges = new EdgeList();
      newEdges.ensureCapacity(edges.size());
      newInputIds = new IntArrayList();
      newInputIds.ensureCapacity(edges.size());
    }

    /**
     * Performs the graph transformations. Walks through the 'inEdges' and 'outEdges' sorted arrays
     * set up by the constructor, and performs a merge join: For each edge, gathers all the
     * duplicate copies of the edge in both directions (outgoing and incoming). Then decides what to
     * do based on 'options' and how many copies of the edge there are in each direction.
     *
     * <p>When complete, the content of the 'edges' EdgeList and 'inputIds' IntArrayList are
     * replaced with the transformed graph, and additional InputEdgeIdSets with merged sets of input
     * edge ids will have been added to the provided 'idSetLexicon' if necessary.
     */
    public void run(S2Error error) {
      int numEdges = edges.size();
      if (numEdges == 0) {
        return;
      }

      // The current position in the outEdges list.
      int out = 0;
      // The current position in the inEdges list.
      int in = 0;

      // 'outEdge' is the current outgoing edge being considered.
      final Edge outEdge = new Edge();
      edges.get(outEdges.getInt(out), outEdge);

      // 'reverseInEdge' is the current incoming edge being considered, but reversed.
      final Edge reverseInEdge = new Edge();
      edges.getReverse(inEdges.getInt(in), reverseInEdge);

      final Edge minEdge = new Edge();
      for (; ; ) {
        // We're advancing through two lists. Determine which one is behind.
        if (EdgeList.compareEdges(outEdge, reverseInEdge) < 0) {
          minEdge.set(outEdge);
        } else {
          minEdge.set(reverseInEdge);
        }

        // If the minimum edge has reached the end we are done.
        if (minEdge.isSentinel()) {
          break;
        }

        // Track where we are starting from.
        int outBegin = out;
        int inBegin = in;

        // If outEdge is the lowest, advance it until it is not. This skips duplicate outEdges.
        while (outEdge.isEqualTo(minEdge)) {
          if (++out == numEdges) {
            outEdge.setSentinel();
          } else {
            edges.get(outEdges.getInt(out), outEdge);
          }
        }
        // If reverse(inEdge) is the lowest, advance it until it is not. This skips duplicate
        // inEdges.
        while (reverseInEdge.isEqualTo(minEdge)) {
          if (++in == numEdges) {
            reverseInEdge.setSentinel();
          } else {
            edges.getReverse(inEdges.getInt(in), reverseInEdge);
          }
        }

        // How many edges were advanced?
        int nOut = out - outBegin;
        int nIn = in - inBegin;

        if (minEdge.srcId == minEdge.dstId) {
          // minEdge is degenerate.
          Preconditions.checkState(nOut == nIn);
          if (options.degenerateEdges() == DegenerateEdges.DISCARD) {
            continue;
          }
          // In the DISCARD_EXCESS case, degenerate edges that are connected to non-degenerate edges
          // should be discarded. Are there any non-degenerate incident edges?
          if (options.degenerateEdges() == DegenerateEdges.DISCARD_EXCESS
              && ((outBegin > 0 && edges.getSrcId(outEdges.getInt(outBegin - 1)) == minEdge.srcId)
                  || (out < numEdges && edges.getSrcId(outEdges.getInt(out)) == minEdge.srcId)
                  || (inBegin > 0 && edges.getDstId(inEdges.getInt(inBegin - 1)) == minEdge.srcId)
                  || (in < numEdges && edges.getDstId(inEdges.getInt(in)) == minEdge.srcId))) {
            // There were non-degenerate incident edges, so discard.
            continue;
          }

          // Are duplicate edges to be merged? DegenerateEdges.DISCARD_EXCESS also merges duplicate
          // degenerate edges, and minEdge is a degenerate edge.
          boolean merge =
              (options.duplicateEdges() == DuplicateEdges.MERGE
                  || options.degenerateEdges() == DegenerateEdges.DISCARD_EXCESS);

          if (options.edgeType() == EdgeType.UNDIRECTED
              && (options.siblingPairs() == SiblingPairs.REQUIRE
                  || options.siblingPairs() == SiblingPairs.CREATE)) {
            // When we have undirected edges and are guaranteed to have siblings, we cut the number
            // of edges in half. See the S2Builder javadoc for details.
            Preconditions.checkState((nOut & 1) == 0); // Number of edges is always even.
            addEdgeCopies(merge ? 1 : (nOut / 2), minEdge, mergeInputIds(outBegin, out));
          } else if (merge) {
            addEdgeCopies(
                options.edgeType() == EdgeType.UNDIRECTED ? 2 : 1,
                minEdge,
                mergeInputIds(outBegin, out));
          } else if (options.siblingPairs() == SiblingPairs.DISCARD
              || options.siblingPairs() == SiblingPairs.DISCARD_EXCESS) {
            // Any SiblingPair option that discards edges causes the labels of all duplicate edges
            // to be merged together. See the S2Builder javadoc for details.
            addEdgeCopies(nOut, minEdge, mergeInputIds(outBegin, out));
          } else {
            copyEdges(outBegin, out);
          }
        } else if (options.siblingPairs() == SiblingPairs.KEEP) {
          if (nOut > 1 && options.duplicateEdges() == DuplicateEdges.MERGE) {
            addEdge(minEdge, mergeInputIds(outBegin, out));
          } else {
            copyEdges(outBegin, out);
          }
        } else if (options.siblingPairs() == SiblingPairs.DISCARD) {
          if (options.edgeType() == EdgeType.DIRECTED) {
            // If nOut == nIn: balanced sibling pairs
            // If nOut < nIn:  unbalanced siblings, in the form AB, BA, BA
            // If nOut > nIn:  unbalanced siblings, in the form AB, AB, BA
            if (nOut <= nIn) {
              continue;
            }
            // Any option that discards edges causes the labels of all duplicate edges to be merged
            // together. See the S2Builder class javadoc for details.
            addEdgeCopies(
                options.duplicateEdges() == DuplicateEdges.MERGE ? 1 : (nOut - nIn),
                minEdge,
                mergeInputIds(outBegin, out));
          } else {
            if ((nOut & 1) == 0) {
              continue;
            }
            addEdge(minEdge, mergeInputIds(outBegin, out));
          }
        } else if (options.siblingPairs() == SiblingPairs.DISCARD_EXCESS) {
          if (options.edgeType() == EdgeType.DIRECTED) {
            // See comments above. The only difference is that if there are balanced sibling pairs,
            // we want to keep one such pair.
            if (nOut < nIn) {
              continue;
            }
            addEdgeCopies(
                options.duplicateEdges() == DuplicateEdges.MERGE ? 1 : max(1, nOut - nIn),
                minEdge,
                mergeInputIds(outBegin, out));
          } else {
            addEdgeCopies(((nOut & 1) == 1) ? 1 : 2, minEdge, mergeInputIds(outBegin, out));
          }
        } else {
          Preconditions.checkState(
              options.siblingPairs() == SiblingPairs.REQUIRE
                  || options.siblingPairs() == SiblingPairs.CREATE);
          if (error.ok()
              && options.siblingPairs() == SiblingPairs.REQUIRE
              && (options.edgeType() == EdgeType.DIRECTED ? (nOut != nIn) : ((nOut & 1) != 0))) {
            error.init(
                S2Error.Code.BUILDER_MISSING_EXPECTED_SIBLING_EDGES,
                "Expected all input edges to have siblings, but some were missing");
          }
          if (options.duplicateEdges() == DuplicateEdges.MERGE) {
            addEdge(minEdge, mergeInputIds(outBegin, out));
          } else if (options.edgeType() == EdgeType.UNDIRECTED) {
            // Convert graph to use directed edges instead (see documentation of REQUIRE/CREATE for
            // undirected edges).
            addEdgeCopies((nOut + 1) / 2, minEdge, mergeInputIds(outBegin, out));
          } else {
            copyEdges(outBegin, out);
            if (nIn > nOut) {
              // Automatically created edges have no input edge ids or labels.
              addEdgeCopies(nIn - nOut, minEdge, IdSetLexicon.EMPTY_SET_ID);
            }
          }
        }
      }

      // Replace the contents of the provided "edges" and their associated input id sets with the
      // processed results.
      edges.swap(newEdges);
      // TODO(torrey): Consider implementing trim on EdgeList.
      // edges.trim();

      inputIds.clear();
      inputIds.addAll(newInputIds);
      inputIds.trim();
      newInputIds.clear();
    }

    /**
     * Adds the given pair of vertex ids as a new edge, with a corresponding set of input edge ids
     * identified by the given IdSet id.
     */
    private void addEdge(Edge newEdge, int inputEdgeIdSetId) {
      newEdges.add(newEdge.srcId, newEdge.dstId);
      newInputIds.add(inputEdgeIdSetId);
    }

    /**
     * Adds 'numEdges' copies of the new edge defined by the given pair of vertex ids (with its
     * input edge ids defined by the given inputEdgeIdSetId) as new edges.
     */
    private void addEdgeCopies(int numEdges, Edge newEdge, int inputEdgeIdSetId) {
      for (int i = 0; i < numEdges; ++i) {
        newEdges.add(newEdge.srcId, newEdge.dstId);
        newInputIds.add(inputEdgeIdSetId);
      }
    }

    /**
     * Copy outgoing edges between outBegin (inclusive) and outEnd (exclusive) to newEdges, also
     * copying the corresponding input edge ids.
     */
    private void copyEdges(int outBegin, int outEnd) {
      for (int i = outBegin; i < outEnd; ++i) {
        int outEdge = outEdges.getInt(i);
        newEdges.add(edges.getSrcId(outEdge), edges.getDstId(outEdge));
        newInputIds.add(inputIds.getInt(outEdge));
      }
    }

    /**
     * Merges all the input edge id sets from the lexicon from outBegin to outEnd into a new idSet
     * and adds it to the lexicon, returning the id of that idSet.
     */
    private int mergeInputIds(int outBegin, int outEnd) {
      if (outEnd - outBegin == 1) {
        return inputIds.getInt(outEdges.getInt(outBegin)); // A single id, so no merging is needed.
      }

      IntArrayList tmpIds = new IntArrayList();
      for (int i = outBegin; i < outEnd; ++i) {
        idSetLexicon.idSet(inputIds.getInt(outEdges.getInt(i))).forEach(id -> tmpIds.add(id));
      }
      return idSetLexicon.add(tmpIds);
    }
  }
}
