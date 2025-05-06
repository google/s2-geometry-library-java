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

import com.google.common.base.Preconditions;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2BuilderGraph.EdgeList;
import com.google.common.geometry.S2BuilderUtil.GraphShape;
import com.google.common.geometry.primitives.IntVector;
import com.google.common.geometry.primitives.Ints.IntList;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.List;

/**
 * S2PolygonDegeneracyFinder provides public static methods {@link
 * #findPolygonDegeneracies(S2BuilderGraph)} to find polygon degeneracies in a graph, and {@link
 * #isFullyDegenerate(S2BuilderGraph)} to determine if a graph consists entirely of degenerate edges
 * and/or sibling pairs.
 */
final class S2PolygonDegeneracyFinder {
  private final S2BuilderGraph g;
  private final S2BuilderGraph.VertexInMap in;
  private final S2BuilderGraph.VertexOutMap out;

  /** Has the vertex (by id) been visited? */
  private boolean[] isVertexUsed;

  /** Does the edge (by id) belong to a degeneracy? */
  private final boolean[] isEdgeDegeneracy;

  /** Does the vertex (by id) have unbalanced sibling pairs? */
  private final boolean[] isVertexUnbalanced;

  /** Allocated once here for reuse in computeDependencies() and getReverse(). */
  private final S2BuilderGraph.Edge reverseInEdge = new S2BuilderGraph.Edge();

  private final S2BuilderGraph.Edge outEdge = new S2BuilderGraph.Edge();

  // TODO(torrey): Is the requirement about degeneracies not coinciding with non-degenerate
  // portions of the polygon boundary described below checked? Is there a test for that check?

  /**
   * Given a graph representing a polygon, finds all degenerate edges and sibling pairs and
   * classifies them as being either shells or holes. The resulting list is sorted by edge id.
   *
   * <p>Degeneracies are not allowed to coincide with any non-degenerate portions of the polygon's
   * boundary (since that would make it impossible to classify the degeneracy as a shell or hole).
   * Specifically, degenerate edges must coincide only with other degenerate edges, and sibling
   * pairs must coincide only with other sibling pairs. (Below we require a slightly stronger
   * condition, namely that sibling pairs cannot coincide with any other edges.)
   *
   * <pre>
   * REQUIRES: g.options().edgeType() == DIRECTED
   * REQUIRES: g.options().siblingPairs() == DISCARD_EXCESS (or DISCARD)
   * REQUIRES: g.options().degenerateEdges() == DISCARD_EXCESS (or DISCARD)
   * </pre>
   *
   * <p>Usually callers will want to specify SiblingPairs.DISCARD_EXCESS and
   * DegenerateEdges.DISCARD_EXCESS in order to remove all redundant degeneracies. DISCARD is also
   * allowed in case you want to keep only one type of degeneracy (i.e., degenerate edges or sibling
   * pairs).
   *
   * <p>If the graph edges cannot be assembled into loops, the result is undefined.
   */
  public static PolygonDegeneracyList findPolygonDegeneracies(S2BuilderGraph g) {
    checkGraphOptions(g);
    if (g.options().degenerateEdges() == DegenerateEdges.DISCARD
        && g.options().siblingPairs() == SiblingPairs.DISCARD) {
      return new PolygonDegeneracyList(); // All degeneracies have already been discarded.
    }
    S2PolygonDegeneracyFinder df = new S2PolygonDegeneracyFinder(g);
    return df.run();
  }

  /**
   * Given a graph representing a polygon, returns true if the graph consists entirely of degenerate
   * edges and/or sibling pairs. Such a graph represents either the empty polygon together with a
   * collection of degenerate shells, or the full polygon together with a collection of degenerate
   * holes.
   *
   * <p><pre>
   * REQUIRES: g.options().edge_type() == DIRECTED
   * REQUIRES: g.options().sibling_pairs() == DISCARD_EXCESS (or DISCARD)
   * REQUIRES: g.options().degenerateEdges() == DISCARD_EXCESS (or DISCARD)
   */
  public static boolean isFullyDegenerate(S2BuilderGraph g) {
    checkGraphOptions(g);
    EdgeList edges = g.edges();
    for (int edgeId = 0; edgeId < g.numEdges(); ++edgeId) {
      if (edges.getSrcId(edgeId) == edges.getDstId(edgeId)) {
        continue;
      }
      if (!edges.containsSibling(edgeId)) {
        return false;
      }
    }
    return true;
  }

  /** Constructs a new S2PolygonDegeneracyFinder for the given graph. */
  private S2PolygonDegeneracyFinder(S2BuilderGraph g) {
    this.g = g;
    this.in = new S2BuilderGraph.VertexInMap(g);
    this.out = new S2BuilderGraph.VertexOutMap(g);
    isEdgeDegeneracy = new boolean[g.numEdges()];
    isVertexUnbalanced = new boolean[g.numVertices()];
  }

  /** Runs the S2PolygonDegeneracyFinder and returns any degeneracies found. */
  private PolygonDegeneracyList run() {
    // Mark all degenerate edges and sibling pairs in the "isEdgeDegeneracy" vector, and mark any
    // vertices with unbalanced edges in the "isVertexUnbalanced" vector.
    int numDegeneracies = computeDegeneracies();
    if (numDegeneracies == 0) {
      return new PolygonDegeneracyList();
    }

    // If all edges are degenerate, then use isFullPolygon() to classify the degeneracies (they are
    // necessarily all the same type).
    if (numDegeneracies == g.numEdges()) {
      return new PolygonDegeneracyList(g.numEdges(), g.isFullPolygon());
    }

    // Otherwise repeatedly build components starting from an unvisited degeneracy. (This avoids
    // building components that don't contain any degeneracies.)  Each component records the
    // "isHole" status of each degeneracy relative to the root vertex of that component. If the
    // component contains any non-degenerate portions, then we also determine whether the root
    // vertex is contained by the component (rootSign). In addition we keep track of the number of
    // components that were completely degenerate (to help us decide whether to build an index).
    ArrayList<Component> components = new ArrayList<>();
    int knownVertexId = -1;
    int knownVertexSign = 0;
    int numUnknownSigns = 0;
    isVertexUsed = new boolean[g.numVertices()];
    for (int edgeId = 0; edgeId < g.numEdges(); ++edgeId) {
      if (isEdgeDegeneracy[edgeId]) {
        int rootVertexId = g.edges().getSrcId(edgeId);
        if (isVertexUsed[rootVertexId]) {
          continue;
        }
        Component component = buildComponent(rootVertexId);
        if (component.rootSign == 0) {
          ++numUnknownSigns;
        } else {
          knownVertexId = rootVertexId;
          knownVertexSign = component.rootSign;
        }
        components.add(component);
      }
    }

    // If some components have an unknown rootSign (i.e., it is unknown whether the root vertex is
    // contained by the polygon or not), we determine the sign of those root vertices by counting
    // crossings starting from a vertex whose sign is known. Depending on how many components we
    // need to do this for, it may be worthwhile to build an index first.
    if (numUnknownSigns > 0) {
      if (knownVertexSign == 0) {
        knownVertexId = findUnbalancedVertex();
        knownVertexSign = containsVertexSign(knownVertexId);
      }
      // TODO(torrey): Tune this using a benchmark. The value of 25 is from the C++ implementation.
      int kMaxUnindexedSignComputations = 25; // Tuned using benchmarks.
      if (numUnknownSigns <= kMaxUnindexedSignComputations) {
        computeUnknownSignsBruteForce(knownVertexId, knownVertexSign, components);
      } else {
        computeUnknownSignsIndexed(knownVertexId, knownVertexSign, components);
      }
    }
    // Finally we convert the "isHole" status of each degeneracy from a relative value (compared to
    // the component's root vertex) to an absolute one, and sort all the degeneracies by edge id.
    return mergeDegeneracies(components);
  }

  private S2BuilderGraph.Edge getReverse(EdgeList edges, IntList inEdgeIds, int inIndex) {
    edges.getReverse(inEdgeIds.get(inIndex), reverseInEdge);
    return reverseInEdge;
  }

  private int computeDegeneracies() {
    int numDegeneracies = 0;

    // EdgeIds sorted in lexicographic order by (destination, origin) vertex ids. All of the
    // incoming edges to each vertex form a contiguous subrange of this ordering.
    IntList inEdgeIds = in.inEdgeIds();
    // All the edges, ordered lexicographically by (origin, destination) vertex ids. All of the
    // outgoing edges from each vertex form a contiguous subrange of this ordering.
    EdgeList edges = g.edges();
    int numEdges = g.numEdges();

    // inIndex iterates through inEdgeIds, while outEdgeId iterates through edges. Both of those
    // are numEdges long.
    for (int inIndex = 0, outEdgeId = 0; outEdgeId < numEdges; ++outEdgeId) {
      edges.get(outEdgeId, outEdge);
      if (outEdge.srcId == outEdge.dstId) {
        isEdgeDegeneracy[outEdgeId] = true;
        ++numDegeneracies;
      } else {
        // Advance inIndex until reverse(inEdge) is equal to or greater than outEdge.
        while (inIndex < numEdges && getReverse(edges, inEdgeIds, inIndex).isLessThan(outEdge)) {
          ++inIndex;
        }
        // If reverse(inEdge) equals outEdge, they are siblings.
        if (inIndex < numEdges && getReverse(edges, inEdgeIds, inIndex).equals(outEdge)) {
          isEdgeDegeneracy[outEdgeId] = true;
          ++numDegeneracies;
        } else {
          // This edge does not have a sibling, which mean that we can determine whether either
          // vertex is contained by the polygon (using semi-open boundaries) by examining only the
          // edges incident to that vertex. We only mark the first vertex since there is no
          // advantage to finding more than one unbalanced vertex per connected component.
          isVertexUnbalanced[outEdge.srcId] = true;
        }
      }
    }
    return numDegeneracies;
  }

  /**
   * For tracking the frontier while building a component. Holds a vertex ID and the containment
   * status of that vertex.
   */
  private static class VertexIdContained {
    final int vertexId;
    final boolean contained;

    public VertexIdContained(int vertexId, boolean contained) {
      this.vertexId = vertexId;
      this.contained = contained;
    }

    public static VertexIdContained of(int vertexId, boolean contained) {
      return new VertexIdContained(vertexId, contained);
    }
  }

  /**
   * Build a connected component starting at the given root vertex. The information returned
   * includes: the root vertex, whether the containment status of the root vertex could be
   * determined using only the edges in this component, and a vector of the edges that belong to
   * degeneracies along with the shell/hole status of each such edge relative to the root vertex.
   */
  private Component buildComponent(int rootVertexId) {
    Component result = new Component();
    result.rootVertexId = rootVertexId;
    // We keep track of the frontier of unexplored vertices as a stack of vertex ids, and whether
    // each vertex is on the same side of the polygon boundary as the root vertex.
    ArrayDeque<VertexIdContained> frontier = new ArrayDeque<>();
    frontier.push(VertexIdContained.of(rootVertexId, true));
    isVertexUsed[rootVertexId] = true;
    while (!frontier.isEmpty()) {
      // Get a vertex and it's containment status from the frontier stack.
      VertexIdContained vertexIdContained = frontier.pop();
      int v0 = vertexIdContained.vertexId;
      // Is v0 on the same side of the polygon boundary as the root vertex?
      boolean v0SameInside = vertexIdContained.contained;
      // If the containment status of the root vertex isn't known, and v0 has unbalanced sibling
      // pairs, compute the sign of v0. Use that to set the containment status of the root vertex.
      if (result.rootSign == 0 && isVertexUnbalanced[v0]) {
        int v0Sign = containsVertexSign(v0);
        Preconditions.checkState(v0Sign != 0); // TODO(torrey): Debug-mode only check?
        result.rootSign = v0SameInside ? v0Sign : -v0Sign;
      }
      // Keep building the connected component, considering out-edges from v0.
      OfInt outEdges = out.edgeIds(v0).intIterator();
      while (outEdges.hasNext()) {
        int edgeId = outEdges.nextInt();
        // Get the endpoint of the out-edge, and determine if the edge is a shell or a hole relative
        // to the root vertex.
        int v1 = g.edges().getDstId(edgeId);
        boolean sameInside = v0SameInside ^ crossingParity(v0, v1, false);
        if (isEdgeDegeneracy[edgeId]) {
          result.degeneracies.add(edgeId, sameInside);
        }
        if (isVertexUsed[v1]) {
          continue;
        }
        sameInside ^= crossingParity(v1, v0, true);
        frontier.push(VertexIdContained.of(v1, sameInside));
        isVertexUsed[v1] = true;
      }
    }
    return result;
  }

  /**
   * Counts the number of times that (v0, v1) crosses the edges incident to v0, and returns the
   * result modulo 2. This is equivalent to calling S2EdgeUtil.vertexCrossing for the edges incident
   * to v0, except that this implementation is more efficient (since it doesn't need to determine
   * which two edge vertices are the same.
   *
   * <p>If "includeSame" is false, then the edge (v0, v1) and its sibling (v1, v0) (if any) are
   * excluded from the parity calculation.
   */
  private boolean crossingParity(int v0, int v1, boolean includeSame) {
    int crossings = 0;
    S2Point p0 = g.vertex(v0);
    S2Point p1 = g.vertex(v1);
    S2Point p0Ref = S2.ortho(p0);

    S2BuilderGraph.VertexOutEdges outEdges = out.edges(v0);
    for (int i = 0; i < outEdges.size(); i++) {
      int dstId = outEdges.getDstId(i);
      if (dstId == v1) {
        if (includeSame) {
          ++crossings;
        }
      } else if (S2Predicates.orderedCCW(p0Ref, g.vertex(dstId), p1, p0)) {
        ++crossings;
      }
    }

    OfInt i = in.edgeIds(v0).intIterator();
    while (i.hasNext()) {
      int srcId = g.edges().getSrcId(i.nextInt());
      if (srcId == v1) {
        if (includeSame) {
          ++crossings;
        }
      } else if (S2Predicates.orderedCCW(p0Ref, g.vertex(srcId), p1, p0)) {
        ++crossings;
      }
    }
    return (crossings & 1) != 0;
  }

  /** Returns the first unbalanced vertex id. */
  private int findUnbalancedVertex() {
    for (int vertexId = 0; vertexId < g.numVertices(); ++vertexId) {
      if (isVertexUnbalanced[vertexId]) {
        return vertexId;
      }
    }
    throw new IllegalStateException("Could not find previously marked unbalanced vertex");
  }

  private int containsVertexSign(int v0) {
    S2ContainsVertexQuery query = new S2ContainsVertexQuery(g.vertex(v0));
    S2BuilderGraph.VertexOutEdges outEdges = out.edges(v0);
    for (int i = 0; i < outEdges.size(); i++) {
      int dstId = outEdges.getDstId(i);
      query.addOutgoing(g.vertex(dstId));
    }
    OfInt i = in.edgeIds(v0).intIterator();
    while (i.hasNext()) {
      int srcId = g.edges().getSrcId(i.nextInt());
      query.addIncoming(g.vertex(srcId));
    }
    return query.containsSign();
  }

  /**
   * Determines any unknown signs of component root vertices by counting crossings starting from a
   * vertex whose sign is known. This version simply tests all edges for crossings.
   */
  private void computeUnknownSignsBruteForce(
      int knownVertexId, int knownVertexSign, List<Component> components) {
    S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser();
    for (Component component : components) {
      if (component.rootSign != 0) {
        continue;
      }
      boolean inside = knownVertexSign > 0;
      crosser.init(g.vertex(knownVertexId), g.vertex(component.rootVertexId));
      for (int edgeId = 0; edgeId < g.numEdges(); ++edgeId) {
        if (isEdgeDegeneracy[edgeId]) {
          continue;
        }
        inside ^=
            crosser.edgeOrVertexCrossing(
                g.vertex(g.edges().getSrcId(edgeId)), g.vertex(g.edges().getDstId(edgeId)));
      }
      component.rootSign = inside ? 1 : -1;
    }
  }

  /**
   * Like computeUnknownSignsBruteForce, except that this method uses an index to find the set of
   * edges that cross a given edge.
   */
  private void computeUnknownSignsIndexed(
      int knownVertexId, int knownVertexSign, List<Component> components) {
    // Add all the edges in the graph to a shape index as a single shape, and build a query.
    S2ShapeIndex index = new S2ShapeIndex();
    index.add(new GraphShape(g));
    S2EdgeQuery query = new S2EdgeQuery(index);

    // ArrayList<ShapeEdgeId> crossingEdges = new ArrayList<>();
    S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser();
    for (Component component : components) {
      if (component.rootSign != 0) {
        continue;
      }
      boolean inside = knownVertexSign > 0;
      crosser.init(g.vertex(knownVertexId), g.vertex(component.rootVertexId));
      // Instead of checking every edge in the graph, get only edges that might cross the edge from
      // the known vertex to the root vertex of the component.
      S2EdgeQuery.Edges crossingEdges =
          query.getCandidates(g.vertex(knownVertexId), g.vertex(component.rootVertexId), 0);

      while (!crossingEdges.isEmpty()) {
        int edgeId = crossingEdges.nextEdge();
        if (isEdgeDegeneracy[edgeId]) {
          continue;
        }
        inside ^=
            crosser.edgeOrVertexCrossing(
                g.vertex(g.edges().getSrcId(edgeId)), g.vertex(g.edges().getDstId(edgeId)));
      }
      component.rootSign = inside ? 1 : -1;
    }
  }

  /**
   * Merges the degeneracies from all components together, and computes the final "isHole" status of
   * each edge (since up to this point, the "isHole" value has been expressed relative to the root
   * vertex of each component). Returns the degeneracies in sorted order.
   */
  private static PolygonDegeneracyList mergeDegeneracies(List<Component> components) {
    PolygonDegeneracyList result = new PolygonDegeneracyList();
    for (Component component : components) {
      // TODO(torrey): Debug-mode only check?
      Preconditions.checkState(component.rootSign != 0);
      boolean invert = component.rootSign < 0;
      for (int d = 0; d < component.degeneracies.size(); ++d) {
        result.add(component.degeneracies.edgeId(d), component.degeneracies.isHole(d) ^ invert);
      }
    }
    result.sort();
    return result;
  }

  private static void checkGraphOptions(S2BuilderGraph g) {
    Preconditions.checkState(g.options().edgeType() == EdgeType.DIRECTED);
    Preconditions.checkState(
        g.options().degenerateEdges() == DegenerateEdges.DISCARD
            || g.options().degenerateEdges() == DegenerateEdges.DISCARD_EXCESS);
    Preconditions.checkState(
        g.options().siblingPairs() == SiblingPairs.DISCARD
            || g.options().siblingPairs() == SiblingPairs.DISCARD_EXCESS);
  }

  /**
   * A list of polygon degeneracies. A polygon degeneracy is either a degenerate edge (an edge from
   * a vertex to itself) or a sibling edge pair (consisting of an edge and its corresponding reverse
   * edge). "isHole" indicates whether the degeneracy corresponds to a polygon hole (as opposed to a
   * polygon shell).
   */
  public static class PolygonDegeneracyList {
    /** Encoded degeneracies, which consist of a non-negative edgeId and isHole boolean. */
    private final IntVector encoded = IntVector.empty();

    /** Constructs an empty PolygonDegeneracyList. */
    public PolygonDegeneracyList() {}

    /**
     * Constructs a PolygonDegeneracyList with 'numEdges' degeneracies on consecutive edge ids
     * starting from zero, where each has the same, given 'isHole' property.
     */
    public PolygonDegeneracyList(int numEdges, boolean isHole) {
      for (int i = 0; i < numEdges; ++i) {
        encoded.add(encode(i, isHole));
      }
    }

    /** Add a degeneracy to the list. */
    public void add(int edgeId, boolean isHole) {
      encoded.add(encode(edgeId, isHole));
    }

    /** Returns the number of degeneracies in the list. */
    public int size() {
      return encoded.size();
    }

    /** Decodes and returns the edgeId for the degeneracy with the given index. */
    public int edgeId(int index) {
      checkIndex(index);
      return decodeEdgeId(encoded.get(index));
    }

    /** Decodes and returns the isHole property of the degeneracy with the given index. */
    public boolean isHole(int index) {
      checkIndex(index);
      return encoded.get(index) < 0;
    }

    /** Sorts the polygon degeneracies by edge id, breaking ties with "isHole". */
    public void sort() {
      encoded.sort(
          (a, b) -> {
            int edgeIdA = decodeEdgeId(a);
            int edgeIdB = decodeEdgeId(b);
            if (edgeIdA < edgeIdB) {
              return -1;
            } else if (edgeIdA > edgeIdB) {
              return 1;
            }

            // Holes (true) should come before shells (false), so the order is reversed.
            return Boolean.compare(decodeIsHole(b), decodeIsHole(a));
          });
    }

    private void checkIndex(int index) {
      if (index < 0 || index >= encoded.size()) {
        throw new IndexOutOfBoundsException(
            "index " + index + " out of bounds on list of size " + encoded.size());
      }
    }

    /** Encodes the given edgeId and isHole into an integer. */
    private int encode(int edgeId, boolean isHole) {
      if (isHole) {
        edgeId |= Integer.MIN_VALUE;
      }
      return edgeId;
    }

    private int decodeEdgeId(int encoded) {
      return encoded & Integer.MAX_VALUE;
    }

    private boolean decodeIsHole(int encoded) {
      return encoded < 0;
    }
  }

  /**
   * The algorithm builds a set of connected components containing all edges that form degeneracies.
   * The shell/hole status of each degeneracy is initially unknown, and is expressed relative to the
   * root vertex: "isHole" means that the degeneracy is a hole if and only if the root vertex turns
   * out to be inside the polygon.
   */
  private static class Component {
    /** The root vertex from which this component was built. */
    public int rootVertexId;

    /** +1 if "rootVertexId" is inside the polygon, -1 if outside, and 0 if unknown. */
    public int rootSign = 0;

    /**
     * The degeneracies found in this component. "isHole" is expressed relative to the root vertex:
     * the degeneracy is a hole iff the root vertex turns out to be inside the polygon (i.e.,
     * rootSign > 0).
     */
    public final PolygonDegeneracyList degeneracies = new PolygonDegeneracyList();
  }
}
