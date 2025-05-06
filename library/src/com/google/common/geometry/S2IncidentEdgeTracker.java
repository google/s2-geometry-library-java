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

import com.google.common.base.Preconditions;
import com.google.common.geometry.S2Shape.MutableEdge;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.TreeMap;

/**
 * A class for detecting and tracking shape edges that are incident on the same vertex. Edges of
 * multiple shapes may be tracked, but lookup is by shape id and vertex: there is no facility to get
 * all edges of all shapes at a vertex. Edge vertices must compare exactly equal to be considered
 * the same vertex, no tolerance is applied as this isn't intended for e.g. snapping shapes
 * together, which S2Builder does better and more robustly.
 *
 * <p>To use, instantiate and then add edges with one or more sequences of calls, where each
 * sequence begins with startShape(), followed by addEdge() calls to add edges for that shape, and
 * ends with finishShape(). Those sequences do not need to visit shapes or edges in order. Then,
 * call incidentEdges() to get the resulting map from IncidentEdgeKeys (which are shapeId, vertex
 * pairs) to a set of edgeIds of the shape that are incident to that vertex.
 *
 * <p>This works on a block of edges at a time, meaning that to detect incident edges on a
 * particular vertex, we must have at least three edges incident at that vertex when finishShape()
 * is called. We don't maintain partial information between calls. However, subject to this
 * constraint, a single shape's edges may be defined with multiple sequences of startShape(),
 * addEdge()... , finishShape() calls.
 *
 * <p>The reason for this is simple: most edges don't have more than two incident edges (the
 * incoming and outgoing edge). If we had to maintain information between calls, we'd end up with a
 * map that contains every vertex, to no benefit. Instead, at each step, we discard vertices that
 * contain two or fewer incident edges.
 *
 * <p>In principle this isn't a real limitation because generally we process an S2ShapeIndex cell at
 * a time, and, if a vertex has multiple edges, we'll see all the edges in the same cell as the
 * vertex, and, in general, it's possible to aggregate edges before calling.
 *
 * <p>The tracker maintains incident edges until it's cleared. If you call it with each cell from an
 * S2ShapeIndex, then at the end you will have all the incident edge information for the whole
 * index. If only a subset is needed, call reset() when you're done.
 */
@SuppressWarnings("Assertion")
public final class S2IncidentEdgeTracker {
  // TODO(torrey): Change the data structures below to something more efficient.

  /** A tuple of (shape id, vertex) that compares by shape id. */
  public static class IncidentEdgeKey implements Comparable<IncidentEdgeKey> {
    /**
     * We need a strict ordering to be a valid key for an ordered container, but we don't actually
     * care about the ordering of the vertices (as long as they're grouped by shape id). Vertices
     * are 3D points so they don't have a natural ordering, so we'll just compare them
     * lexicographically.
     */
    static final Comparator<IncidentEdgeKey> COMPARATOR =
        Comparator.comparingInt(IncidentEdgeKey::shapeId).thenComparing(IncidentEdgeKey::vertex);

    int shapeId;
    S2Point vertex;

    /** Default constructor does not define shapeId or vertex. */
    public IncidentEdgeKey() {}

    /** Constructor that leaves the vertex undefined. */
    public IncidentEdgeKey(int shapeId) {
      this.shapeId = shapeId;
    }

    /** Constructor that sets the shapeId and vertex. */
    public IncidentEdgeKey(int shapeId, S2Point vertex) {
      this.shapeId = shapeId;
      this.vertex = vertex;
    }

    /** Returns the shape id. */
    public int shapeId() {
      return shapeId;
    }

    /** Returns the vertex. */
    public S2Point vertex() {
      return vertex;
    }

    /** Returns true if this IncidentEdgeKey is equal to the other IncidentEdgeKey. */
    public boolean isEqualTo(IncidentEdgeKey other) {
      return shapeId == other.shapeId && vertex.equalsPoint(other.vertex);
    }

    @Override
    public int compareTo(IncidentEdgeKey other) {
      return COMPARATOR.compare(this, other);
    }
  }

  /**
   * An ordered map from keys of (shape id, vertex) to a set of edgeIds. We can and do encounter the
   * same edges multiple times, so we need to use a set to deduplicate edges as they're inserted.
   */
  public static class IncidentEdgeMap extends TreeMap<IncidentEdgeKey, Set<Integer>> {
    // TODO(torrey): This implementation is convenient but very inefficient. It appears that there
    // is no need to keep the map sorted as updates are done, so instead we could use an unordered
    // map or some other structure and sort it before returning it. However, need to look at how it
    // is used by clients: do they call incidentEdges() multiple times and add more shapes and edges
    // between calls?
  }

  /** Container for temporary storage of incident edges. */
  private static class VertexEdge {
    final S2Point vertex;
    final int edgeId;

    /** Constructs a VertexEdge with the given values. */
    public VertexEdge(S2Point vertex, int edgeId) {
      this.vertex = vertex;
      this.edgeId = edgeId;
    }

    /** Factory constructor. */
    public static VertexEdge of(S2Point vertex, int edgeId) {
      return new VertexEdge(vertex, edgeId);
    }
  }

  /**
   * Scratch space to temporarily store incident edges in. Having a separate space for these lets us
   * remove vertices with only one incident edge easily.
   */
  private final List<VertexEdge> nursery = new ArrayList<>();

  /** The current shape ID that edges are being added for. */
  private int currentShapeId = -1;

  /** Actual incident edge map. */
  private final IncidentEdgeMap incidentEdgeMap = new IncidentEdgeMap();

  /** Creates a new incident edge tracker. */
  public S2IncidentEdgeTracker() {}

  /** Starts a new shape. Must be called before a sequence of addEdge calls. */
  public void startShape(int shapeId) {
    nursery.clear();
    currentShapeId = shapeId;
  }

  /** Adds an edge for the current shape to the edge tracker. */
  public void addEdge(int edgeId, MutableEdge edge) {
    addEdge(edgeId, edge.a, edge.b);
  }

  /** Adds an edge for the current shape to the edge tracker. */
  public void addEdge(int edgeId, S2Point a, S2Point b) {
    Preconditions.checkState(currentShapeId >= 0);

    // Push non-degenerate edges into the nursery.
    nursery.add(VertexEdge.of(a, edgeId));
    if (!a.equalsPoint(b)) {
      nursery.add(VertexEdge.of(b, edgeId));
    }
  }

  /**
   * finishShape() is called after all edges of a shape have been added. After calling, any vertices
   * with more than 2 incident edges will appear in the incident edge map.
   */
  public void finishShape() {
    Preconditions.checkState(currentShapeId >= 0);

    // We want to keep any vertices with more than two incident edges. We could sort the array by
    // vertex and remove any with fewer, but that would require shifting the array and could turn
    // quadratic quickly.
    //
    // Instead we'll scan forward from each vertex, swapping entries with the same vertex into a
    // contiguous range. Once we've done all the swapping we can just make sure that we have at
    // least three edges in the range.
    int nurserySize = nursery.size();
    for (int start = 0; start < nurserySize; ) {
      int end = start + 1;

      // Scan to the end of the array, swap entries so that entries with the same vertex as the
      // start are adjacent.
      int next = start;
      S2Point currVertex = nursery.get(start).vertex;
      VertexEdge tmp;
      while (++next < nurserySize) {
        if (nursery.get(next).vertex.equalsPoint(currVertex)) {
          tmp = nursery.get(next);
          nursery.set(next, nursery.get(end));
          nursery.set(end++, tmp);
        }
      }

      // Most vertices will have two incident edges (the incoming edge and the outgoing edge), which
      // aren't interesting, skip them.
      int numEdges = end - start;
      if (numEdges <= 2) {
        start = end;
        continue;
      }

      IncidentEdgeKey key = new IncidentEdgeKey(currentShapeId, nursery.get(start).vertex);
      Set<Integer> edges = incidentEdgeMap.get(key);
      // If we don't have this key yet, create it manually with a pre-sized hash set to avoid
      // rehashes. Start with a size of 8, which is four edges with a load factor of 50%.
      // TODO(user): When AndroidJdkLibsChecker is fixed, use computeIfAbsent() here. Until
      // then it must be avoided to avoid breaking Android builds.
      if (edges == null) {
        edges = new HashSet<>(8, 0.5f);
        incidentEdgeMap.put(key, edges);
      }

      while (start != end) {
        edges.add(nursery.get(start++).edgeId);
      }
    }

    // We're done with this shape, clear the nursery.
    nursery.clear();
  }

  /** Clear any accumulated state. */
  public void reset() {
    incidentEdgeMap.clear(); // Note this is linear in the table size, not constant time.
  }

  /** Returns the incident edge map. */
  public IncidentEdgeMap incidentEdges() {
    return incidentEdgeMap;
  }
}
