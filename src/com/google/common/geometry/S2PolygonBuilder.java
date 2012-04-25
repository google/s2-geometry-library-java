/*
 * Copyright 2006 Google Inc.
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
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import com.google.common.collect.TreeMultimap;

import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.Stack;

/**
 * This is a simple class for assembling polygons out of edges. It requires that
 * no two edges cross. It can handle both directed and undirected edges, and
 * optionally it can also remove duplicate edge pairs (consisting of two
 * identical edges or an edge and its reverse edge). This is useful for
 * computing seamless unions of polygons that have been cut into pieces.
 *
 *  Here are some of the situations this class was designed to handle:
 *
 *  1. Computing the union of disjoint polygons that may share part of their
 * boundaries. For example, reassembling a lake that has been split into two
 * loops by a state boundary.
 *
 *  2. Constructing polygons from input data that does not follow S2
 * conventions, i.e. where loops may have repeated vertices, or distinct loops
 * may share edges, or shells and holes have opposite or unspecified
 * orientations.
 *
 *  3. Computing the symmetric difference of a set of polygons whose edges
 * intersect only at vertices. This can be used to implement a limited form of
 * polygon intersection or subtraction as well as unions.
 *
 *  4. As a tool for implementing other polygon operations by generating a
 * collection of directed edges and then assembling them into loops.
 *
 */
public strictfp class S2PolygonBuilder {
  private Options options;

  /**
   * The current set of edges, grouped by origin. The set of destination
   * vertices is a multiset so that the same edge can be present more than once.
   */
  private final Map<S2Point, Multiset<S2Point>> edges = Maps.newHashMap();

  /**
   * Unique collection of the starting (first) vertex of all edges, in the order
   * they are added to {@code edges}.
   */
  private final List<S2Point> startingVertices = Lists.newArrayList();

  /**
   * Default constructor for well-behaved polygons. Uses the DIRECTED_XOR
   * options.
   */
  public S2PolygonBuilder() {
    this(Options.DIRECTED_XOR);

  }

  public S2PolygonBuilder(Options options) {
    this.options = options;
  }

  public static class Options {
    /**
     * These are the options that should be used for assembling well-behaved
     * input data into polygons. All edges should be directed such that "shells"
     * and "holes" have opposite orientations (typically CCW shells and
     * clockwise holes), unless it is known that shells and holes do not share
     * any edges.
     */
    public static final Options DIRECTED_XOR = new Options(false, true);

    /**
     * These are the options that should be used for assembling polygons that do
     * not follow the conventions above, e.g. where edge directions may vary
     * within a single loop, or shells and holes are not oppositely oriented.
     */
    public static final Options UNDIRECTED_XOR = new Options(true, true);

    /**
     * These are the options that should be used for assembling edges where the
     * desired output is a collection of loops rather than a polygon, and edges
     * may occur more than once. Edges are treated as undirected and are not
     * XORed together, in particular, adding edge A->B also adds B->A.
     */
    public static final Options UNDIRECTED_UNION = new Options(true, false);

    /**
     * Finally, select this option when the desired output is a collection of
     * loops rather than a polygon, but your input edges are directed and you do
     * not want reverse edges to be added implicitly as above.
     */
    public static final Options DIRECTED_UNION = new Options(false, false);

    private boolean undirectedEdges;
    private boolean xorEdges;
    private boolean validate;
    private S1Angle mergeDistance;

    /**
     * If positive, this is the fraction of {@code mergeDistance} that a point
     * can be from an edge before it will be snapped to the edge.
     */
    private double edgeSpliceFraction;

    public Options(boolean undirectedEdges, boolean xorEdges) {
      this.undirectedEdges = undirectedEdges;
      this.xorEdges = xorEdges;
      this.validate = false;
      this.mergeDistance = S1Angle.radians(0);
      this.edgeSpliceFraction = 0.866;
    }

    /**
     * If "undirected_edges" is false, then the input is assumed to consist of
     * edges that can be assembled into oriented loops without reversing any of
     * the edges. Otherwise, "undirected_edges" should be set to true.
     */
    public boolean getUndirectedEdges() {
      return undirectedEdges;
    }

    /**
     * If "xor_edges" is true, then any duplicate edge pairs are removed. This
     * is useful for computing the union of a collection of polygons whose
     * interiors are disjoint but whose boundaries may share some common edges
     * (e.g. computing the union of South Africa, Lesotho, and Swaziland).
     *
     *  Note that for directed edges, a "duplicate edge pair" consists of an
     * edge and its corresponding reverse edge. This means that either (a)
     * "shells" and "holes" must have opposite orientations, or (b) shells and
     * holes do not share edges. Otherwise undirected_edges() should be
     * specified.
     *
     *  There are only two reasons to turn off xor_edges():
     *
     *  (1) assemblePolygon() will be called, and you want to assert that there
     * are no duplicate edge pairs in the input.
     *
     *  (2) assembleLoops() will be called, and you want to keep abutting loops
     * separate in the output rather than merging their regions together (e.g.
     * assembling loops for Kansas City, KS and Kansas City, MO simultaneously).
     */
    public boolean getXorEdges() {
      return xorEdges;
    }

    /**
     * Default value: false
     */
    public boolean getValidate() {
      return validate;
    }

    /**
     * Default value: 0
     */
    public S1Angle getMergeDistance() {
      return mergeDistance;
    }

    /**
     * If true, isValid() is called on all loops and polygons before
     * constructing them. If any loop is invalid (e.g. self-intersecting), it is
     * rejected and returned as a set of "unused edges". Any remaining valid
     * loops are kept. If the entire polygon is invalid (e.g. two loops
     * intersect), then all loops are rejected and returned as unused edges.
     */
    public void setValidate(boolean validate) {
      this.validate = validate;
    }

    /**
     * <p>If set to a positive value, all vertex pairs that are separated by
     * less than this distance will be merged together. Note that vertices can
     * move arbitrarily far if there is a long chain of vertices separated by
     * less than this distance.
     *
     * <p>
     * This method is useful for assembling polygons out of input data where
     * vertices and/or edges may not be perfectly aligned.
     *
     * Default value: 0.
     */
    public void setMergeDistance(S1Angle mergeDistance) {
      this.mergeDistance = mergeDistance;
    }

    // Used for testing only
    void setUndirectedEdges(boolean undirectedEdges) {
      this.undirectedEdges = undirectedEdges;
    }

    // Used for testing only
    void setXorEdges(boolean xorEdges) {
      this.xorEdges = xorEdges;
    }

    /**
     * Returns the edge splice fraction, which defaults to 0.866 (approximately
     * sqrt(3)/2).
     *
     * <p>
     * The edge splice radius is automatically set to this fraction of the
     * vertex merge radius. If the edge splice radius is positive, then all
     * vertices that are closer than this distance to an edge are spliced into
     * that edge. Note that edges can move arbitrarily far if there is a long
     * chain of vertices in just the right places.
     *
     * <p>
     * You can turn off edge splicing by setting this value to zero. This will
     * save some time if you don't need this feature, or you don't want vertices
     * to be spliced into nearby edges for some reason.
     *
     * <p>
     * Note that the edge splice fraction must be less than sqrt(3)/2 in order
     * to avoid infinite loops in the merge algorithm. The default value is very
     * close to this maximum and therefore results in the maximum amount of edge
     * splicing for a given vertex merge radius.
     *
     * <p>
     * The only reason to reduce the edge splice fraction is if you want to
     * limit changes in edge direction due to splicing. The direction of an edge
     * can change by up to asin(edge_splice_fraction) due to each splice. Thus
     * by default, edges are allowed to change direction by up to 60 degrees per
     * splice. However, note that most direction changes are much smaller than
     * this: the worst case occurs only if the vertex being spliced is just
     * outside the vertex merge radius from one of the edge endpoints.
     */
    public double getEdgeSpliceFraction() {
      return edgeSpliceFraction;
    }

    public void setEdgeSpliceFraction(double edgeSpliceFraction) {
      Preconditions.checkState(edgeSpliceFraction < Math.sqrt(3) / 2,
          "Splice fraction must be at least sqrt(3)/2 to ensure termination " +
          "of edge splicing algorithm.");
      this.edgeSpliceFraction = edgeSpliceFraction;
    }
  }

  public Options options() {
    return options;
  }

  /**
   * Add the given edge to the polygon builder and returns true if the edge was
   * actually added to the edge graph.
   *
   * <p>This method should be used for input data that may not follow S2 polygon
   * conventions. Note that edges are not allowed to cross each other. Also note
   * that as a convenience, edges where v0 == v1 are ignored.
   */
  public boolean addEdge(S2Point v0, S2Point v1) {
    if (v0.equals(v1)) {
      return false;
    }

    // If xor_edges is true, we look for an existing edge in the opposite
    // direction. We either delete that edge or insert a new one.
    if (options.getXorEdges() && hasEdge(v1, v0)) {
      eraseEdge(v1, v0);
      return false;
    }

    if (edges.get(v0) == null) {
      edges.put(v0, HashMultiset.<S2Point>create());
      startingVertices.add(v0);
    }

    edges.get(v0).add(v1);
    if (options.getUndirectedEdges()) {
      if (edges.get(v1) == null) {
        edges.put(v1, HashMultiset.<S2Point>create());
        startingVertices.add(v1);
      }
      edges.get(v1).add(v0);
    }

    return true;
  }

  /**
   * Add all edges in the given loop. If the sign() of the loop is negative
   * (i.e. this loop represents a hole), the reverse edges are added instead.
   * This implies that "shells" are CCW and "holes" are CW, as required for the
   * directed edges convention described above.
   *
   * This method does not take ownership of the loop.
   */
  public void addLoop(S2Loop loop) {
    int sign = loop.sign();
    for (int i = loop.numVertices(); i > 0; --i) {
      // Vertex indices need to be in the range [0, 2*num_vertices()-1].
      addEdge(loop.vertex(i), loop.vertex(i + sign));
    }
  }

  /**
   * Add all loops in the given polygon. Shells and holes are added with
   * opposite orientations as described for AddLoop(). This method does not take
   * ownership of the polygon.
   */
  public void addPolygon(S2Polygon polygon) {
    for (int i = 0; i < polygon.numLoops(); ++i) {
      addLoop(polygon.loop(i));
    }
  }

  /**
   * Assembles the given edges into as many non-crossing loops as possible. When
   * there is a choice about how to assemble the loops, then CCW loops are
   * preferred. Returns true if all edges were assembled. If "unused_edges" is
   * not NULL, it is initialized to the set of edges that could not be assembled
   * into loops.
   *
   *  Note that if xor_edges() is false and duplicate edge pairs may be present,
   * then undirected_edges() should be specified unless all loops can be
   * assembled in a counter-clockwise direction. Otherwise this method may not
   * be able to assemble all loops due to its preference for CCW loops.
   *
   * This method resets the S2PolygonBuilder state so that it can be reused.
   */
  public boolean assembleLoops(List<S2Loop> loops, List<S2Edge> unusedEdges) {
    if (options.getMergeDistance().radians() > 0) {
      PointIndex index = new PointIndex(
          options.getMergeDistance().radians(),
          options.getEdgeSpliceFraction());
      Map<S2Point, S2Point> mergeMap = buildMergeMap(index);
      moveVertices(mergeMap);
      if (options.getEdgeSpliceFraction() > 0) {
        spliceEdges(index);
      }
    }

    List<S2Edge> dummyUnusedEdges = Lists.newArrayList();
    if (unusedEdges == null) {
      unusedEdges = dummyUnusedEdges;
    }

    // We repeatedly choose an edge and attempt to assemble a loop
    // starting from that edge.  (This is always possible unless the
    // input includes extra edges that are not part of any loop.)  To
    // maintain a consistent scanning order over edges between
    // different machine architectures (e.g. 'clovertown' vs. 'opteron'),
    // we follow the order they were added to the builder.
    unusedEdges.clear();
    for (int i = 0; i < startingVertices.size(); ) {
      S2Point v0 = startingVertices.get(i);
      Multiset<S2Point> candidates = edges.get(v0);
      if (candidates == null) {
        i++;
        continue;
      }
      S2Point v1 = candidates.iterator().next();
      S2Loop loop = assembleLoop(v0, v1, unusedEdges);
      if (loop != null) {
        eraseLoop(loop, loop.numVertices());
        loops.add(loop);
      }
    }
    return unusedEdges.isEmpty();
  }

  /**
   * Like AssembleLoops, but normalizes all the loops so that they enclose less
   * than half the sphere, and then assembles the loops into a polygon.
   *
   *  For this method to succeed, there should be no duplicate edges in the
   * input. If this is not known to be true, then the "xor_edges" option should
   * be set (which is true by default).
   *
   *  Note that S2Polygons cannot represent arbitrary regions on the sphere,
   * because of the limitation that no loop encloses more than half of the
   * sphere. For example, an S2Polygon cannot represent a 100km wide band around
   * the equator. In such cases, this method will return the *complement* of the
   * expected region. So for example if all the world's coastlines were
   * assembled, the output S2Polygon would represent the land area (irrespective
   * of the input edge or loop orientations).
   */
  public boolean assemblePolygon(S2Polygon polygon, List<S2Edge> unusedEdges) {
    List<S2Loop> loops = Lists.newArrayList();
    boolean success = assembleLoops(loops, unusedEdges);

    // If edges are undirected, then all loops are already CCW. Otherwise we
    // need to make sure the loops are normalized.
    if (!options.getUndirectedEdges()) {
      for (S2Loop loop : loops) {
        loop.normalize();
      }
    }
    if (options.getValidate() && !S2Polygon.isValid(loops)) {
      if (unusedEdges != null) {
        for (S2Loop loop : loops) {
          rejectLoop(loop, loop.numVertices(), unusedEdges);
        }
      }
      return false;
    }
    polygon.init(loops);
    return success;
  }

  /**
   * Convenience method for when you don't care about unused edges.
   */
  public S2Polygon assemblePolygon() {
    S2Polygon polygon = new S2Polygon();
    List<S2Edge> unusedEdges = Lists.newArrayList();

    assemblePolygon(polygon, unusedEdges);

    return polygon;
  }

  private void eraseEdge(S2Point v0, S2Point v1) {
    // Note that there may be more than one copy of an edge if we are not XORing
    // them, so a VertexSet is a multiset.

    Multiset<S2Point> vset = edges.get(v0);
    // assert (vset.count(v1) > 0);
    vset.remove(v1);
    if (vset.isEmpty()) {
      edges.remove(v0);
    }

    if (options.getUndirectedEdges()) {
      vset = edges.get(v1);
      // assert (vset.count(v0) > 0);
      vset.remove(v0);
      if (vset.isEmpty()) {
        edges.remove(v1);
      }
    }
  }

  private void eraseLoop(List<S2Point> v, int n) {
    for (int i = n - 1, j = 0; j < n; i = j++) {
      eraseEdge(v.get(i), v.get(j));
    }
  }

  private void eraseLoop(S2Loop v, int n) {
    for (int i = n - 1, j = 0; j < n; i = j++) {
      eraseEdge(v.vertex(i), v.vertex(j));
    }
  }

  /**
   * We start at the given edge and assemble a loop taking left turns whenever
   * possible. We stop the loop as soon as we encounter any vertex that we have
   * seen before *except* for the first vertex (v0). This ensures that only CCW
   * loops are constructed when possible.
   */
  private S2Loop assembleLoop(S2Point v0, S2Point v1, List<S2Edge> unusedEdges) {

    // The path so far.
    List<S2Point> path = Lists.newArrayList();

    // Maps a vertex to its index in "path".
    Map<S2Point, Integer> index = Maps.newHashMap();
    path.add(v0);
    path.add(v1);

    index.put(v1, 1);

    while (path.size() >= 2) {
      // Note that "v0" and "v1" become invalid if "path" is modified.
      v0 = path.get(path.size() - 2);
      v1 = path.get(path.size() - 1);

      S2Point v2 = null;
      boolean v2Found = false;
      Multiset<S2Point> vset = edges.get(v1);
      if (vset != null) {
        for (S2Point v : vset) {
          // We prefer the leftmost outgoing edge, ignoring any reverse edges.
          if (v.equals(v0)) {
            continue;
          }
          if (!v2Found || S2.orderedCCW(v0, v2, v, v1)) {
            v2 = v;
          }
          v2Found = true;
        }
      }
      if (!v2Found) {
        // We've hit a dead end. Remove this edge and backtrack.
        unusedEdges.add(new S2Edge(v0, v1));
        eraseEdge(v0, v1);
        index.remove(v1);
        path.remove(path.size() - 1);
      } else if (index.get(v2) == null) {
        // This is the first time we've visited this vertex.
        index.put(v2, path.size());
        path.add(v2);
      } else {
        // We've completed a loop. In a simple case last edge is the same as the
        // first one (since we did not add the very first vertex to the index we
        // would not know that the loop is completed till the second vertex is
        // examined). In this case we just remove the last edge to preserve the
        // original vertex order. In a more complicated case the edge that
        // closed the loop is different and we should remove initial vertices
        // that are not part of the loop.
        if (index.get(v2) == 1 && path.get(0).equals(path.get(path.size() - 1))) {
          path.remove(path.size() - 1);
        } else {
          // We've completed a loop. Throw away any initial vertices that
          // are not part of the loop.
          path = path.subList(index.get(v2), path.size());
        }

        S2Loop loop = new S2Loop(path);
        if (options.getValidate() && !loop.isValid()) {
          // We've constructed a loop that crosses itself, which can only happen
          // if there is bad input data. Throw away the whole loop.
          rejectLoop(path, path.size(), unusedEdges);
          eraseLoop(path, path.size());
          return null;
        }

        // In the case of undirected edges, we may have assembled a clockwise
        // loop while trying to assemble a CCW loop.  To fix this, we assemble
        // a new loop starting with an arbitrary edge in the reverse direction.
        // This is guaranteed to assemble a loop that is interior to the previous
        // one and will therefore eventually terminate.
        if (options.getUndirectedEdges() && !loop.isNormalized()) {
          return assembleLoop(path.get(1), path.get(0), unusedEdges);
        }

        return loop;
      }
    }
    return null;
  }

  /** Erases all edges of the given loop and marks them as unused. */
  private void rejectLoop(S2Loop v, int n, List<S2Edge> unusedEdges) {
    for (int i = n - 1, j = 0; j < n; i = j++) {
      unusedEdges.add(new S2Edge(v.vertex(i), v.vertex(j)));
    }
  }

  /** Erases all edges of the given loop and marks them as unused. */
  private void rejectLoop(List<S2Point> v, int n, List<S2Edge> unusedEdges) {
    for (int i = n - 1, j = 0; j < n; i = j++) {
      unusedEdges.add(new S2Edge(v.get(i), v.get(j)));
    }
  }

  /** Moves a set of vertices from old to new positions. */
  private void moveVertices(Map<S2Point, S2Point> mergeMap) {
    if (mergeMap.isEmpty()) {
      return;
    }

    // We need to copy the set of edges affected by the move, since
    // this.edges_could be reallocated when we start modifying it.
    List<S2Edge> edgesCopy = Lists.newArrayList();
    for (Map.Entry<S2Point, Multiset<S2Point>> edge : this.edges.entrySet()) {
      S2Point v0 = edge.getKey();
      Multiset<S2Point> vset = edge.getValue();
      for (S2Point v1 : vset) {
        if (mergeMap.get(v0) != null || mergeMap.get(v1) != null) {

          // We only need to modify one copy of each undirected edge.
          if (!options.getUndirectedEdges() || v0.lessThan(v1)) {
            edgesCopy.add(new S2Edge(v0, v1));
          }
        }
      }
    }

    // Now erase all the old edges, and add all the new edges. This will
    // automatically take care of any XORing that needs to be done, because
    // EraseEdge also erases the sibiling of undirected edges.
    for (int i = 0; i < edgesCopy.size(); ++i) {
      S2Point v0 = edgesCopy.get(i).getStart();
      S2Point v1 = edgesCopy.get(i).getEnd();
      eraseEdge(v0, v1);
      S2Point new0 = mergeMap.get(v0);
      if (new0 != null) {
        v0 = new0;
      }
      S2Point new1 = mergeMap.get(v1);
      if (new1 != null) {
        v1 = new1;
      }
      addEdge(v0, v1);
    }
  }

  /**
   * Clusters vertices that are separated by at most merge_distance() and
   * returns a map of each one to a single representative vertex for all the
   * vertices in the cluster.
   */
  private Map<S2Point, S2Point> buildMergeMap(PointIndex index) {
    // The overall strategy is to start from each vertex and grow a maximal
    // cluster of mergeable vertices. In graph theoretic terms, we find the
    // connected components of the undirected graph whose edges connect pairs of
    // vertices that are separated by at most vertex_merge_radius().
    //
    // We then choose a single representative vertex for each cluster, and
    // update "merge_map" appropriately. We choose an arbitrary existing
    // vertex rather than computing the centroid of all the vertices to avoid
    // creating new vertex pairs that need to be merged. (We guarantee that all
    // vertex pairs are separated by at least the merge radius in the output.)

    // First, we build the set of all the distinct vertices in the input.
    // We need to include the source and destination of every edge.
    Set<S2Point> vertices = Sets.newHashSet();
    for (Map.Entry<S2Point, Multiset<S2Point>> edge : edges.entrySet()) {
      vertices.add(edge.getKey());
      for (S2Point v : edge.getValue().elementSet()) {
        vertices.add(v);
      }
    }

    // Build a spatial index containing all the distinct vertices.
    for (S2Point p : vertices) {
      index.add(p);
    }

    // Next, we loop through all the vertices and attempt to grow a maximal
    // mergeable group starting from each vertex.
    Map<S2Point, S2Point> mergeMap = Maps.newHashMap();
    Stack<S2Point> frontier = new Stack<S2Point>();
    List<S2Point> mergeable = Lists.newArrayList();
    for (S2Point vstart : vertices) {
      // Skip any vertices that have already been merged with another vertex.
      if (mergeMap.containsKey(vstart)) {
        continue;
      }

      // Grow a maximal mergeable component starting from "vstart", the
      // canonical representative of the mergeable group.
      frontier.add(vstart);
      while (!frontier.isEmpty()) {
        // Pop the top frontier point and get all points nearby
        index.queryCap(frontier.pop(), mergeable);
        for (S2Point vnear : mergeable) {
          if (!vstart.equals(vnear)) {
            // Erase from the index any vertices that will be merged.  This
            // ensures that we won't try to merge the same vertex twice.
            index.remove(vnear);
            frontier.push(vnear);
            mergeMap.put(vnear, vstart);
          }
        }
      }
    }

    return mergeMap;
  }

  /**
   * Returns true if the given directed edge [v0 -> v1] is present in the
   * directed edge graph.
   */
  public boolean hasEdge(S2Point v0, S2Point v1) {
    Multiset<S2Point> vset = edges.get(v0);
    return vset == null ? false : vset.count(v1) > 0;
  }

  /** Uses the point index to help splice vertices that are near an edge onto the edge. */
  private void spliceEdges(PointIndex index) {
    // We keep a stack of unprocessed edges.  Initially all edges are
    // pushed onto the stack.
    List<S2Edge> pendingEdges = Lists.newArrayList();
    for (Map.Entry<S2Point, Multiset<S2Point>> edge : edges.entrySet()) {
      S2Point v0 = edge.getKey();
      for (S2Point v1 : edge.getValue().elementSet()) {
        // We only need to modify one copy of each undirected edge.
        if (!options.getUndirectedEdges() || v0.compareTo(v1) < 0) {
          pendingEdges.add(new S2Edge(v0, v1));
        }
      }
    }

    // For each edge, we check whether there are any nearby vertices that should
    // be spliced into it.  If there are, we choose one such vertex, split
    // the edge into two pieces, and iterate on each piece.
    while (!pendingEdges.isEmpty()) {
      // Must remove last edge before pushing new edges.
      S2Edge lastPair = pendingEdges.remove(pendingEdges.size() - 1);
      S2Point v0 = lastPair.getStart();
      S2Point v1 = lastPair.getEnd();

      // If we are XORing, edges may be erased before we get to them.
      if (options.getXorEdges() && !hasEdge(v0, v1)) {
        continue;
      }

      // Find nearest point to [v0,v1], or null if no point is within range.
      S2Point vmid = index.findNearbyPoint(v0, v1);
      if (vmid != null) {
        // Replace [v0,v1] with [v0,vmid] and [vmid,v1], and then add the new
        // edges to the set of edges to test for adjacent points.
        eraseEdge(v0, v1);
        if (addEdge(v0, vmid)) {
          pendingEdges.add(new S2Edge(v0, vmid));
        }
        if (addEdge(vmid, v1)) {
          pendingEdges.add(new S2Edge(vmid, v1));
        }
      }
    }
  }

  /**
   * A PointIndex is a cheap spatial index to help us find mergeable vertices.
   * Given a set of points, it can efficiently find all of the points within a
   * given search radius of an arbitrary query location. It is essentially just
   * a hash map from cell ids at a given fixed level to the set of points
   * contained by that cell id.
   *
   * <p>This class is not suitable for general use because it only supports
   * fixed-radius queries and has various special-purpose operations to avoid
   * the need for additional data structures.
   *
   * <p>This class is <strong>not thread-safe</strong>.
   */
  private class PointIndex {
    private double vertexRadius;
    private double edgeFraction;
    private int level;
    private final TreeMultimap<S2CellId, S2Point> delegate = TreeMultimap.create();
    private final List<S2CellId> ids = Lists.newArrayList();

    public PointIndex(double searchRadius, double edgeFraction) {
      this.vertexRadius = searchRadius;
      this.edgeFraction = edgeFraction;

      // We choose a cell level such that if dist(A,B) <= search_radius, the
      // S2CellId at that level containing A is a vertex neighbor of B (see
      // S2CellId.getVertexNeighbors). This turns out to be the highest
      // level such that a spherical cap (i.e. "disc") of the given radius
      // fits completely inside all cells at that level.
      this.level =
          Math.min(S2Projections.MIN_WIDTH.getMaxLevel(2 * searchRadius), S2CellId.MAX_LEVEL - 1);
    }

    /** Add a point to the index in each cell neighbor at the index level. */
    public void add(S2Point p) {
      S2CellId.fromPoint(p).getVertexNeighbors(level, ids);
      for (S2CellId id : ids) {
        delegate.put(id, p);
      }
      ids.clear();
    }

    public void remove(S2Point p) {
      S2CellId.fromPoint(p).getVertexNeighbors(level, ids);
      for (S2CellId id : ids) {
        delegate.remove(id, p);
      }
      ids.clear();
    }

    /**
     * Return the set of points whose distance to "axis" is less than
     * {@code searchRadius}.
     */
    public void queryCap(S2Point axis, List<S2Point> output) {
      output.clear();
      S2CellId id = S2CellId.fromPoint(axis).parent(level);
      for (Map.Entry<S2CellId, Collection<S2Point>> entry :
          delegate.asMap().tailMap(id).entrySet()) {
        if (entry.getKey().id() != id.id()) {
          break;
        }
        for (S2Point p : entry.getValue()) {
          if (axis.angle(p) < vertexRadius) {
            output.add(p);
          }
        }
      }
    }

    /**
     * Return a point whose distance from the edge (v0,v1) is less than
     * searchRadius, and which is not equal to v0 or v1. The current
     * implementation returns the closest such point.
     */
    public S2Point findNearbyPoint(S2Point v0, S2Point v1) {
      // Strategy: we compute a very cheap covering by approximating the edge as
      // two spherical caps, one around each endpoint, and then computing a
      // 4-cell covering of each one. We could improve the quality of the
      // covering by using some intermediate points along the edge as well.
      double length = v0.angle(v1);
      S2Point normal = S2.robustCrossProd(v0, v1);
      int queryLevel = Math.min(level, S2Projections.MIN_WIDTH.getMaxLevel(length));
      S2CellId.fromPoint(v0).getVertexNeighbors(queryLevel, ids);
      S2CellId.fromPoint(v1).getVertexNeighbors(queryLevel, ids);
      // Sort the cell ids so that we can skip duplicates in the loop below.
      Collections.sort(ids);
      double bestDistance = 2 * vertexRadius;
      S2Point nearby = null;
      for (int i = ids.size(); --i >= 0; ) {
        // Skip duplicates.
        if (i > 0 && ids.get(i - 1).id() == ids.get(i).id()) {
          continue;
        }
        S2CellId maxId = ids.get(i).rangeMax();
        for (Map.Entry<S2CellId, Collection<S2Point>> cellPoints :
            delegate.asMap().tailMap(ids.get(i).rangeMin()).entrySet()) {
          if (cellPoints.getKey().compareTo(maxId) > 0) {
            break;
          }
          for (S2Point p : cellPoints.getValue()) {
            if (!p.equals(v0) && !p.equals(v1)) {
              double dist = S2EdgeUtil.getDistance(p, v0, v1, normal).radians();
              if (dist < bestDistance) {
                bestDistance = dist;
                nearby = p;
              }
            }
          }
        }
      }
      ids.clear();
      if (bestDistance < edgeFraction * vertexRadius) {
        return nearby;
      } else {
        return null;
      }
    }
  }
}
