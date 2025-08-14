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

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.base.Strings;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2BuilderGraph.LabelFetcher;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.IdentityHashMap;
import java.util.List;

/**
 * A layer type that assembles edges (directed or undirected) into an S2Polygon. Returns an error if
 * the edges cannot be assembled into loops.
 *
 * <p>Note that by default, S2PolygonLayer does not validate the output polygon. To enable
 * validation, use {@link Options#setValidate(boolean)}, which will cause {@link
 * #build(S2BuilderGraph, S2Error)} to return false and set the error if the constructed polygon is
 * invalid.
 *
 * <p>If the input edges are directed, they must be oriented such that the polygon interior is to
 * the left of all edges. Directed edges are always preferred (see {@link S2Builder.EdgeType}).
 *
 * <p>Before the edges are assembled into loops, "sibling pairs" consisting of an edge and its
 * reverse edge are automatically removed. Such edge pairs represent zero-area degenerate regions,
 * which S2Polygon does not allow. If you need to build polygons with degeneracies, use {@link
 * S2LaxPolygonLayer} instead.
 *
 * <p>S2PolygonLayer is implemented such that if the input to S2Builder is a polygon and is not
 * modified, then the output has the same cyclic ordering of loop vertices and the same loop
 * ordering as the input polygon.
 *
 * <p>If the polygon has no edges, then the graph's IsFullPolygonPredicate is called to determine
 * whether the output polygon should be empty (containing no points) or full (containing all
 * points). This predicate can be specified as part of the S2Builder input geometry.
 */
public class S2PolygonLayer implements S2BuilderShapesLayer {
  /** The options for building the S2Polygon. */
  private final Options options;

  /**
   * The output polygon. Either provided by {@link #setPolygon(S2Polygon)} before build() is called,
   * or constructed by build().
   */
  private S2Polygon polygon = null;

  /**
   * If provided to the constructor, maps shape edge ids to label set ids after build(). Otherwise
   * null.
   */
  private final IntVector labelSetIds;

  /**
   * When labelSetIds are required, they are first fetched by loop, and then reordered to match the
   * constructed polygon loop ordering. These ArrayLists are temporary storage for that process.
   */
  private final ArrayList<IntVector> labelSetIdsForLoops;

  private final ArrayList<IntVector> newLabelSetIdsForLoops;

  /**
   * If provided to the constructor, labelSetLexicon provides label sets from label set ids after
   * build(). Otherwise null.
   */
  private final IdSetLexicon labelSetLexicon;

  /** The label fetcher is reused for each build() call, but constructed only if needed. */
  private LabelFetcher fetcher = null;

  /** This list of loops is reused for each build() call. */
  private final ArrayList<S2Loop> loops = new ArrayList<>();

  /** This list of loops is reused for each build() call. */
  private final ArrayList<int[]> edgeLoops = new ArrayList<>();

  /** Constructs an S2PolygonLayer using default {@link Options}. */
  public S2PolygonLayer() {
    this.options = new Options();
    labelSetIds = null;
    labelSetIdsForLoops = null;
    newLabelSetIdsForLoops = null;
    labelSetLexicon = null;
  }

  /** Constructs an S2PolygonLayer with the provided options. */
  public S2PolygonLayer(Options options) {
    Preconditions.checkNotNull(options);
    this.options = options;
    labelSetIds = null;
    labelSetIdsForLoops = null;
    newLabelSetIdsForLoops = null;
    labelSetLexicon = null;
  }

  /**
   * Constructs an S2PolygonLayer with the provided options. Each time build() is called, the
   * provided IdSetLexicon and IntVector will be filled in so that labels attached to the input
   * edges can be obtained for the edges of the output polygon that the input edges were snapped to.
   */
  public S2PolygonLayer(Options options, IdSetLexicon labelSetLexicon, IntVector labelSetIds) {
    Preconditions.checkNotNull(options);
    Preconditions.checkNotNull(labelSetLexicon);
    Preconditions.checkNotNull(labelSetIds);
    this.options = options;
    this.labelSetIds = labelSetIds;
    this.labelSetLexicon = labelSetLexicon;

    labelSetIdsForLoops = new ArrayList<>();
    newLabelSetIdsForLoops = new ArrayList<>();
    fetcher = new LabelFetcher();
  }

  /**
   * Provides a polygon that will be initialized with the built polygon. If not called, a new
   * polygon will be constructed and can be obtained with {@link #getPolygon()}.
   */
  public void setPolygon(S2Polygon polygon) {
    this.polygon = polygon;
  }

  /**
   * Returns the new output polygon after build() has been called, if build() returned true. If
   * {@link #setPolygon(S2Polygon)} was called before build(), this returns the same polygon.
   */
  public S2Polygon getPolygon() {
    Preconditions.checkNotNull(polygon, "Must call build() first.");
    return polygon;
  }

  @Override
  public String toString() {
    return "S2PolygonLayer with options " + options;
  }

  @Override
  public Iterable<S2Polygon.Shape> shapes() {
    return ImmutableList.of(getPolygon().shape());
  }

  /** Returns the mapping from loop number -> edge id -> IdSet id after the layer has been built. */
  @VisibleForTesting
  ArrayList<IntVector> getLabelSetIdsForLoops() {
    Preconditions.checkNotNull(labelSetIdsForLoops, "Must call build() first.");
    return labelSetIdsForLoops;
  }

  @Override
  public GraphOptions graphOptions() {
    // Prevent degenerate edges and sibling edge pairs. There should not be any duplicate edges
    // if the input is valid, but if there are then we keep them since this tends to produce more
    // comprehensible errors.
    return new GraphOptions(
        options.edgeType(),
        GraphOptions.DegenerateEdges.DISCARD,
        GraphOptions.DuplicateEdges.KEEP,
        GraphOptions.SiblingPairs.DISCARD);
  }

  @Override
  public boolean build(S2BuilderGraph g, S2Error error) {
    if (labelSetIdsForLoops != null) {
      labelSetIds.clear();
      labelSetIdsForLoops.clear();
      labelSetLexicon.clear();
    }

    // It's tricky to compute the edge labels for S2Polygons because the S2Polygon.init methods
    // can reorder and/or invert the loops. We handle this by remembering the original vector
    // index of each loop and whether or not the loop contained S2.origin(). By comparing this with
    // the final S2Polygon loops we can fix up the edge labels appropriately.
    //
    // See {@link #appendEdgeLabels} and {@link #reorderEdgeLabels} for the details.
    LoopMap loopMap = new LoopMap();
    loops.clear();
    edgeLoops.clear();
    if (polygon == null) {
      polygon = new S2Polygon();
    }
    if (g.numEdges() == 0) {
      // The polygon is either full or empty.
      ArrayList<S2Loop> loops = new ArrayList<>();
      if (g.isFullPolygonPredicate().test(g)) {
        loops.add(S2Loop.full());
      }
      polygon.init(loops);
    } else if (g.options().edgeType() == S2Builder.EdgeType.DIRECTED) {
      if (!g.getDirectedLoops(S2BuilderGraph.LoopType.SIMPLE, edgeLoops, error)) {
        return false;
      }
      appendS2Loops(g, edgeLoops, loops);
      appendEdgeLabels(g, edgeLoops);
      edgeLoops.clear(); // Release memory
      initLoopMap(loops, loopMap);
      polygon.initOriented(loops);
    } else {
      ArrayList<S2BuilderGraph.UndirectedComponent> components = new ArrayList<>();
      if (!g.getUndirectedComponents(S2BuilderGraph.LoopType.SIMPLE, components, error)) {
        return false;
      }
      // It doesn't really matter which complement of each component we use, since below we
      // normalize all the loops so that they enclose at most half of the sphere (to ensure that
      // the loops can always be nested).
      //
      // The only reason to prefer one over the other is that when there are multiple loops that
      // touch, only one of the two complements matches the structure of the input loops.
      // getUndirectedComponents() tries to ensure that this is always complement 0 of each
      // component.
      for (S2BuilderGraph.UndirectedComponent component : components) {
        appendS2Loops(g, component.getComplement(0), loops);
        appendEdgeLabels(g, component.getComplement(0));
      }
      components.clear(); // Release memory
      initLoopMap(loops, loopMap);
      for (S2Loop loop : loops) {
        loop.normalize();
      }
      polygon.initNested(loops);
    }
    reorderEdgeLabels(loopMap);

    if (labelSetIdsForLoops != null) {
      // labelSetIdsForLoops maps (loop number -> loop edgeId -> IdSetId). We want a flat map from
      // from shape edgeId -> IdSetId.
      for (int loopId = 0; loopId < polygon.numLoops(); loopId++) {
        S2Loop loop = polygon.loop(loopId);
        for (int offset = 0; offset < loop.numVertices(); offset++) {
          // Iterate the edges of inverted loops in reverse.
          labelSetIds.add(labelSetIdsForLoops.get(loopId).get(orientedEdge(loop, offset)));
        }
      }
    }

    if (options.validate()) {
      return !polygon.findValidationError(error);
    }

    return true;
  }

  /**
   * If the loop is a hole and has N three vertices & edges, with N at least 3, then the first N-2
   * edges have to be visited in reverse order. See the comments on reorderEdgeLabels() below.
   */
  private static int orientedEdge(S2Loop loop, int offset) {
    if (loop.isHole() && loop.numVertices() > 2 && offset < (loop.numVertices() - 1)) {
      offset = loop.numVertices() - offset - 2;
    }
    return offset;
  }

  /**
   * Given an S2BuilderGraph and a list of edge loops as arrays of edge ids in the graph, fills in
   * the provided 'outLoops' with a new S2Loop for each provided int[]. The S2Loops are not
   * validated.
   */
  private void appendS2Loops(S2BuilderGraph g, List<int[]> edgeLoops, List<S2Loop> outLoops) {
    for (int[] edgeLoop : edgeLoops) {
      S2Point[] vertices = new S2Point[edgeLoop.length];
      int j = 0;
      for (int edgeId : edgeLoop) {
        vertices[j++] = g.vertex(g.edgeSrcId(edgeId));
      }
      S2Loop l = new S2Loop(vertices);
      outLoops.add(l);
    }
  }

  /**
   * Given an S2BuilderGraph and a list of edge loops as arrays of edge ids in the graph, fills in
   * the internal 'loopLabelSetIds'. When complete, then for a given edge loop L with edge ids
   * (e1..en), {@code int id = loopLabelSetIds.get(L).get(ei)} is the id of the set of labels on
   * that loop edge, which are stored in the labelSetLexicon.
   */
  private void appendEdgeLabels(S2BuilderGraph g, List<int[]> edgeLoops) {
    if (labelSetIdsForLoops == null) {
      return;
    }
    fetcher.init(g, options.edgeType());

    IntVector edgeLabelSet = new IntVector(); // Temporary storage for labels.
    for (int[] edgeLoop : edgeLoops) {
      // A new map from loop edge to label set ids is constructed for each loop.
      IntVector loopLabelSetIds = new IntVector();
      loopLabelSetIds.ensureCapacity(edgeLoop.length);

      for (int edgeId : edgeLoop) {
        // Get the labels for this edge of the loop.
        fetcher.fetch(edgeId, edgeLabelSet);
        // Store the labels for this edge in the labelSetLexicon, and store the id of the stored
        // label set for this edge.
        loopLabelSetIds.add(labelSetLexicon.add(edgeLabelSet));
      }

      // Add the vector mapping edge id to labelset id for this loop to the list.
      labelSetIdsForLoops.add(loopLabelSetIds);
    }
  }

  /** For some loop, an integer index into 'loopLabelSetIds' and a boolean 'containsOrigin'. */
  private static class LoopData {
    final int index;
    final boolean containsOrigin;

    LoopData(int index, boolean containsOrigin) {
      this.index = index;
      this.containsOrigin = containsOrigin;
    }

    static LoopData of(int index, boolean containsOrigin) {
      return new LoopData(index, containsOrigin);
    }
  }

  /**
   * A LoopMap maps each S2Loop to an index into 'loopLabelSetIds' and a boolean 'containsOrigin'.
   * It is an IdentityHashMap because loops may be modified (inverted) without changing identity.
   */
  private static class LoopMap extends IdentityHashMap<S2Loop, LoopData> {}

  /**
   * Fills in the given 'loopMap', mapping each of the given S2Loops to a pair of its index in the
   * given 'loops' list, and if it contains the S2.origin.
   */
  private void initLoopMap(List<S2Loop> loops, LoopMap loopMap) {
    if (labelSetIdsForLoops == null) {
      return;
    }
    for (int i = 0; i < loops.size(); ++i) {
      S2Loop loop = loops.get(i);
      loopMap.put(loop, LoopData.of(i, loop.containsOrigin()));
    }
  }

  /**
   * Updates loopLabelSetIds with reordered maps from loop edge ids to the associated label set ids,
   * after polygon loops have been possibly inverted and reordered.
   */
  private void reorderEdgeLabels(LoopMap loopMap) {
    if (labelSetIdsForLoops == null) {
      return;
    }

    // The new loopLabelSetIds.
    newLabelSetIdsForLoops.clear();
    newLabelSetIdsForLoops.ensureCapacity(labelSetIdsForLoops.size());

    for (int i = 0; i < polygon.numLoops(); ++i) {
      S2Loop loop = polygon.loop(i);

      // Get the old (index into loopLabelSetIds, containsOrigin) for loop 'i'.
      LoopData old = loopMap.get(loop);
      Preconditions.checkState(old != null, "Loop #%s", i);

      // Get the IntVector mapping edgeId->labelSetId for this loop.
      IntVector edgeLabelSets = labelSetIdsForLoops.get(old.index);
      Preconditions.checkState(edgeLabelSets != null, "Loop #%s with old index %s", i, old.index);

      // If the loop was inverted, the edgeId to labelSetId ordering must be updated.
      if (loop.containsOrigin() != old.containsOrigin) {
        // S2Loop.invert() reverses the order of the vertices, which leaves the last edge unchanged.
        // For example, the loop ABCD (with edges AB, BC, CD, DA) becomes the loop DCBA (with edges
        // DC, CB, BA, AD).
        if (edgeLabelSets.size() > 2) {
          edgeLabelSets.reverse(0, edgeLabelSets.size() - 1);
        }
      }
      // Add the new vector to the list.
      newLabelSetIdsForLoops.add(edgeLabelSets);
    }

    labelSetIdsForLoops.clear();
    labelSetIdsForLoops.addAll(newLabelSetIdsForLoops);
  }

  /** Options for the S2PolygonLayer. */
  public static class Options {
    private S2Builder.EdgeType edgeType;
    private boolean validate;

    /** Options constructor that sets default options: directed edges and no polygon validation. */
    public Options() {
      edgeType = S2Builder.EdgeType.DIRECTED;
      validate = false;
    }

    /** Options constructor that specifies the edge type but otherwise sets the defaults. */
    public Options(S2Builder.EdgeType edgeType) {
      Preconditions.checkNotNull(edgeType);
      this.edgeType = edgeType;
      validate = false;
    }

    @Override
    public String toString() {
      return Strings.lenientFormat("EdgeType %s, Validate=%s", edgeType, validate);
    }

    /**
     * Sets the option for considering input edges provided to S2Builder as directed or undirected.
     * Directed edges should be used whenever possible to avoid ambiguity.See {@link
     * S2Builder.EdgeType} for more details.
     *
     * <p>If the input edges are directed, they should be oriented so that the polygon interior is
     * to the left of all edges. This means that for a polygon with holes, the outer loops
     * ("shells") should be directed counter-clockwise while the inner loops ("holes") should be
     * directed clockwise. Note that {@link S2Builder#addPolygon(S2Polygon)} does this
     * automatically.
     *
     * <p>The default is EdgeType.DIRECTED.
     */
    @CanIgnoreReturnValue
    public Options setEdgeType(S2Builder.EdgeType edgeType) {
      Preconditions.checkNotNull(edgeType);
      this.edgeType = edgeType;
      return this;
    }

    /** Returns the edge type. */
    public S2Builder.EdgeType edgeType() {
      return edgeType;
    }

    /**
     * Sets whether the output polygon should be validated. If true, findValidationError() is called
     * on the output polygon. Any error found will be returned by {@link S2Builder#build(S2Error)}.
     * The default is false.
     */
    @CanIgnoreReturnValue
    public Options setValidate(boolean validate) {
      this.validate = validate;
      return this;
    }

    /** Returns true if the output polygon should be validated. */
    public boolean validate() {
      return validate;
    }
  }
}
