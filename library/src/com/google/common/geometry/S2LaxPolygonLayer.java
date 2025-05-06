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
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2BuilderGraph.EdgeList;
import com.google.common.geometry.S2BuilderGraph.LabelFetcher;
import com.google.common.geometry.S2BuilderGraph.LoopType;
import com.google.common.geometry.S2PolygonDegeneracyFinder.PolygonDegeneracyList;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.List;

/**
 * A layer type that assembles edges (directed or undirected) into a S2LaxPolygonShape. Returns an
 * error if the edges cannot be assembled into loops.
 *
 * <p>If the input edges are directed, they must be oriented such that the polygon interior is to
 * the left of all edges. Directed edges are always preferred (see S2Builder.EdgeType).
 *
 * <p>LaxPolygonLayer is implemented such that if the input to S2Builder is a polygon and is not
 * modified, then the output has the same cyclic ordering of loop vertices and the same loop
 * ordering as the input polygon.
 *
 * <p>If the given edge graph is degenerate (i.e., it consists entirely of degenerate edges and
 * sibling pairs), then the IsFullPolygonPredicate associated with the edge graph is called to
 * determine whether the output polygon should be empty (possibly with degenerate shells) or full
 * (possibly with degenerate holes). This predicate can be specified as part of the S2Builder input
 * geometry.
 *
 * <p>Note that support for undirected edges, and validation of the constructed polygons is not yet
 * implemented. See b/281698554.
 */
@SuppressWarnings("Assertion")
public class S2LaxPolygonLayer implements S2BuilderShapesLayer {
  /** Options for building the lax polygon. */
  private final Options options;

  /** The new polygon constructed during each call to build(). */
  private S2LaxPolygonShape polygon = null;

  /**
   * Maps loop to edgeId to IdSetId. The IdSet Id can be used with the labelSetLexicon to get the
   * labels for that loop edge.
   */
  private ArrayList<IntVector> labelSetIds = new ArrayList<>();

  /**
   * If provided to the constructor, labelSetLexicon provides label sets from label set ids after
   * build(). Otherwise null.
   */
  private final IdSetLexicon labelSetLexicon;

  /** The label fetcher is reused for each build() call, but constructed only if needed. */
  private LabelFetcher fetcher = null;

  /** Constructs an S2LaxPolygonLayer using default {@link Options}. */
  public S2LaxPolygonLayer() {
    this.options = new Options();
    labelSetIds = null;
    labelSetLexicon = null;
  }

  /** Constructs an S2LaxPolygonLayer with the provided options. */
  public S2LaxPolygonLayer(Options options) {
    this.options = options;
    labelSetIds = null;
    labelSetLexicon = null;
  }

  /**
   * Constructs an S2LaxPolygonLayer with the provided options. Each time build() is called, the
   * provided IdSetLexicon and IntVector will be filled in so that labels attached to the input
   * edges can be obtained for the edges of the output polygon that the input edges were snapped to.
   */
  public S2LaxPolygonLayer(
      Options options, IdSetLexicon labelSetLexicon, ArrayList<IntVector> labelSetIds) {
    Preconditions.checkNotNull(options);
    Preconditions.checkNotNull(labelSetLexicon);
    Preconditions.checkNotNull(labelSetIds);
    this.options = options;
    this.labelSetIds = labelSetIds;
    this.labelSetLexicon = labelSetLexicon;
  }

  /** Returns the output lax polygon after build() has been called, if build() returned true. */
  public S2LaxPolygonShape getPolygon() {
    Preconditions.checkNotNull(polygon, "Must call build() first.");
    return polygon;
  }

  @Override
  public String toString() {
    return "S2LaxPolygonLayer with options " + options;
  }

  @Override
  public Iterable<S2Shape> shapes() {
    return ImmutableList.of(getPolygon());
  }

  /** Returns the mapping from loop number to edgeId to IdSetId after the layer has been built. */
  @VisibleForTesting
  ArrayList<IntVector> getLabelSetIdsForLoops() {
    Preconditions.checkNotNull(labelSetIds, "Must call build() first.");
    return labelSetIds;
  }

  @Override
  public GraphOptions graphOptions() {
    if (options.degenerateBoundaries() == DegenerateBoundaries.DISCARD) {
      // There should not be any duplicate edges, but if there are then we keep them since this
      // yields more comprehensible error messages.
      return new GraphOptions(
          options.edgeType(), DegenerateEdges.DISCARD, DuplicateEdges.KEEP, SiblingPairs.DISCARD);
    } else {
      // Keep at most one copy of each sibling pair and each isolated vertex.
      return new GraphOptions(
          options.edgeType(),
          DegenerateEdges.DISCARD_EXCESS,
          DuplicateEdges.KEEP,
          SiblingPairs.DISCARD_EXCESS);
    }
  }

  /**
   * Returns true if building an S2LaxPolylineShape succeeded. Otherwise returns false and fills in
   * 'error'.
   */
  @Override
  public boolean build(S2BuilderGraph g, S2Error error) {
    if (labelSetIds != null) {
      labelSetIds.clear();
      labelSetLexicon.clear();
    }
    if (g.options().edgeType() == EdgeType.DIRECTED) {
      return buildDirected(g, error);
    } else {
      error.init(S2Error.Code.UNIMPLEMENTED, "Undirected edges not supported yet");
      return false;
    }
  }

  /** Returns all edges of "g" except for those identified by "edgesToDiscard". */
  private static void discardEdges(
      S2BuilderGraph g,
      IntVector edgesToDiscard,
      EdgeList newEdges,
      IntVector newInputEdgeIdSetIds) {
    assert edgesToDiscard.isSorted();
    newEdges.clear();
    newInputEdgeIdSetIds.clear();
    newEdges.ensureCapacity(g.numEdges());
    newInputEdgeIdSetIds.ensureCapacity(g.numEdges());

    int i = 0; // Index into edgesToDiscard
    for (int edgeId = 0; edgeId < g.numEdges(); ++edgeId) {
      if (i != edgesToDiscard.size() && edgeId == edgesToDiscard.get(i)) {
        ++i;
      } else {
        newEdges.copyEdge(g.edges(), edgeId);
        newInputEdgeIdSetIds.add(g.inputEdgeIdSetId(edgeId));
      }
    }
    assert i == edgesToDiscard.size();
  }

  private static void maybeAddFullLoop(S2BuilderGraph g, ArrayList<List<S2Point>> loops) {
    if (g.isFullPolygonPredicate().test(g)) {
      loops.add(ImmutableList.of()); // Full loop.
    }
  }

  /** Returns true if successful, otherwise sets S2Error. */
  private boolean buildDirected(S2BuilderGraph g, S2Error error) {
    // Some cases are implemented by constructing a new graph with certain degenerate edges removed
    // (overwriting "g"). "newEdges" is where the edges for the new graph are stored.
    EdgeList newEdges = new EdgeList();
    IntVector newInputEdgeIdSetIds = new IntVector();
    ArrayList<List<S2Point>> loops = new ArrayList<>();
    DegenerateBoundaries degenerateBoundaries = options.degenerateBoundaries();

    switch (degenerateBoundaries) {
      case DISCARD:
        // This is the easiest case, since there are no degeneracies.
        if (g.numEdges() == 0) {
          maybeAddFullLoop(g, loops);
        }
        break;
      case KEEP:
        // S2LaxPolygonShape doesn't need to distinguish degenerate shells from holes except when
        // the entire graph is degenerate, in which case we need to decide whether it represents an
        // empty polygon possibly with degenerate shells, or a full polygon possibly with degenerate
        // holes.
        if (S2PolygonDegeneracyFinder.isFullyDegenerate(g)) {
          maybeAddFullLoop(g, loops);
        }
        break;
      default:
        // For DISCARD_SHELLS and DISCARD_HOLES we first determine whether any degenerate loops of
        // the given type exist, and if so we construct a new graph with those edges removed
        // (overwriting "g").
        boolean discardHoles = degenerateBoundaries == DegenerateBoundaries.DISCARD_HOLES;
        PolygonDegeneracyList degeneracies = S2PolygonDegeneracyFinder.findPolygonDegeneracies(g);
        if (degeneracies.size() == g.numEdges()) {
          if (degeneracies.size() == 0) {
            maybeAddFullLoop(g, loops);
          } else if (degeneracies.isHole(0)) {
            loops.add(ImmutableList.of()); // Full loop.
          }
        }
        IntVector edgesToDiscard = new IntVector();
        for (int i = 0; i < degeneracies.size(); ++i) {
          if (degeneracies.isHole(i) == discardHoles) {
            edgesToDiscard.add(degeneracies.edgeId(i));
          }
        }
        if (!edgesToDiscard.isEmpty()) {
          // Construct a new graph that discards the unwanted edges.
          edgesToDiscard.sort();
          discardEdges(g, edgesToDiscard, newEdges, newInputEdgeIdSetIds);
          g =
              new S2BuilderGraph(
                  g.options(),
                  g.vertices(),
                  newEdges,
                  newInputEdgeIdSetIds,
                  g.inputEdgeIdSetLexicon(),
                  g.labelSetIds(),
                  g.labelSetLexicon(),
                  g.isFullPolygonPredicate());
        }
    }

    ArrayList<int[]> edgeLoops = new ArrayList<>();
    if (!g.getDirectedLoops(LoopType.CIRCUIT, edgeLoops, error)) {
      return error.ok();
    }

    appendPolygonLoops(g, edgeLoops, loops);
    appendEdgeLabels(g, edgeLoops);
    polygon = S2LaxPolygonShape.create(loops);
    return true;
  }

  private void appendPolygonLoops(
      S2BuilderGraph g, ArrayList<int[]> edgeLoops, ArrayList<List<S2Point>> loops) {
    for (int[] edgeLoop : edgeLoops) {
      ArrayList<S2Point> vertices = new ArrayList<>(edgeLoop.length);
      for (int edgeId : edgeLoop) {
        vertices.add(g.vertex(g.edges().getSrcId(edgeId)));
      }
      loops.add(vertices);
    }
  }

  private void appendEdgeLabels(S2BuilderGraph g, List<int[]> edgeLoops) {
    if (labelSetIds == null) {
      return;
    }

    if (fetcher == null) {
      fetcher = new LabelFetcher(g, options.edgeType());
    } else {
      fetcher.init(g, options.edgeType());
    }

    // Temporary storage for the labels on one edge.
    IntVector labels = new IntVector();
    for (int[] edgeLoop : edgeLoops) {
      IntVector loopLabelSetIds = IntVector.ofCapacity(edgeLoop.length);
      for (int edgeId : edgeLoop) {
        // TODO(torrey): This pattern occurs in many layers. Consider making it possible for
        // the fetcher to add directly to an IdSetLexicon and IntVector.
        fetcher.fetch(edgeId, labels);
        loopLabelSetIds.add(labelSetLexicon.add(labels));
      }
      labelSetIds.add(loopLabelSetIds);
    }
  }

  /**
   * Specifies whether degenerate boundaries should be discarded or kept. (A degenerate boundary
   * consists of either a sibling edge pair or an edge from a vertex to itself.) Optionally,
   * degenerate boundaries may be kept only if they represent shells, or only if they represent
   * holes.
   *
   * <p>This option is useful for normalizing polygons with various boundary conditions. For
   * example, DISCARD_HOLES can be used to normalize closed polygons (those that include their
   * boundary), since degenerate holes do not affect the set of points contained by such polygons.
   * Similarly, DISCARD_SHELLS can be used to normalize polygons with open boundaries. DISCARD is
   * used to normalize polygons with semi-open boundaries (since degenerate loops do not affect
   * point containment in that case), and finally KEEP is useful for working with any type of
   * polygon where degeneracies are assumed to contain an infinitesmal interior. (This last model is
   * the most useful for working with simplified geometry, since it maintains the closest fidelity
   * to the original geometry.)
   */
  public enum DegenerateBoundaries {
    DISCARD,
    DISCARD_HOLES,
    DISCARD_SHELLS,
    KEEP
  }

  /**
   * Options for the S2LaxPolygonLayer are if the edges are directed or undirected, and how to
   * handle degenerate loop boundaries.
   */
  public static class Options {
    private EdgeType edgeType;
    private DegenerateBoundaries degenerateBoundaries;

    /**
     * Constructor that uses the default options, which are EdgeType.DIRECTED and
     * DegenerateBoundaries.KEEP.
     */
    public Options() {
      edgeType = EdgeType.DIRECTED;
      degenerateBoundaries = DegenerateBoundaries.KEEP;
    }

    /** Constructor that specifies the edge type. */
    public Options(EdgeType edgeType) {
      Preconditions.checkNotNull(edgeType);
      this.edgeType = edgeType;
    }

    /** Indicates whether the input edges provided to S2Builder are directed or undirected. */
    public EdgeType edgeType() {
      return edgeType;
    }

    /**
     * Sets the option for considering input edges provided to S2Builder as directed or undirected.
     * Directed edges should be used whenever possible to avoid ambiguity. See {@link
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
    public Options setEdgeType(EdgeType edgeType) {
      Preconditions.checkNotNull(edgeType);
      this.edgeType = edgeType;
      return this;
    }

    /** Indicates whether degenerate boundaries should be discarded or kept. */
    public DegenerateBoundaries degenerateBoundaries() {
      return degenerateBoundaries;
    }

    /**
     * Sets the option for keeping or discarding degenerate boundaries. See {@link
     * DegenerateBoundaries} for documentation of the options.
     *
     * <p>The default is DegenerateBoundaries.KEEP.
     */
    @CanIgnoreReturnValue
    public Options setDegenerateBoundaries(DegenerateBoundaries degenerateBoundaries) {
      Preconditions.checkNotNull(degenerateBoundaries);
      this.degenerateBoundaries = degenerateBoundaries;
      return this;
    }

    @Override
    public String toString() {
      return "EdgeType "
          + edgeType.name()
          + ", DegenerateBoundaries "
          + degenerateBoundaries.name();
    }
  }
}
