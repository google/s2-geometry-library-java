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
import com.google.common.base.Strings;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2BuilderGraph.LabelFetcher;
import com.google.common.geometry.S2BuilderGraph.PolylineBuilder;
import com.google.common.geometry.S2BuilderGraph.PolylineType;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.List;

/**
 * A layer type that assembles edges (directed or undirected) into multiple S2Polylines. The {@link
 * #build(S2BuilderGraph, S2Error)} method returns false and sets the provided S2Error if S2Builder
 * finds any problems with the input edges or (if validation is enabled) an invalid S2Polyline is
 * generated.
 *
 * <p>Note that by default, S2PolylineVectorLayer does not validate the output polylines. To enable
 * validation, use {@link Options#setValidate(boolean)}, which will cause {@link
 * #build(S2BuilderGraph, S2Error)} to return false and set the error if any output polyline is
 * invalid.
 *
 * <p>Duplicate edges are handled correctly (e.g., if a polyline backtracks on itself, or loops
 * around and retraces some of its previous edges.) The implementation attempts to preserve the
 * order of input edges whenever possible, so that if the input is a polyline and it is not modified
 * by S2Builder, then the output will be the same polyline even if the polyline forms a loop.
 * However, note that this is not guaranteed when undirected edges are used: for example, if the
 * input consists of a single undirected edge, then either directed edge may be returned.
 */
public class S2PolylineVectorLayer implements S2BuilderShapesLayer {
  /** The options for building S2Polylines. */
  private final Options options;

  /** The set of polylines constructed during each call to build(). */
  private List<S2Polyline> polylines = null;

  /**
   * If provided to the constructor, labelSetIds maps shape edge ids to label set ids after build(),
   * with each IntVector in the labelSetIds list providing label set ids for the corresponding
   * S2Polyline in the polylines list. Otherwise null.
   */
  private final List<IntVector> labelSetIds;

  /**
   * If provided to the constructor, the labelSetLexicon provides label sets from label set ids
   * after build(), for all the constructed polylines. Otherwise null.
   */
  private final IdSetLexicon labelSetLexicon;

  /** The polyline builder is reused for each build() call. */
  private final PolylineBuilder polylineBuilder = new PolylineBuilder();

  /** The label fetcher is reused for each build() call, but constructed only if needed. */
  private LabelFetcher fetcher = null;

  /** Constructs an S2PolylineVectorLayer using default {@link Options}. */
  public S2PolylineVectorLayer() {
    this.options = new Options();
    labelSetIds = null;
    labelSetLexicon = null;
  }

  /** Constructs an S2PolylineVectorLayer with the provided options. */
  public S2PolylineVectorLayer(Options options) {
    this.options = options;
    labelSetIds = null;
    labelSetLexicon = null;
  }

  /**
   * Constructs an S2PolylineVectorLayer with the provided options. Each time build() is called, the
   * provided IdSetLexicon and List of IntVectors will be filled in so that the labels attached to
   * the input edges can be obtained for the edges of the output polylines that the input edges
   * where snapped to.
   */
  public S2PolylineVectorLayer(
      Options options, IdSetLexicon labelSetLexicon, List<IntVector> labelSetIds) {
    Preconditions.checkNotNull(options);
    Preconditions.checkNotNull(labelSetLexicon);
    Preconditions.checkNotNull(labelSetIds);
    this.options = options;
    this.labelSetIds = labelSetIds;
    this.labelSetLexicon = labelSetLexicon;
  }

  @Override
  public String toString() {
    return "S2PolylineVectorLayer with Options " + options;
  }

  /** Returns the output polylines after build() has been called, if build() returned true. */
  public List<S2Polyline> getPolylines() {
    Preconditions.checkNotNull(polylines, "Must call build() first.");
    return polylines;
  }

  @Override
  public Iterable<S2Polyline> shapes() {
    return getPolylines();
  }

  @Override
  public GraphOptions graphOptions() {
    return new GraphOptions(
        options.edgeType(),
        GraphOptions.DegenerateEdges.DISCARD_EXCESS,
        options.duplicateEdges(),
        options.siblingPairs());
  }

  @Override
  public boolean build(S2BuilderGraph g, S2Error error) {
    // Use the PolylineBuilder to get the paths or walks, as arrays of edge ids.
    polylineBuilder.init(g);
    List<int[]> edgePolylines =
        (options.polylineType() == PolylineType.PATH)
            ? polylineBuilder.buildPaths()
            : polylineBuilder.buildWalks();

    // The polylines to be returned.
    polylines = new ArrayList<>(edgePolylines.size());
    // Temporary storage for labels.
    IntVector labels = null;

    if (labelSetIds != null) {
      if (fetcher == null) {
        fetcher = new LabelFetcher();
      }
      fetcher.init(g, options.edgeType());
      // TODO(torrey): Reuse the IntVectors in labelSetIds, along with reworking LabelFetcher to
      // fill in the lexicon and label set ids directly.
      labelSetIds.clear();
      labelSetLexicon.clear();
      labels = new IntVector();
    }

    for (int[] edgePolyline : edgePolylines) {
      // Create a new array to hold the vertices, and fill it with the polyline vertices.
      S2Point[] vertices = new S2Point[edgePolyline.length + 1];
      vertices[0] = g.vertex(g.edgeSrcId(edgePolyline[0]));
      int i = 1;
      for (int edgeId : edgePolyline) {
        vertices[i++] = g.vertex(g.edgeDstId(edgeId));
      }

      // S2Polyline represents a degenerate edge as a polyline with a single vertex, and does not
      // support adjacent identical vertices, so that case must be detected and handled specially.
      if (vertices.length == 2 && vertices[0].equalsPoint(vertices[1])) {
        vertices = new S2Point[] { vertices[0] };
      }
      // Create a polyline. This constructor takes ownership of the array of vertices.
      S2Polyline polyline = new S2Polyline(vertices);
      if (options.validate()) {
        if (polyline.findValidationError(error)) {
          return false;
        }
      }
      polylines.add(polyline);

      if (labelSetIds != null) {
        IntVector polylineLabels = new IntVector();
        for (int edgeId : edgePolyline) {
          fetcher.fetch(edgeId, labels);
          polylineLabels.add(labelSetLexicon.add(labels));
        }
        labelSetIds.add(polylineLabels);
      }
    }
    return true;
  }

  /** Options for the S2PolylineVectorLayer. */
  public static class Options {
    private EdgeType edgeType;
    private PolylineType polylineType;
    private DuplicateEdges duplicateEdges;
    private SiblingPairs siblingPairs;
    private boolean validate;

    /**
     * Options constructor that sets the default options: directed edges forming paths, keeping
     * duplicate edges and sibling pairs, no validation, and no edge labels.
     */
    public Options() {
      edgeType = EdgeType.DIRECTED;
      polylineType = PolylineType.PATH;
      duplicateEdges = DuplicateEdges.KEEP;
      siblingPairs = SiblingPairs.KEEP;
      validate = false;
    }

    /** Options constructor that specifies the edge type but otherwise uses the defaults. */
    public Options(EdgeType edgeType) {
      this.edgeType = edgeType;
      polylineType = PolylineType.PATH;
      duplicateEdges = DuplicateEdges.KEEP;
      siblingPairs = SiblingPairs.KEEP;
      validate = false;
    }

    @Override
    public String toString() {
      return Strings.lenientFormat(
          "EdgeType %s, PolylineType %s, DuplicateEdges %s, SiblingPairs %s, Validate=%s",
          edgeType, polylineType, duplicateEdges, siblingPairs, validate);
    }

    /**
     * Sets the option for treating input edges provided to S2Builder as directed or undirected.
     * Directed edges should be used whenever possible to avoid ambiguity. See {@link
     * S2Builder.EdgeType} for more details.
     *
     * <p>The implementation attempts to preserve the structure of directed input edges whenever
     * possible, so that if the input is a list of disjoint polylines and none of them need to be
     * modified then the output will be the same polylines in the same order. With undirected edges,
     * there are no such guarantees.
     *
     * <p>The default is EdgeType.DIRECTED.
     */
    @CanIgnoreReturnValue
    public Options setEdgeType(EdgeType edgeType) {
      this.edgeType = edgeType;
      return this;
    }

    /** Returns the edge type. */
    public EdgeType edgeType() {
      return edgeType;
    }

    /**
     * Controls how polylines are constructed. If the polyline type is PATH, then only vertices of
     * indegree and outdegree 1 (or degree 2 in the case of undirected edges) will appear in the
     * interior of polylines. Use this option if you want to split polylines into separate pieces
     * whenever they self-intersect or cross each other.
     *
     * <p>If "polylineType" is WALK, then each polyline will be as long as possible. Polylines may
     * pass through the same vertex or even the same edge multiple times (if duplicate edges are
     * present).
     *
     * <p>The default is PolylineType.PATH.
     */
    @CanIgnoreReturnValue
    public Options setPolylineType(PolylineType polylineType) {
      this.polylineType = polylineType;
      return this;
    }

    /** Returns the polyline type. */
    public PolylineType polylineType() {
      return polylineType;
    }

    /**
     * Indicates whether duplicate edges in the input should be kept (KEEP) or merged together
     * (MERGE). Note you can use edge labels to determine which input edges were merged into a given
     * output edge.
     *
     * <p>The default is DuplicateEdges.KEEP.
     */
    @CanIgnoreReturnValue
    public Options setDuplicateEdges(DuplicateEdges duplicateEdges) {
      this.duplicateEdges = duplicateEdges;
      return this;
    }

    /** Returns if duplicate edges should be kept or merged. */
    public DuplicateEdges duplicateEdges() {
      return duplicateEdges;
    }

    /**
     * Indicates whether sibling edge pairs (i.e., pairs consisting of an edge and its reverse edge)
     * should be kept (KEEP) or discarded (DISCARD). For example, if a polyline backtracks on
     * itself, the DISCARD option would cause this section of the polyline to be removed. Note that
     * this option may cause a single polyline to split into several pieces (e.g., if a polyline has
     * a "lollipop" shape).
     *
     * <p>siblingPairs must be either DISCARD or KEEP. The CREATE and REQUIRE options are not
     * allowed. The default is KEEP.
     */
    @CanIgnoreReturnValue
    public Options setSiblingPairs(SiblingPairs siblingPairs) {
      Preconditions.checkArgument(
          siblingPairs == SiblingPairs.KEEP || siblingPairs == SiblingPairs.DISCARD);
      this.siblingPairs = siblingPairs;
      return this;
    }

    /** Returns if sibling pairs should be kept or discarded. */
    public SiblingPairs siblingPairs() {
      return siblingPairs;
    }

    /**
     * If true, calls findValidationError() on the output polylines. If any error is found, it will
     * be returned by {@link S2Builder#build(S2Error)}. The default is false.
     */
    @CanIgnoreReturnValue
    public Options setValidate(boolean validate) {
      this.validate = validate;
      return this;
    }

    /** Returns true if polylines are to be validated. */
    public boolean validate() {
      return validate;
    }
  }
}
