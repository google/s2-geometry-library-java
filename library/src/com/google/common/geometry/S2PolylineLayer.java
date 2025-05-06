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
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2BuilderGraph.LabelFetcher;
import com.google.common.geometry.S2BuilderGraph.PolylineBuilder;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.List;

/**
 * A layer type that assembles edges (directed or undirected) into a single S2Polyline. Returns an
 * error if the edges cannot be assembled into a single unbroken polyline.
 *
 * <p>Duplicate edges are handled correctly (e.g., if a polyline backtracks on itself, or loops
 * around and retraces some of its previous edges.) The implementation attempts to preserve the
 * order of input edges whenever possible, so that if the input is a polyline and it is not modified
 * by S2Builder, then the output will be the same polyline (even if the polyline backtracks on
 * itself or forms a loop).
 *
 * <p>Note that by default, S2PolylineLayer does not validate the output polyline. To enable
 * validation, use {@link Options#setValidate(boolean)}, which will cause {@link
 * #build(S2BuilderGraph, S2Error)} to return false and set the error if the output polyline is
 * invalid.
 *
 * <p>S2PolylineLayer does not support options such as discarding sibling pairs or merging duplicate
 * edges because these options can split the polyline into several pieces. Use S2PolylineVectorLayer
 * if you need these features.
 */
public class S2PolylineLayer implements S2BuilderShapesLayer {
  /** The options for building the S2Polyline. */
  private final Options options;

  /** The new polyline constructed during each call to build(). */
  private S2Polyline polyline = null;

  /**
   * If provided to the constructor, labelsetMaps shape edge ids to label set ids after build().
   * Otherwise null.
   */
  private final IntVector labelSetIds;

  /**
   * If provided to the constructor, labelSetLexicon provides label sets from label set ids after
   * build(). Otherwise null.
   */
  private final IdSetLexicon labelSetLexicon;

  /** The polyline builder is reused for each build() call. */
  private final PolylineBuilder polylineBuilder = new PolylineBuilder();

  /** The label fetcher is reused for each build() call, but constructed only if needed. */
  private LabelFetcher fetcher = null;

  /** Constructs an S2PolylineLayer using default {@link Options}. */
  public S2PolylineLayer() {
    this.options = new Options();
    labelSetIds = null;
    labelSetLexicon = null;
  }

  /** Constructs an S2PolylineLayer with the provided options. */
  public S2PolylineLayer(Options options) {
    Preconditions.checkNotNull(options);
    this.options = options;
    labelSetIds = null;
    labelSetLexicon = null;
  }

  /**
   * Constructs an S2PolylineLayer with the provided options. Each time build() is called, the
   * provided IdSetLexicon and IntVector will be filled in so that labels attached to the input
   * edges can be obtained for the edges of the output polygon that the input edges were snapped to.
   */
  public S2PolylineLayer(Options options, IdSetLexicon labelSetLexicon, IntVector labelSetIds) {
    Preconditions.checkNotNull(options);
    Preconditions.checkNotNull(labelSetLexicon);
    Preconditions.checkNotNull(labelSetIds);
    this.options = options;
    this.labelSetIds = labelSetIds;
    this.labelSetLexicon = labelSetLexicon;
  }

  /** Returns a copy of the current options in a new Options object. */
  public Options options() {
    return new Options(options);
  }

  /** Returns the new output polyline after build() has been called, if build() returned true. */
  public S2Polyline getPolyline() {
    Preconditions.checkNotNull(polyline, "Must call build() first.");
    return polyline;
  }

  @Override
  public String toString() {
    return "S2PolylineLayer with options " + options;
  }

  @Override
  public Iterable<S2Polyline> shapes() {
    return ImmutableList.of(getPolyline());
  }

  @Override
  public GraphOptions graphOptions() {
    // Remove edges that collapse to a single vertex, but keep duplicate and sibling edges, since
    // merging duplicates or discarding siblings can make it impossible to assemble the edges into a
    // single polyline.
    return new GraphOptions(
        options.edgeType(),
        GraphOptions.DegenerateEdges.DISCARD,
        GraphOptions.DuplicateEdges.KEEP,
        GraphOptions.SiblingPairs.KEEP);
  }

  /**
   * Returns true if building an S2Polyline succeeded. Otherwise returns false and fills in 'error'.
   */
  @Override
  public boolean build(S2BuilderGraph g, S2Error error) {
    if (labelSetIds != null) {
      labelSetIds.clear();
      labelSetLexicon.clear();
    }

    if (g.numEdges() == 0) {
      polyline = new S2Polyline(ImmutableList.of());
      return true;
    }

    polylineBuilder.init(g);
    List<int[]> edgePolylines = polylineBuilder.buildWalks();
    if (edgePolylines.size() != 1) {
      error.init(
          S2Error.Code.BUILDER_EDGES_DO_NOT_FORM_POLYLINE,
          "Input edges cannot be assembled into polyline");
      return false;
    }

    int[] edgePolyline = edgePolylines.get(0);
    S2Point[] vertices = new S2Point[edgePolyline.length + 1];
    vertices[0] = g.vertex(g.edgeSrcId(edgePolyline[0]));
    int i = 1;
    for (int edgeId : edgePolyline) {
      vertices[i++] = g.vertex(g.edgeDstId(edgeId));
    }

    if (labelSetIds != null) {
      if (fetcher == null) {
        fetcher = new LabelFetcher(g, options.edgeType());
      } else {
        fetcher.init(g, options.edgeType());
      }
      IntVector labels = new IntVector(); // Temporary storage for labels.
      for (int edgeId : edgePolyline) {
        fetcher.fetch(edgeId, labels);
        labelSetIds.add(labelSetLexicon.add(labels));
      }
    }

    polyline = new S2Polyline(vertices);
    if (options.validate()) {
      return !polyline.findValidationError(error);
    }
    return true;
  }

  /** Options for the S2PolylineLayer. */
  public static class Options {
    private S2Builder.EdgeType edgeType;
    private boolean validate;

    /** Options constructor that sets default options: directed edges and no validation. */
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

    /** Options copy constructor. */
    public Options(Options other) {
      edgeType = other.edgeType;
      validate = other.validate;
    }

    @Override
    public String toString() {
      return Strings.lenientFormat("EdgeType %s, Validate=%s", edgeType, validate);
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
     * Sets whether output polylines should be validated. If true, calls findValidationError() on
     * the output polylines. Any error found will be returned by {@link S2Builder#build(S2Error)}.
     * The default is false.
     */
    @CanIgnoreReturnValue
    public Options setValidate(boolean validate) {
      this.validate = validate;
      return this;
    }

    /** Returns if polylines should be validated. */
    public boolean validate() {
      return validate;
    }
  }
}
