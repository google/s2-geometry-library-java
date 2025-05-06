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
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2Builder.EdgeType;
import com.google.common.geometry.S2Builder.GraphOptions;
import com.google.common.geometry.S2Builder.GraphOptions.DegenerateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.DuplicateEdges;
import com.google.common.geometry.S2Builder.GraphOptions.SiblingPairs;
import com.google.common.geometry.S2BuilderGraph.LabelFetcher;
import com.google.common.geometry.primitives.IdSetLexicon;
import com.google.common.geometry.primitives.IntVector;
import java.util.ArrayList;
import java.util.List;

/**
 * A layer type that collects degenerate edges as points. This layer expects all edges to be
 * degenerate. In case of finding non-degenerate edges it sets S2Error but it still generates the
 * output with degenerate edges.
 */
@SuppressWarnings("Assertion")
public class S2PointVectorLayer implements S2BuilderShapesLayer {
  /** Options for building the point vector. */
  private final Options options;

  /** The new list of points constructed during each call to build(). */
  private List<S2Point> points = null;

  /**
   * If provided to the constructor, labelSetIds maps edgeId to label set ids after build().
   * Otherwise null.
   */
  private final IntVector labelSetIds;

  /**
   * If provided to the constructor, labelSetLexicon provides label sets from label set ids after
   * build(). Otherwise null.
   */
  private final IdSetLexicon labelSetLexicon;

  /** The label fetcher is reused for each build() call, but constructed only if needed. */
  private LabelFetcher fetcher = null;

  /**
   * Options for building the point vector. There is one option: should duplicate edges (i.e.
   * points) be merged or not. The default is to merge them.
   */
  public static class Options {
    private DuplicateEdges duplicateEdges;

    public Options() {
      duplicateEdges = DuplicateEdges.MERGE;
    }

    public Options(DuplicateEdges duplicateEdges) {
      this.duplicateEdges = duplicateEdges;
    }

    public void setDuplicateEdges(DuplicateEdges duplicateEdges) {
      this.duplicateEdges = duplicateEdges;
    }

    public DuplicateEdges duplicateEdges() {
      return duplicateEdges;
    }
  }

  /** Constructs a new S2PointVectorLayer with default options. */
  public S2PointVectorLayer() {
    labelSetIds = null;
    labelSetLexicon = null;
    this.options = new Options();
  }

  /** Constructs a new S2PointVectorLayer with the given options. */
  public S2PointVectorLayer(Options options) {
    labelSetIds = null;
    labelSetLexicon = null;
    this.options = options;
  }

  /**
   * Constructs a new S2PointVectorLayer with the given options, which will fill in the given
   * labelSetIds and labelSetLexicon during each call to build().
   */
  public S2PointVectorLayer(IntVector labelSetIds, IdSetLexicon labelSetLexicon, Options options) {
    this.labelSetIds = labelSetIds;
    this.labelSetLexicon = labelSetLexicon;
    this.options = options;
  }

  @Override
  public String toString() {
    return "S2PointVectorLayer with Options.duplicateEdges=" + options.duplicateEdges();
  }

  /** Returns the output points as a list after build() has been called. */
  public List<S2Point> getPointVector() {
    Preconditions.checkNotNull(points, "Must call build() first.");
    return points;
  }

  /** Returns the output points as an S2Point.Shape after build() has been called. */
  public S2Point.Shape getPointVectorShape() {
    Preconditions.checkNotNull(points, "Must call build() first.");
    return S2Point.Shape.fromList(points); // does not copy
  }

  @Override
  public Iterable<S2Point.Shape> shapes() {
    S2Point.Shape shape = getPointVectorShape();
    if (shape.numEdges() == 0) {
      return ImmutableList.of();
    }
    return ImmutableList.of(shape);
  }

  @Override
  public GraphOptions graphOptions() {
    return new GraphOptions(
        EdgeType.DIRECTED, DegenerateEdges.KEEP, options.duplicateEdges(), SiblingPairs.KEEP);
  }

  /**
   * Returns true if collecting a list of points (degenerate edges) from the given graph succeeded.
   * Otherwise returns false and fills in 'error'.
   */
  @Override
  public boolean build(S2BuilderGraph g, S2Error error) {
    if (labelSetIds != null) {
      labelSetIds.clear();
      labelSetLexicon.clear();

      if (fetcher == null) {
        fetcher = new LabelFetcher(g, EdgeType.DIRECTED);
      } else {
        fetcher.init(g, EdgeType.DIRECTED);
      }
    }

    // The output list of points.
    points = new ArrayList<>();

    // Temporary storage for the labels on one edge.
    IntVector labels = new IntVector();

    for (int edgeId = 0; edgeId < g.edges().size(); edgeId++) {
      int srcId = g.edgeSrcId(edgeId);
      int dstId = g.edgeDstId(edgeId);

      // If the edge endpoint vertex ids are different, we need to check the points.
      if (srcId != dstId) {
        S2Point src = g.vertex(srcId);
        S2Point dst = g.vertex(dstId);
        if (!src.equalsPoint(dst)) {
          error.init(S2Error.Code.INVALID_ARGUMENT, "Found non-degenerate edges");
          continue;
        }
      }

      // A degenerate edge representing a point.
      points.add(g.vertex(srcId));
      if (labelSetIds != null) {
        fetcher.fetch(edgeId, labels);
        labelSetIds.add(labelSetLexicon.add(labels));
      }
    }
    return error.ok();
  }
}
