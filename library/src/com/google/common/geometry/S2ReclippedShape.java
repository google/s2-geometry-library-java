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

import com.google.common.geometry.S2ContainsPointQuery.S2VertexModel;
import com.google.common.geometry.S2IndexCellData.EdgeAndIdChain;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.List;
import jsinterop.annotations.JsConstructor;

/**
 * A class representing a portion of an S2ClippedShape that has been re-clipped into a new cell.
 * Each edge is stored explicitly as two S2Points along with its chain and chain offset, to help
 * reduce the cost of edge id lookup.
 *
 * <p>Initializing an S2ReclippedShape requires an S2RobustCellClipper configured for the cell to
 * re-clip to, which must be contained by the S2ClippedShape's cell (equality is allowed).
 *
 * <p>If the S2RobustCellClipper is configured to actually compute crossings (rather than merely
 * test for edge/cell intersection), then the crossings are also stored and accessible through the
 * crossings() method.
 *
 * <p>This class stores more state than we otherwise would if this was meant to be e.g. a
 * serializable format. But this is intended to make join and filter processing simpler, so some
 * redundancy is useful for convenience.
 */
@SuppressWarnings("Assertion")
public class S2ReclippedShape {

  private final List<S2RobustCellClipper.Crossing> crossings = new ArrayList<>();
  private final List<Edge> edges = new ArrayList<>();

  private S2CellId cellId;
  private int dimension = -1;
  private int shapeId = -1;
  private boolean containsCenter;

  /**
   * Resets the shape id of the associated shape. The next call to init() after calling will always
   * clip edges.
   */
  public void reset() {
    shapeId = -1;
  }

  /**
   * Reclip a given clipped shape, using data from the given S2IndexCellData instance. The clipper
   * must be set to the cell to clip to, which must be contained by the index cell (equality is
   * allowed).
   *
   * <p>If the clipped shape id matches shapeId(), then it's presumed the shape has already been
   * clipped and reprocessing the edges is skipped. This allows calling init() in the inner loop
   * during join processing without duplicating work. To force processing of edges, call reset()
   * first.
   *
   * <p>When the S2RobustCellClipper is configured to compute crossings, then the crossings are
   * saved and accessible through the crossings() method. This behavior may be overridden by setting
   * save_crossings to false. The crossings are only valid until the clipper is reset or the clipper
   * is set to a new cell.
   *
   * <p>Returns false if processing was skipped, true otherwise.
   */
  @CanIgnoreReturnValue
  public boolean init(
      S2RobustCellClipper clipper,
      S2ClippedShape clipped,
      S2IndexCellData indexCellData,
      boolean saveCrossings) {
    assert indexCellData.id().contains(clipper.cell().id());

    if (shapeId == clipped.shapeId()) {
      return false;
    }

    shapeId = clipped.shapeId();
    dimension = indexCellData.shape(shapeId).dimension();

    S2Cell clipCell = clipper.cell();
    cellId = clipCell.id();

    edges.clear();
    crossings.clear();

    // Recompute the center containment bit if we're clipping to a smaller cell.
    containsCenter = clipped.containsCenter();
    if (dimension == 2 && !clipCell.id().equals(indexCellData.id())) {
      containsCenter =
          indexCellData.shapeContains(clipped, clipCell.getCenter(), S2VertexModel.SEMI_OPEN);
    }

    // Clip edges to the cell, save any that intersect it.
    clipper.reset(); // Invalidates any existing crossings.
    for (EdgeAndIdChain edge : indexCellData.shapeEdges(shapeId)) {
      S2RobustCellClipper.RobustClipResult result = clipper.clipEdge(edge);
      if (result.hit()) {
        Edge newEdge = new Edge(edge);
        edges.add(newEdge);
        newEdge.v0Contained = result.v0Inside();
        newEdge.v1Contained = result.v1Inside();
      }
    }

    // Optionally save the crossings that were generated.
    if (saveCrossings && clipper.options().enableCrossings()) {
      // clipper.getCrossings() returns crossings that will be reused after clipper.reset() is
      // called.
      for (S2RobustCellClipper.Crossing crossing : clipper.getCrossings()) {
        crossings.add(crossing);
      }
    }

    return true;
  }

  /** As above, but defaults saveCrossings to true. */
  @CanIgnoreReturnValue
  public boolean init(S2RobustCellClipper clipper, S2ClippedShape clipped, S2IndexCellData data) {
    return init(clipper, clipped, data, true);
  }

  /** Returns the S2CellId for the cell that the shape was clipped to. */
  public S2CellId cellId() {
    return cellId;
  }

  /** Returns the shape id for the shape that was clipped. */
  public int shapeId() {
    return shapeId;
  }

  /**
   * Returns the dimension of the shape for convenience. Returns -1 if no shape has been clipped.
   */
  public int dimension() {
    if (shapeId() >= 0) {
      return dimension;
    }
    return -1;
  }

  /** Returns true if the reclipped shape contains the cell center. */
  public boolean containsCenter() {
    return containsCenter;
  }

  /** Returns the list of edges. */
  public List<Edge> edges() {
    return edges;
  }

  /**
   * Returns the list of crossings (if enabled) that were computed. Crossings are reused when the
   * clipper provided to the init() method is reset or set to a new cell.
   */
  public List<S2RobustCellClipper.Crossing> crossings() {
    return crossings;
  }

  /** Checks that the given point is contained by the bounds of the clipped cell. */
  private boolean pointInCell(S2Point point) {
    R2Rect bounds = new S2PaddedCell(cellId(), S2ShapeIndex.CELL_PADDING).bound();

    R2Vector uv = new R2Vector();
    S2Cell cell = new S2Cell(cellId());
    S2Projections.validFaceXyzToUv(cell.face(), point, uv);
    return bounds.contains(uv);
  }

  /**
   * Returns true if the reclipped shape contains the given point under the given vertex model. This
   * is a faster alternative to a full S2ContainsPointQuery as no additional iterator lookup or cell
   * decoding is needed.
   *
   * <p>The tradeoff for this performance is that only points within the clipped cell may be tested.
   * This is enforced via an assertion. When assertions are disabled, points outside the cell may
   * return incorrect results.
   *
   * <p>Note that we don't store the cell center in this class, so it must be passed in.
   */
  public boolean contains(S2Point center, S2Point point, S2VertexModel model) {
    assert pointInCell(point);

    assert center.equalsPoint(cellId().toPoint());
    assert shapeId() >= 0;

    // Points and polylines don't contain anything except under the CLOSED model.
    if (dimension() < 2) {
      if (model != S2VertexModel.CLOSED) {
        return false;
      }

      // Point/polyline only contains point if it's a vertex.
      for (Edge edge : edges()) {
        if (edge.v0.equalsPoint(point) || edge.v1.equalsPoint(point)) {
          return true;
        }
      }
      return false;
    }

    // Test containment by drawing a line segment from the cell center to the given point and
    // counting edge crossings.
    S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(center, point);

    boolean inside = containsCenter();
    for (Edge edge : edges()) {
      int sign = crosser.robustCrossing(edge.v0, edge.v1);
      if (sign < 0) {
        continue;
      }
      if (sign == 0) {
        // For the OPEN and CLOSED models, check whether point is a vertex.
        if (model != S2VertexModel.SEMI_OPEN
            && (edge.v0.equalsPoint(point) || edge.v1.equalsPoint(point))) {
          return (model == S2VertexModel.CLOSED);
        }
        sign = S2EdgeUtil.vertexCrossing(crosser.a(), crosser.b(), edge.v0, edge.v1) ? 1 : 0;
      }
      inside ^= (sign != 0);
    }
    return inside;
  }

  /** As above, but defaults the vertex model to OPEN. */
  public boolean contains(S2Point center, S2Point point) {
    return contains(center, point, S2VertexModel.OPEN);
  }

  /**
   * A reclipped shape edge, as start and end points, edge chain, chain offset and a flag for each
   * vertex indicating whether it was contained in the cell.
   *
   * <p>Note that equality and hashing consider only the start and end points. This is intentional:
   * we want to consider two edges from different shapes as equal if they have the same start and
   * end points, even though they will have different chainId and offset values.
   */
  public static class Edge {

    // TODO(user): Create or reuse an EdgeList, or just have reclipped results in new
    // S2ClippedShapes. Not all of the data stored in these edges is used, and storing Edges in
    // an ArrayList is memory-inefficient. Also, many points are stored twice, once at each end of
    // an edge.
    public final S2Point v0;
    public final S2Point v1;
    public final int chainId;
    public final int offset;
    public boolean v0Contained;
    public boolean v1Contained;

    /** Constructor from EdgeAndIdChain. Sets all fields except v0Contained and v1Contained. */
    public Edge(EdgeAndIdChain edge) {
      this(edge.start(), edge.end(), edge.chainId(), edge.offset());
    }

    /** Constructor that sets only the edge endpoints. */
    public Edge(S2Point v0, S2Point v1) {
      this(v0, v1, -1, -1);
    }

    /** Constructor that sets all fields except v0Contained and v1Contained. */
    @JsConstructor
    public Edge(S2Point v0, S2Point v1, int chainId, int offset) {
      this.v0 = v0;
      this.v1 = v1;
      this.chainId = chainId;
      this.offset = offset;
    }

    /** Edges are equal if their endpoints are equal. Other fields are ignored. */
    @Override
    public boolean equals(Object o) {
      if (!(o instanceof Edge)) {
        return false;
      }
      Edge other = (Edge) o;
      return v0.equalsPoint(other.v0) && v1.equalsPoint(other.v1);
    }

    @Override
    public int hashCode() {
      return v0.hashCode() * 13 + v1.hashCode();
    }

    /** Creates an Edge that is the reverse of this edge. Only sets the start and end points. */
    public Edge reversed() {
      return new Edge(v1, v0);
    }

    /** Returns true if this edge starts at the given point, i.e. is outgoing from the point. */
    public boolean outgoing(S2Point point) {
      return v0.equalsPoint(point);
    }

    /** Returns true if this edge ends at the given point, i.e. is incoming to the point. */
    public boolean incoming(S2Point point) {
      return v1.equalsPoint(point);
    }

    /**
     * Returns true if this edge is incident on the given point, i.e. the point is either one of its
     * endpoints.
     */
    public boolean incidentOn(S2Point point) {
      return v0.equalsPoint(point) || v1.equalsPoint(point);
    }

    /** Returns true if this edge has both endpoints exactly equal. */
    public boolean isDegenerate() {
      return v0.equalsPoint(v1);
    }

    /** Returns true if this edge has the same endpoints as the given 'other' edge. */
    public boolean hasSameEndpoints(Edge other) {
      return v0.equalsPoint(other.v0) && v1.equalsPoint(other.v1);
    }

    /** Returns true if this edge has the opposite endpoints as the given 'other' edge. */
    public boolean isSiblingOf(Edge other) {
      return v0.equalsPoint(other.v1) && v1.equalsPoint(other.v0);
    }
  }
}
