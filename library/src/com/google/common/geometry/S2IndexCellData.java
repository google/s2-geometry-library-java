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

import static java.lang.Math.min;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2ContainsPointQuery.S2VertexModel;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2Shape.ChainPosition;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import com.google.common.geometry.primitives.IntVector;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.List;

/**
 * A class for working with {@link S2ShapeIndex.Cell} data. For larger queries like validation, we
 * often look up edges multiple times, and sometimes need to work with the edges themselves, their
 * edge ids, or their chain and offset. This class avoids the costs of repeated calls to virtual
 * methods, S2Shape decoding that may not be aligned with the required edges, and repeated chain
 * offset computation. The access methods here do no work, so they inline, which is essential for
 * looping over cell data multiple times.
 *
 * <p>This is meant to support larger querying and validation operations such as S2ValidationQuery
 * that have to proceed cell-by cell through an index.
 *
 * <p>To use, simply call {@link #loadCell(S2ShapeIndex, S2CellId, S2ShapeIndex.Cell)} to decode the
 * contents of a cell. Or, call {@link #loadCell(S2ShapeIndex, S2CellId, S2ShapeIndex.Cell, int)} to
 * decode the edges of a single clipped shape in the cell.
 *
 * <p>The class promises that the edges will be looked up once when loadCell() or loadShape() is
 * called, and the edges, edgeIds, chain, and chain offsets are loaded into contiguous arrays that
 * we can serve requests from efficiently. Since the chain and offset are computed anyways when
 * looking up an edge via the shape.edge() API, we simply cache those values so the cost is minimal.
 *
 * <pre>
 * The order of edges looks like this:
 *
 *   |     0D Shapes     |     1D Shapes     |     2D Shapes     |  Dimensions
 *   |  5  |   1   |  3  |  2  |   7   |  0  |  6  |   4   |  8  |  Shapes
 *   [ ......................... Edges ..........................]  Edges
 *
 * This allows us to look up individual shapes very quickly, as well as all shapes in a given
 * dimension or contiguous range of dimensions:
 *
 *   edges()        - Return list of all edges.
 *   shapeEdges()   - Return sublist of edges of a given shape.
 *   dimEdges()     - Return sublist of all edges of a given dimension.
 *   dimRangeEges() - Return sublist of all edges of a range of dimensions.
 *
 * </pre>
 *
 * We use a stable sort, so similarly to S2ShapeIndex.Cell, we promise that shapes _within a
 * dimension_ are in the same order they are in the index itself, and the edges _within a shape_ are
 * similarly in the same order.
 *
 * <p>The clipped shapes in a cell are exposed through the shapes() method.
 */
@SuppressWarnings("Assertion")
public class S2IndexCellData {

  /** A visitor for shape edges with their endpoints, edge ids, chain id, and offset. */
  public interface Visitor {
    boolean visit(int edgeId, int chainId, int offset, S2Point start, S2Point end);
  }

  /** An interface to a shape edge with its endpoints, edge id, chain id, and offset. */
  public interface EdgeAndIdChain {
    /** Returns the start vertex of the edge. */
    public S2Point start();

    /** Returns the end vertex of the edge. */
    public S2Point end();

    /** Returns the id of the edge within its shape. */
    public int edgeId();

    /** Returns the id of the chain the edge belongs to. */
    public int chainId();

    /** Returns the offset of the edge within the chain. */
    public int offset();

    /**
     * Two EdgeAndChainId instances are considered equal if they have the same endpoint vertices in
     * the same order.
     */
    default boolean isEqualTo(EdgeAndIdChain other) {
      return start().equalsPoint(other.start()) && end().equalsPoint(other.end());
    }

    /** Returns true if this edge has the opposite endpoints as the other edge. */
    default boolean isReverseOf(EdgeAndIdChain other) {
      return start().equalsPoint(other.end()) && end().equalsPoint(other.start());
    }

    /**
     * Two EdgeAndChainId instances are considered equal if they have the same endpoint vertices.
     */
    public static boolean equals(EdgeAndIdChain x, EdgeAndIdChain y) {
      return x.isEqualTo(y);
    }

    /**
     * EdgeAndIdChain instances are compared by their starting vertices, breaking ties with their
     * end vertices.
     */
    public static boolean lessThan(EdgeAndIdChain x, EdgeAndIdChain y) {
      return x.start().compareTo(y.start()) < 0
          || (x.start().equalsPoint(y.start()) && x.end().compareTo(y.end()) < 0);
    }

    /** Returns true if the given point is an endpoint of this edge. */
    default boolean hasEndpoint(S2Point point) {
      return start().equalsPoint(point) || end().equalsPoint(point);
    }
  }

  /** Simple pair for defining an integer value range. A Range may be empty, with size zero. */
  private static class Range {
    int start = 0;
    int size = 0;

    /** Constructs an empty Range. */
    public Range() {}

    /** Constructs a Range with the given start and size. */
    public Range(int start, int size) {
      Preconditions.checkArgument(start >= 0);
      Preconditions.checkArgument(size >= 0);
      this.start = start;
      this.size = size;
    }

    /** Clears the range. */
    public void clear() {
      start = 0;
      size = 0;
    }
  }

  /** Essentially a Pair for mapping a shape id to a Range. */
  private static class ShapeRange {
    public int shapeId;
    public Range range;

    public static ShapeRange of(int shapeId, Range range) {
      ShapeRange shapeRange = new ShapeRange();
      shapeRange.shapeId = shapeId;
      shapeRange.range = range;
      return shapeRange;
    }
  }

  /** The shape index that the currently loaded cell belongs to. */
  private S2ShapeIndex index = null;

  /** The currently loaded cell. */
  private S2ShapeIndex.Cell cell = null;

  /** The id of the currently loaded cell. */
  private S2CellId cellId;

  /** The id of the currently loaded shape, or -1 if all shapes in the cell are loaded. */
  private int loadedShapeId;

  // Computing the cell center and S2Cell can cost as much as looking up the edges themselves, so
  // defer doing it until needed.
  private volatile boolean s2cellSet = false;
  private volatile boolean centerSet = false;
  private S2Cell s2cell;
  private S2Point cellCenter;

  /** Dimensions that we wish to decode. The default is all of them. */
  private final boolean[] dimWanted = {true, true, true};

  // Parallel resizable vectors storing data about the edges of the current cell.

  /** Ids of the edges within their shapes. */
  private final IntVector edgeIds = new IntVector();

  /** Ids of the shape chains the edges belongs to. */
  private final IntVector chainIds = new IntVector();

  /** Offset of the edges within their shape chains. */
  private final IntVector offsets = new IntVector();

  /** The source vertices of the edges. */
  private final ArrayList<S2Point> starts = new ArrayList<>();

  /** The destination vertices of the edges. */
  private final ArrayList<S2Point> ends = new ArrayList<>();

  /** Maps each shape id to the range of the parallel arrays that shape's edges are stored in. */
  private final List<ShapeRange> shapeRanges = new ArrayList<>();

  /** Maps each dimension to the range of the parallel arrays those edges are stored in. */
  private final Range[] dimRanges = new Range[3];

  /** Constructs an empty S2IndexCellData. A cell must be loaded via loadCell() before use. */
  public S2IndexCellData() {
    dimRanges[0] = new Range();
    dimRanges[1] = new Range();
    dimRanges[2] = new Range();
  }

  /** Constructs an S2IndexCellData from the given index, cell id, and cell. */
  public S2IndexCellData(S2ShapeIndex index, S2CellId id, S2ShapeIndex.Cell cell) {
    loadCell(index, id, cell);
  }

  /**
   * Resets internal state to defaults without shrinking the parallel arrays. The next call to
   * loadCell() or loadShape() will process the cell regardless of whether it was already loaded.
   * Should be called when processing a new index.
   */
  public void reset() {
    index = null;
    cell = null;

    edgeIds.clear();
    chainIds.clear();
    offsets.clear();
    starts.clear();
    ends.clear();

    shapeRanges.clear();
    dimWanted[0] = true;
    dimWanted[1] = true;
    dimWanted[2] = true;
  }

  /** Returns true if a dimension (0, 1, or 2) is set to be decoded. */
  public boolean dimWanted(int dim) {
    Preconditions.checkArgument(0 <= dim && dim <= 2);
    return dimWanted[dim];
  }

  /**
   * Configures whether a particular dimension of shape should be decoded.
   *
   * <p>TODO(torrey): When additional dimensions are added, loadCell() should reload the same cell
   * if needed. Also in C++.
   */
  public void setDimWanted(int dim, boolean wanted) {
    Preconditions.checkArgument(0 <= dim && dim <= 2);
    dimWanted[dim] = wanted;
  }

  /** Returns the id of the currently loaded cell. */
  public S2CellId id() {
    return cellId;
  }

  /** Returns the index of the currently loaded cell. */
  public S2ShapeIndex index() {
    return index;
  }

  /** Returns the id of the currently loaded shape, or -1 if all shapes in the cell are loaded. */
  public int loadedShapeId() {
    return loadedShapeId;
  }

  /** Returns an S2Cell instance for the currently loaded cell. */
  public S2Cell cell() {
    // TODO(torrey): Consider a mutex or similar here, as C++ uses, to make this thread safe.
    if (!s2cellSet) {
      s2cell = new S2Cell(cellId);
      s2cellSet = true;
    }
    return s2cell;
  }

  /** Returns the center point of the currently loaded cell. */
  public S2Point center() {
    // TODO(torrey): Consider a mutex or similar here, as C++ uses, to make this thread safe.
    if (!centerSet) {
      cellCenter = cellId.toPoint();
      centerSet = true;
    }
    return cellCenter;
  }

  /**
   * Loads the data from the given cell. Previous cell data is cleared. If the index, id and cell
   * pointer are the same as in the previous call to loadCell, loading is not performed since we
   * already have the data decoded.
   */
  public void loadCell(S2ShapeIndex index, S2CellId id, S2ShapeIndex.Cell cell) {
    loadCell(index, id, cell, -1);
  }

  /**
   * Loads the data from the given cell. Previous cell data is cleared. If the index, S2CellId, and
   * loaded shapeId are the same as in the previous call to loadCell, loading is not performed since
   * we already have the data decoded. Only shapes with the required dimension are loaded, all
   * others are ignored. If shapeId is -1, all shapes in the cell are loaded, otherwise only the
   * shape with the specified id is loaded.
   */
  public void loadCell(S2ShapeIndex index, S2CellId id, S2ShapeIndex.Cell cell, int shapeId) {
    assert index != null;
    if (this.index == index && this.cellId.equals(id) && this.loadedShapeId == shapeId) {
      return;
    }
    this.index = index;

    // Cache cell information.
    this.loadedShapeId = shapeId;
    this.cell = cell;
    this.cellId = id;

    // Reset flags so we'll recompute cached values.
    // TODO(torrey): Consider a mutex or similar here, as C++ uses, to make this thread safe.
    s2cellSet = false;
    centerSet = false;

    // Clear previous edges.
    starts.clear();
    ends.clear();
    edgeIds.clear();
    chainIds.clear();
    offsets.clear();

    shapeRanges.clear();

    // Reset per-dimension range information.
    for (Range range : dimRanges) {
      range.clear();
    }

    int minDim = 0;
    while (minDim <= 2 && !dimWanted(minDim)) {
      ++minDim;
    }

    int maxDim = 2;
    while (maxDim >= 0 && !dimWanted(maxDim)) {
      --maxDim;
    }

    // No dimensions wanted, we're done.
    if (minDim > 2 || maxDim < 0) {
      return;
    }

    for (int dim = minDim; dim <= maxDim; ++dim) {
      int dimStart = edgeIds.size();

      for (int clippedShapeIndex = 0; clippedShapeIndex < cell.numShapes(); ++clippedShapeIndex) {
        S2ClippedShape clipped = cell.clipped(clippedShapeIndex);
        int clippedShapeId = clipped.shapeId();
        // If we're only loading a single shape, skip any others.
        if (loadedShapeId != -1 && clippedShapeId != loadedShapeId) {
          continue;
        }
        S2Shape shape = index.getShapes().get(clippedShapeId);

        // Only process the current dimension.
        if (shape.dimension() != dim) {
          continue;
        }

        // In the event we wanted dimensions 0 and 2, but not 1.
        if (!dimWanted(dim)) {
          continue;
        }

        // Materialize clipped shape edges into the edges list. Track where we start so we can add
        // information about the range for this shape.
        int shapeStart = edgeIds.size();
        // Reusable ChainPosition and MutableEdge instances.
        ChainPosition position = new ChainPosition();
        MutableEdge edge = new MutableEdge();
        for (int i = 0; i < clipped.numEdges(); ++i) {
          int edgeId = clipped.edge(i);

          // Looking up an edge requires looking up which chain it's in, which is often a binary
          // search. So let's manually lookup the chain information and use that to find the edge,
          // so we only have to do that search once.
          shape.getChainPosition(edgeId, position);
          shape.getChainEdge(position.chainId, position.offset, edge);

          starts.add(edge.a);
          ends.add(edge.b);
          edgeIds.add(edgeId);
          chainIds.add(position.chainId);
          offsets.add(position.offset);
        }

        // Record the range of edges belonging to the shape.
        shapeRanges.add(
            ShapeRange.of(clippedShapeId, new Range(shapeStart, edgeIds.size() - shapeStart)));
      }

      // Record the range of edges belonging to the current dimension.
      dimRanges[dim] = new Range(dimStart, edgeIds.size() - dimStart);
    }
  }

  /** Returns the number of clipped shapes in the currently loaded cell. */
  public int numClipped() {
    return cell.numShapes();
  }

  /** Returns the S2Shape for the given clipped shape. */
  public S2Shape shape(S2ClippedShape clipped) {
    S2Shape shape = index.getShapes().get(clipped.shapeId());
    assert shape != null;
    return shape;
  }

  /** Returns the S2Shape for the given shape id from the index. */
  public S2Shape shape(int shapeId) {
    S2Shape shape = index.getShapes().get(shapeId);
    assert shape != null;
    return shape;
  }

  /** Returns the clipped shapes in the current cell. */
  public List<S2ClippedShape> clippedShapes() {
    return cell.clippedShapes();
  }

  /**
   * Visits the subset of loaded edges defined by the given range, until the visitor returns false,
   * in which case false is returned, or all edges in the range have been visited, in which case
   * true is returned.
   */
  private boolean visitRange(int start, int size, Visitor visitor) {
    for (int i = start; i < start + size; ++i) {
      if (!visitor.visit(
          edgeIds.get(i), chainIds.get(i), offsets.get(i), starts.get(i), ends.get(i))) {
        return false;
      }
    }
    return true;
  }

  /** Returns an abstract list of the edges in the given range for the current cell. */
  private List<EdgeAndIdChain> edgeRange(int start, int size) {
    return new AbstractList<EdgeAndIdChain>() {

      @Override
      public EdgeAndIdChain get(int index) {
        return new EdgeAndIdChain() {
          @Override
          public S2Point start() {
            return starts.get(start + index);
          }

          @Override
          public S2Point end() {
            return ends.get(start + index);
          }

          @Override
          public int edgeId() {
            return edgeIds.get(start + index);
          }

          @Override
          public int chainId() {
            return chainIds.get(start + index);
          }

          @Override
          public int offset() {
            return offsets.get(start + index);
          }
        };
      }

      @Override
      public int size() {
        return size;
      }
    };
  }

  /** Returns an abstract list of all the edges in the current cell. */
  public List<EdgeAndIdChain> edges() {
    return edgeRange(0, edgeIds.size());
  }

  /**
   * Visits all the edges in the current cell until the visitor returns false, in which case false
   * is returned, or all edges have been visited, in which case true is returned.
   */
  @CanIgnoreReturnValue
  public boolean visitEdges(Visitor visitor) {
    return visitRange(0, edgeIds.size(), visitor);
  }

  /**
   * Returns an abstract list of the edges in the current cell for a given shape. This will be a
   * sub-list of the edges.
   *
   * <p>REQUIRES: id is an id for one of the shapes in shapes(). It may be that the shape has no
   * edges in the current cell, in which case an empty list is returned.
   */
  public List<EdgeAndIdChain> shapeEdges(int shapeId) {
    for (ShapeRange shapeRange : shapeRanges) {
      if (shapeRange.shapeId == shapeId) {
        Range range = shapeRange.range;
        if (range.start >= edgeIds.size()) {
          return ImmutableList.of();
        }
        return edgeRange(range.start, range.size);
      }
    }
    // This shape isn't present in the current cell.
    return ImmutableList.of();
  }

  /**
   * Visits all the edges in the current cell until the visitor returns false, in which case false
   * is returned, or all edges have been visited, in which case true is returned.
   *
   * <p>REQUIRES: id is an id for one of the shapes in shapes(). It may be that the shape has no
   * edges in the current cell, in which the visitor will not be called.
   */
  @CanIgnoreReturnValue
  public boolean visitShapeEdges(int shapeId, Visitor visitor) {
    for (ShapeRange shapeRange : shapeRanges) {
      if (shapeRange.shapeId == shapeId) {
        Range range = shapeRange.range;
        if (range.start >= edgeIds.size()) {
          return true;
        }
        return visitRange(range.start, range.size, visitor);
      }
    }
    // This shape isn't present in the current cell.
    return true;
  }

  /**
   * Returns the edges in the current cell for a given dimension.
   *
   * <p>REQUIRES: dim is a valid dimension (0, 1, or 2). It may be that there are no edges with the
   * given dimension in the current cell, in which case an empty list is returned.
   */
  public List<EdgeAndIdChain> dimEdges(int dim) {
    Preconditions.checkArgument(0 <= dim && dim <= 2);
    Range range = dimRanges[dim];
    if (range.start < edgeIds.size()) {
      return edgeRange(range.start, range.size);
    }
    return ImmutableList.of();
  }

  /**
   * Visits the edges in the current cell for a given dimension until the visitor returns false, in
   * which case false is returned, or all edges have been visited, in which case true is returned.
   *
   * <p>REQUIRES: dim is a valid dimension (0, 1, or 2). It may be that there are no edges with the
   * given dimension in the current cell, in which case the visitor will not be called.
   */
  @CanIgnoreReturnValue
  public boolean visitDimEdges(int dim, Visitor visitor) {
    Preconditions.checkArgument(0 <= dim && dim <= 2);
    Range range = dimRanges[dim];
    if (range.start < edgeIds.size()) {
      return visitRange(range.start, range.size, visitor);
    }
    return true;
  }

  /**
   * Returns the edges in the current cell for a range of dimensions.
   *
   * <p>REQUIRES: dimensions are valid (0, 1, or 2)
   */
  public List<EdgeAndIdChain> dimRangeEdges(int dim0, int dim1) {
    Preconditions.checkArgument(dim0 <= dim1);
    Preconditions.checkArgument(0 <= dim0 && dim0 <= 2);
    Preconditions.checkArgument(0 <= dim1 && dim1 <= 2);

    int start = dimRanges[dim0].start;
    int size = 0;
    for (int dim = dim0; dim <= dim1; ++dim) {
      start = min(start, dimRanges[dim].start);
      size += dimRanges[dim].size;
    }

    if (start < edgeIds.size()) {
      return edgeRange(start, size);
    }
    return ImmutableList.of();
  }

  /**
   * Tests whether a shape in the current cell contains a point. The logic here is is the same as
   * {@link S2ContainsPointQuery#shapeContains(int, S2Point)} but only applies to the current cell
   * and doesn't have to lookup edges or cell centers again.
   *
   * <p>Since this _only_ operates on edges within the current cell, it must not be used to test
   * points that are outside of the current cell, as there may be intervening edges between the
   * point and the cell center that we can't see.
   *
   * <p>REQUIRES: The cell must contain the given point.
   */
  public boolean shapeContains(S2ClippedShape clipped, S2Point point, S2VertexModel model) {
    // When assertions are enabled, verify that the cell contains the point.
    assert new S2PaddedCell(id(), S2ShapeIndex.CELL_PADDING)
        .bound()
        .contains(S2Projections.validFaceXyzToUv(id().face(), point));

    // Points and polylines don't contain anything except under the CLOSED model.
    S2Shape shape = index.getShapes().get(clipped.shapeId());
    if (shape.dimension() < 2) {
      if (model != S2VertexModel.CLOSED) {
        return false;
      }

      // Point/polyline only contains point if it's a vertex.
      for (EdgeAndIdChain edge : shapeEdges(clipped.shapeId())) {
        // TODO(torrey): For consecutive edge chains, this checks each point twice. Also in C++.
        if (edge.hasEndpoint(point)) {
          return true;
        }
      }
      return false;
    }

    // Test containment by drawing a line segment from the cell center to the given point and
    // counting edge crossings.
    S2Point center = this.center();
    EdgeCrosser crosser = new EdgeCrosser(center, point);

    boolean[] inside = {clipped.containsCenter()};
    // for (EdgeAndIdChain edge : shapeEdges(clipped.shapeId())) {
    visitShapeEdges(
        clipped.shapeId(),
        (edgeId, chain, offset, start, end) -> {
          int sign = crosser.robustCrossing(start, end);
          if (sign < 0) {
            return true;
          }
          boolean crossing = (sign > 0);
          if (!crossing) {
            // No edge crossing, but the edge endpoint is equal to either the cell center or the
            // point. For the OPEN and CLOSED models, determine which of those it is. If it's the
            // point, then the point is a vertex of the shape, and we're done.
            if (model != S2VertexModel.SEMI_OPEN
                && (start.equalsPoint(point) || end.equalsPoint(point))) {
              return (model == S2VertexModel.CLOSED);
            }
            // Otherwise, for the SEMI_OPEN model, use vertexCrossing() to decide if it crosses or
            // not, and keep counting crossings.
            crossing = S2EdgeUtil.vertexCrossing(crosser.a(), crosser.b(), start, end);
          }
          inside[0] ^= crossing;
          return true;
        });

    return inside[0];
  }

  /**
   * As {@link shapeContains(S2ClippedShape, S2Point, S2VertexModel)} above, tests whether a shape
   * in the current cell contains a point. This version uses the OPEN vertex model.
   */
  public boolean shapeContains(S2ClippedShape clipped, S2Point point) {
    return shapeContains(clipped, point, S2VertexModel.OPEN);
  }
}
