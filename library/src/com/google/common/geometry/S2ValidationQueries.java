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

import static com.google.common.geometry.S2Predicates.orderedCCW;
import static java.lang.Math.max;

import com.google.common.base.Preconditions;
import com.google.common.geometry.S2ContainsPointQuery.S2VertexModel;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2IncidentEdgeTracker.IncidentEdgeKey;
import com.google.common.geometry.S2IncidentEdgeTracker.IncidentEdgeMap;
import com.google.common.geometry.S2IndexCellData.EdgeAndIdChain;
import com.google.common.geometry.S2Shape.ChainPosition;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import com.google.common.geometry.primitives.IntPairVector;
import com.google.common.geometry.primitives.IntVector;
import com.google.common.geometry.primitives.Ints.OfInt;
import com.google.common.geometry.primitives.PooledList;
import com.google.common.geometry.primitives.Sorter.SortableCollection;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

/**
 * S2ValidationQueries defines the following queries, which each have slightly different correctness
 * semantics for their particular domain. See the Javadoc on each class for details.
 *
 * <ul>
 *   <li>{@link S2ValidQuery}
 *   <li>{@link S2LegacyValidQuery}
 * </ul>
 */
@SuppressWarnings("Assertion")
public final class S2ValidationQueries {
  /**
   * Base class for validating geometry contained in an index according to a configurable model of
   * correctness. There are several different notions of "valid" geometry that we have to address,
   * including the basic requirements for S2BooleanOperation, S2Polygon specific checks, etc.
   *
   * <p>Classes of geometry can be thought of as sets of all shapes that have certain invariants
   * (such as interior-on-left or no-degenerate-edges). Classes then naturally build on other
   * classes by adding more rules to further restrict the allowed shapes.
   *
   * <p>This extends-upon structure naturally lends itself to an inheritance based implementation.
   * Beginning with the most permissive class of geometry which just meets the criteria for
   * S2BooleanOperation (S2ValidQuery), we can then build subsequent classes of geometry by
   * inheriting from and extending the checks that are performed.
   *
   * <p>Validation queries can overload the methods in the subclass API below. The validate() driver
   * function will call these functions in this order:
   *
   * {@snippet :
   * - start()
   *   - checkShape()
   *   - startCell()
   *     - startShape()
   *       - checkEdge()
   *     - finishShape()
   *  - finish()
   * }
   *
   * <p>Hoisting the loops into the base class like this allows us to fuse all the inner loops of
   * the various queries so that we only have to iterate over the index and its geometry once.
   *
   * <p>Several protected member functions are defined to give subclasses access to the data being
   * validated:
   *
   * {@snippet :
   *   // Returns the current index being operated on.
   *   S2ShapeIndex index();
   *
   *   // Returns the decoded data for the current cell.
   *   S2IndexCellData currentCell();
   *
   *   // Returns the current incident edge set. We promise that this set is updated with the
   *   // current cell's edges before startCell() is called.
   *   IncidentEdgeMap incidentEdges();
   * }
   *
   * <p>A query then has a {@code boolean validate(S2ShapeIndex index, S2Error error)} method which
   * validates the index and returns true if it's valid, otherwise false, with the validation
   * failure details provided through the error parameter.
   */
  public static class S2ValidationQueryBase {
    private final S2IndexCellData cellBuffer = new S2IndexCellData();
    private final S2IncidentEdgeTracker incidentEdgeTracker = new S2IncidentEdgeTracker();

    /** The index we're validating. */
    private S2ShapeIndex index = null;

    /**
     * Validate the given index by calling the hooks in the derived class. If an error is found, the
     * returned error is set. Otherwise, the returned error will have error.ok() == true.
     */
    public S2Error validate(S2ShapeIndex index) {
      S2Error error = new S2Error();
      validate(index, error);
      return error;
    }

    /**
     * Validate the given index by calling the hooks in the derived class. Returns false if any
     * error is found, and sets the error parameter. Otherwise, returns true.
     */
    @CanIgnoreReturnValue
    public boolean validate(S2ShapeIndex index, S2Error error) {
      this.index = index;
      boolean result = validate(error);
      this.index = null;
      return result;
    }

    /**
     * Validate the index by calling the hooks in the derived class. Returns false if any error is
     * found, and sets the error parameter. Otherwise, returns true.
     */
    private boolean validate(S2Error error) {
      incidentEdgeTracker.reset();
      cellBuffer.reset();

      if (!start(error)) {
        return false;
      }

      // Run basic checks on individual shapes in the index.
      S2Iterator<Cell> iter = index.iterator();
      for (int shapeId = 0; shapeId < index.getShapes().size(); ++shapeId) {
        S2Shape shape = index.getShapes().get(shapeId);
        if (shape != null) {
          if (!checkShape(iter, shape, shapeId, error)) {
            return false;
          }
        }
      }

      for (iter.restart(); !iter.done(); iter.next()) {
        setCurrentCell(iter);

        // Add two dimensional shape edges to the incident edge tracker to support checks for things
        // like crossing polygon chains and split interiors.
        for (S2ClippedShape clipped : currentCell().clippedShapes()) {
          int shapeId = clipped.shapeId();
          S2Shape shape = currentCell().shape(clipped);
          if (shape.dimension() < 2) {
            continue;
          }

          incidentEdgeTracker.startShape(shapeId);
          for (EdgeAndIdChain edge : currentCell().shapeEdges(shapeId)) {
            incidentEdgeTracker.addEdge(edge.edgeId(), edge.start(), edge.end());
          }
          incidentEdgeTracker.finishShape();
        }

        // Now notify that we're starting this cell.
        if (!startCell(error)) {
          return false;
        }

        // Iterate the shapes and edges of the cell.
        for (S2ClippedShape clipped : currentCell().clippedShapes()) {
          int shapeId = clipped.shapeId();
          S2Shape shape = currentCell().shape(clipped);
          if (!startShape(shape, clipped, error)) {
            return false;
          }

          for (EdgeAndIdChain edge : currentCell().shapeEdges(shapeId)) {
            if (!checkEdge(shape, clipped, edge, error)) {
              return false;
            }
          }

          if (!finishShape(shape, clipped, error)) {
            return false;
          }
        }
      }

      // Run any final checks and finish validation
      return finish(error);
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    // Subclass API

    /**
     * Starts the validation process; called once per query. Returns false if a validation error is
     * found.
     */
    protected boolean start(S2Error error) {
      return true;
    }

    /**
     * Validates each individual shape in the index; called once per shape. A reference to the
     * current iterator state is passed in. The function may reposition the iterator in order to do
     * shape checking. Returns false if a validation error is found.
     */
    protected boolean checkShape(S2Iterator<Cell> iter, S2Shape shape, int shapeId, S2Error error) {
      return true;
    }

    /**
     * Starts processing of a cell in the index; called once per cell. Returns false if a validation
     * error is found.
     */
    protected boolean startCell(S2Error error) {
      return true;
    }

    /**
     * Marks start of a clipped shape; called once per clipped shape in a cell. Returns false if a
     * validation error is found.
     */
    protected boolean startShape(S2Shape shape, S2ClippedShape clippedShape, S2Error error) {
      return true;
    }

    /**
     * Validates a single edge of a given shape; called at least once per edge of each shape.
     * Returns false if a validation error is found.
     */
    protected boolean checkEdge(
        S2Shape shape, S2ClippedShape clippedShape, EdgeAndIdChain chain, S2Error error) {
      return true;
    }

    /**
     * Marks end of a clipped shape; called once per clipped shape in a cell. Returns false if a
     * validation error is found.
     */
    protected boolean finishShape(S2Shape shape, S2ClippedShape clippedShape, S2Error error) {
      return true;
    }

    /**
     * Marks end of validation; called once per query. Returns false if a validation error is found.
     */
    protected boolean finish(S2Error error) {
      return true;
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////

    /** Default constructor. */
    public S2ValidationQueryBase() {}

    /** Returns the index we're validating. */
    protected S2ShapeIndex index() {
      assert index != null;
      return index;
    }

    /** Returns current incident edge map. */
    protected IncidentEdgeMap incidentEdges() {
      return incidentEdgeTracker.incidentEdges();
    }

    /** Returns a reference to the data for the current cell. */
    protected S2IndexCellData currentCell() {
      return cellBuffer;
    }

    /** Sets the current cell that we're operating on and loads its data. */
    private void setCurrentCell(S2Iterator<Cell> iter) {
      cellBuffer.loadCell(index, iter.id(), iter.entry());
    }
  }

  /**
   * This class represents the least strict class of geometry and corresponds to the requirements
   * for compatibility with S2BooleanOperation, with these specific constraints:
   *
   * <p>General
   *
   * <ul>
   *   <li>
   *   <li>Points must be unit magnitude (according to S2.isUnitLength()).
   *   <li>Points must be finite (neither inf or nan).
   *   <li>Edges must not have antipodal vertices.
   *   <li>Degenerate edges of the form {A,A} are allowed by default.
   *   <li>Reverse-duplicate edges of the form {A,B},{B,A} are allowed by default.
   *   <li>Shape chains must have at least one edge, excluding the full and empty polygon shapes,
   *       which may have at most one empty chain.
   * </ul>
   *
   * <p>Polygons
   *
   * <ul>
   *   <li>Polygon interiors must be disjoint from all other geometry. Polygon edges thus cannot
   *       cross any other edges, including at vertices. No geometry may test as contained in
   *       another polygon.
   *   <li>Polygons can never have duplicate edges even among different polygons.
   *   <li>Polygon edges must be connected and form a closed loop per chain.
   *   <li>Polygon chains must oriented so that the interior is always on the left.
   * </ul>
   *
   * <p>Polylines
   *
   * <ul>
   *   <li>Polyline edges can cross by default.
   *   <li>Polylines can have duplicate edges by default.
   * </ul>
   */
  public static class S2ValidQuery extends S2ValidationQueryBase {

    /**
     * Information on one vertex of an edge, including whether it's on the boundary of its shape.
     * TestVertex instances are mutable, and stored in a PooledList to be reused if possible.
     */
    private static class TestVertex {
      S2Point vertex = null;
      int edgeId = 0;
      int shapeId = 0;
      int dim = 0;
      boolean onBoundary = false;

      /** Sets the fields of this instance. */
      public void set(S2Point vertex, int edgeId, int shapeId, int dim, boolean onBoundary) {
        this.vertex = vertex;
        this.edgeId = edgeId;
        this.shapeId = shapeId;
        this.dim = dim;
        this.onBoundary = onBoundary;
      }
    }

    /** Used only by checkTouchesAreValid(). Declared once here and reused. */
    private final PooledList<TestVertex> testVertices = new PooledList<>(TestVertex::new);

    /**
     * Like MutableEdge, but also tracks the edge's chain id, edge id, previous id and incoming/
     * outgoing direction. EdgeWithInfo instances are mutable, and stored in a PooledList to be
     * reused if possible.
     */
    private static class EdgeWithInfo {
      S2Point startPoint = null;
      S2Point endPoint = null;
      int edgeId;
      int chainId;
      int prevEdgeId;
      int sign; // +1 for incoming and -1 for outgoing edges.

      /** Sets the fields of this instance. */
      public void set(S2Point start, S2Point end, int edgeId, int chainId, int prev, int sign) {
        this.startPoint = start;
        this.endPoint = end;
        this.edgeId = edgeId;
        this.chainId = chainId;
        this.prevEdgeId = prev;
        this.sign = sign;
      }

      /** Returns the start point of this edge. */
      public S2Point getStart() {
        return startPoint;
      }

      /** Returns the end point of this edge. */
      public S2Point getEnd() {
        return endPoint;
      }

      /** Returns true if this edge has the given point as either endpoint. */
      public boolean hasEndpoint(S2Point point) {
        return point.equalsPoint(startPoint) || point.equalsPoint(endPoint);
      }

      /** Returns true if this edge has the same start and end points as the other edge. */
      public boolean hasSamePoints(EdgeWithInfo other) {
        return startPoint.equalsPoint(other.startPoint) && endPoint.equalsPoint(other.endPoint);
      }

      /** Returns true if this edge has the reversed start and end points as the other edge. */
      public boolean hasReversePoints(EdgeWithInfo other) {
        return startPoint.equalsPoint(other.endPoint) && endPoint.equalsPoint(other.startPoint);
      }
    }

    /** Used only in checkVertexCrossings(). Declared once here and reused. */
    private final PooledList<EdgeWithInfo> edges = new PooledList<>(EdgeWithInfo::new);

    /**
     * Used only in checkVertexCrossings(). Declared once here and reused. Contains a mapping from
     * chainId => count. Typically, only a couple of chains will be incident at a given vertex, so a
     * linear search is acceptable.
     */
    private final IntPairVector chainSums = new IntPairVector();

    /**
     * Types of single-vertex touches allowed between shapes. This is configurable for each
     * combination of dimension. For example, we can configure polyline-polyline (1D-1D) and
     * polyline-polygon (1D-2D) touches separately.
     */
    protected enum TouchType {
      /** Allow no touches between shapes. */
      NONE(0b00),
      /** Interior point may touch the other shape. */
      INTERIOR(0b01),
      /** Boundary point may touch the other shape. */
      BOUNDARY(0b10),
      /** Allow any touches between shapes. */
      ANY(0b11);

      private final int value;

      private TouchType(int value) {
        this.value = value;
      }

      public boolean interiorMayTouch() {
        return (value & 0b01) != 0;
      }

      public boolean boundaryMayTouch() {
        return (value & 0b10) != 0;
      }

      public boolean matches(TouchType other) {
        return (value & other.value) != 0;
      }
    }

    /** A pair of TouchTypes, defining how two shapes are allowed to touch each other. */
    protected static class TouchTypePair {
      public static final TouchTypePair ANY_TO_ANY = TouchTypePair.of(TouchType.ANY, TouchType.ANY);

      /** How the first shape is allowed to touch the second shape. */
      final TouchType first;

      /** How the second shape is allowed to touch the first shape. */
      final TouchType second;

      /** Constructs a TouchTypePair with the given values. */
      public TouchTypePair(TouchType first, TouchType second) {
        this.first = first;
        this.second = second;
      }

      /** Constructs a TouchTypePair with the given values. */
      public static TouchTypePair of(TouchType first, TouchType second) {
        return new TouchTypePair(first, second);
      }

      /** Returns true if this TouchTypePair is equal to the other TouchTypePair. */
      public boolean isEqualTo(TouchTypePair other) {
        return (first.equals(other.first) && second.equals(other.second));
      }
    }

    /** Returns true if both touches match the given allowed touches. */
    private static boolean permittedTouches(
        TouchTypePair allowed, TouchType typeA, TouchType typeB) {
      return allowed.first.matches(typeA) && allowed.second.matches(typeB);
    }

    /**
     * Options for S2ValidQuery. Protected because options are only for subclasses to configure
     * behavior. Mutable, because they are created by the base class but subclasses may modify them.
     */
    protected static class Options {
      private boolean allowDegenerateEdges = true;
      private boolean allowDuplicatePolylineEdges = true;
      private boolean allowReverseDuplicates = true;
      private boolean allowPolylineIntgeriorCrossings = true;

      /** A matrix of allowed touches between geometry of different dimensions. */
      private final ArrayList<ArrayList<TouchTypePair>> allowedTouches;

      /** Constructs default options, which are maximally permissive. */
      public Options() {
        allowedTouches = new ArrayList<>();
        for (int i = 0; i < 3; ++i) {
          allowedTouches.add(new ArrayList<>());
          ArrayList<TouchTypePair> row = allowedTouches.get(i);
          for (int j = 0; j < 3; ++j) {
            row.add(TouchTypePair.ANY_TO_ANY);
          }
        }
      }

      /**
       * Returns a TouchType with fields set indicating what types of vertex from dimension A are
       * allowed to touch vertices from dimension B.
       */
      public TouchTypePair allowedTouches(int dima, int dimb) {
        Preconditions.checkArgument(dima >= 0);
        Preconditions.checkArgument(dima <= 2);
        Preconditions.checkArgument(dimb >= 0);
        Preconditions.checkArgument(dimb <= 2);

        // Just get the top half of the matrix.
        if (dima > dimb) {
          int tmp = dima;
          dima = dimb;
          dimb = tmp;
        }

        return allowedTouches.get(dima).get(dimb);
      }

      /**
       * Set the allowed touches based on geometry dimension. We require that the matrix of
       * combinations be symmetric, so setAllowedTouches(1, 2, X, Y) and setAllowedTouches(2, 1, Y,
       * X) are equivalent.
       */
      @CanIgnoreReturnValue
      public Options setAllowedTouches(int dima, int dimb, TouchTypePair types) {
        assert (0 <= dima && dima <= 2);
        assert (0 <= dimb && dimb <= 2);

        // Just set the top half of the matrix.
        if (dima > dimb) {
          int tmp = dima;
          dima = dimb;
          dimb = tmp;
        }

        allowedTouches.get(dima).set(dimb, types);
        return this;
      }

      /** Set allowed touches such that points are not allowed to touch any other geometry. */
      @CanIgnoreReturnValue
      public Options setNoPointTouchesAllowed() {
        for (int i = 0; i < 3; ++i) {
          setAllowedTouches(0, i, TouchTypePair.of(TouchType.NONE, TouchType.NONE));
        }
        return this;
      }

      /** Returns whether polylines can have duplicate edges. */
      public boolean allowDuplicatePolylineEdges() {
        return allowDuplicatePolylineEdges;
      }

      /** Sets whether polylines can have duplicate edges. */
      @CanIgnoreReturnValue
      public Options setAllowDuplicatePolylineEdges(boolean flag) {
        allowDuplicatePolylineEdges = flag;
        return this;
      }

      /** Returns whether polyline edges can cross. */
      public boolean allowPolylineInteriorCrossings() {
        return allowPolylineIntgeriorCrossings;
      }

      /** Sets whether polyline edges can cross. */
      @CanIgnoreReturnValue
      public Options setAllowPolylineInteriorCrossings(boolean flag) {
        allowPolylineIntgeriorCrossings = flag;
        return this;
      }

      /** Returns whether reverse duplicate edges are allowed. */
      public boolean allowReverseDuplicates() {
        return allowReverseDuplicates;
      }

      /** Sets whether reverse duplicate edges are allowed. */
      @CanIgnoreReturnValue
      public Options setAllowReverseDuplicates(boolean flag) {
        allowReverseDuplicates = flag;
        return this;
      }

      /** Returns whether degenerate edges are allowed. */
      public boolean allowDegenerateEdges() {
        return allowDegenerateEdges;
      }

      /** Sets whether degenerate edges are allowed. */
      @CanIgnoreReturnValue
      public Options setAllowDegenerateEdges(boolean flag) {
        allowDegenerateEdges = flag;
        return this;
      }
    }

    /** The options for this instance of S2ValidQuery. Mutable, but only by subclasses. */
    protected Options options;

    /** Constructor creates default options, which may be modified by subclasses. */
    public S2ValidQuery() {
      this.options = new Options();
    }

    /** Returns the options for this instance of S2ValidQuery. */
    public Options options() {
      return options;
    }

    @Override
    protected boolean checkShape(S2Iterator<Cell> iter, S2Shape shape, int shapeId, S2Error error) {
      // Verify that shape isn't outright lying to us about its dimension.
      int dim = shape.dimension();
      if (dim < 0 || dim > 2) {
        error.init(
            S2Error.Code.INVALID_DIMENSION, "Shape %d has invalid dimension: %d", shapeId, dim);
        return false;
      }

      IntVector chainsToCheck = new IntVector();
      MutableEdge prevEdge = new MutableEdge();
      MutableEdge currEdge = new MutableEdge();
      MutableEdge edge = new MutableEdge();
      MutableEdge last = new MutableEdge();

      for (int chainId = 0; chainId < shape.numChains(); ++chainId) {
        List<S2Point> chain = shape.chain(chainId);
        int chainStart = shape.getChainStart(chainId);
        int chainLength = shape.getChainLength(chainId);

        // Check that the first and last edges in a polygon chain connect to close the chain. This
        // is true even for degenerate chains with one edge that's a point, or two edges that are
        // reverse duplicates. Both are considered closed.
        if (dim == 2 && chain.size() > 0) {
          int edgeId = shape.getChainStart(chainId);
          shape.getEdge(edgeId, currEdge);
          shape.getEdge(shape.prevEdgeWrap(edgeId), prevEdge);

          if (!prevEdge.b.equalsPoint(currEdge.a)) {
            error.init(
                S2Error.Code.LOOP_NOT_ENOUGH_VERTICES,
                "Chain %d of shape %d isn't closed",
                chainId,
                shapeId);
            return false;
          }
        }

        for (int offset = 0; offset < chainLength; ++offset) {
          shape.getChainEdge(chainId, offset, edge);

          // Check that coordinates aren't inf/nan.
          if (!validPoint(edge.a) || !validPoint(edge.b)) {
            error.init(S2Error.Code.INVALID_VERTEX, "Shape %d has invalid coordinates", shapeId);
            return false;
          }

          // Check that vertices are unit length.
          if (!S2.isUnitLength(edge.a) || !S2.isUnitLength(edge.b)) {
            error.init(
                S2Error.Code.NOT_UNIT_LENGTH, "Shape %d has non-unit length vertices", shapeId);
            return false;
          }

          // (Optional) check polyline and polygon edges for degeneracy.
          if (dim > 0 && !options().allowDegenerateEdges()) {
            if (edge.isDegenerate()) {
              error.init(
                  S2Error.Code.DUPLICATE_VERTICES,
                  "Shape %d: chain %d, edge %d is degenerate",
                  shapeId,
                  chainId,
                  chainStart + offset);
              return false;
            }
          }

          // Check that edge doesn't have antipodal vertices.
          if (edge.a.equalsPoint(edge.b.neg())) {
            error.init(
                S2Error.Code.ANTIPODAL_VERTICES,
                "Shape %d has adjacent antipodal vertices",
                shapeId);
            return false;
          }

          // Check that chain edges are connected for polylines and polygons.
          if (dim > 0 && chain.size() >= 2 && offset > 0) {
            shape.getChainEdge(chainId, offset - 1, last);
            if (!last.b.equalsPoint(edge.a)) {
              error.init(
                  S2Error.Code.NOT_CONTINUOUS,
                  "Chain %d of shape %d has neighboring edges that don't connect.",
                  chainId,
                  shapeId);
              return false;
            }
          }
        }

        // The rest of the checks are for non-empty polygon chains only.
        if (dim != 2 || chain.size() == 0) {
          continue;
        }

        // We need at least two distinct points in a chain before we can check its orientation vs
        // the cell center. Scan until we find a vertex different than the first.
        int uniqueCount = 1;
        shape.getChainEdge(chainId, 0, edge);
        S2Point first = edge.a;
        for (int i = 0; i < shape.getChainLength(chainId); i++) {
          S2Point vertex = shape.getChainVertex(chainId, i);
          if (!vertex.equalsPoint(first)) {
            ++uniqueCount;
            break;
          }
        }

        // Only a single unique point. A degenerate edge will never test as a vertex crossing
        // (because 3 out of 4 vertices to S2EdgeUtil.vertexCrossing() would be equal making it
        // false), so they can't toggle interior state and we can ignore them.
        if (uniqueCount == 1) {
          continue;
        }
        chainsToCheck.add(chainId);
      }

      OfInt chainsIter = chainsToCheck.intIterator();
      while (chainsIter.hasNext()) {
        int chainId = chainsIter.nextInt();
        if (!checkChainOrientation(iter, shape, shapeId, chainId, error)) {
          return false;
        }
      }

      return true;
    }

    @Override
    protected boolean startCell(S2Error error) {
      if (!checkForDuplicateEdges(error) || !checkForInteriorCrossings(error)) {
        return false;
      }

      if (!checkTouchesAreValid(error)) {
        return false;
      }

      return true;
    }

    @Override
    protected boolean checkEdge(
        S2Shape shape, S2ClippedShape clipped, EdgeAndIdChain edge, S2Error error) {
      int dim = shape.dimension();

      // For points, we can check that they're not contained in any other polygons locally within
      // the cell by crossing edges to the cell center.
      if (dim == 0 && pointContained(clipped.shapeId(), edge.start(), error)) {
        return false;
      }

      // Edge is OK
      return true;
    }

    @Override
    protected boolean finish(S2Error error) {
      // We've checked edges having interiors on the right, and for crossings at interior points.
      // The only case left is to check for chains that cross at a vertex.
      for (Map.Entry<IncidentEdgeKey, Set<Integer>> entry : incidentEdges().entrySet()) {
        S2Shape shape = index().getShapes().get(entry.getKey().shapeId);
        if (shape.dimension() == 2) {
          if (!checkVertexCrossings(
              entry.getKey().vertex, shape, entry.getKey().shapeId, entry.getValue(), error)) {
            return false;
          }
        }
      }

      // If we get to this point we know that polygon edges don't cross any other edges and that
      // edges are properly oriented with the interior on the left.
      //
      // Since edges don't cross, any given chain must be entirely inside or outside any other
      // polygons. Thus, to determine that polygon interiors are disjoint, we only have to check one
      // vertex of each chain in each shape for containment.
      //
      // We use the OPEN containment model because we check elsewhere if the vertex lands on another
      // vertex.
      S2ContainsPointQuery.Options options = new S2ContainsPointQuery.Options(S2VertexModel.OPEN);
      S2ContainsPointQuery query = new S2ContainsPointQuery(index(), options);
      for (int shapeId = 0; shapeId < index().getShapes().size(); ++shapeId) {
        S2Shape shape = index().getShapes().get(shapeId);
        if (shape == null || shape.dimension() == 0) {
          continue;
        }

        MutableEdge edge = new MutableEdge();
        for (int chain = 0; chain < shape.numChains(); ++chain) {
          if (shape.getChainLength(chain) < 1) {
            continue;
          }

          shape.getChainEdge(chain, 0, edge);
          S2Point vertex = edge.a;
          if (query.contains(vertex)) {
            error.init(
                S2Error.Code.OVERLAPPING_GEOMETRY,
                "Shape %d has one or more edges contained in another shape.",
                shapeId);
            return false;
          }
        }
      }
      return true;
    }

    /**
     * Checks that edges of the given shape incident on vertex are ordered such that the incident
     * chains do not cross.
     *
     * <p>We can check this by looking at all the incident edges and making sure that, for each
     * incoming edge, as we move counter-clockwise around the vertex, we encounter matching pairs of
     * incoming/outgoing edges for each chain.
     *
     * <p>Returns true if chains do not cross at the vertex, false otherwise.
     */
    private boolean checkVertexCrossings(
        S2Point vertex, S2Shape shape, int shapeId, Iterable<Integer> edgeIds, S2Error error) {

      // Aggregate edges of the current shape in "edges", a PooledList<EdgeWithInfo>, and sort them
      // CCW around the vertex.
      edges.clear();
      MutableEdge edge = new MutableEdge();
      ChainPosition pos = new ChainPosition();
      for (int edgeId : edgeIds) {
        shape.getChainPosition(edgeId, pos);
        shape.getEdge(edgeId, edge);

        edges.add().set(
            edge.a,
            edge.b,
            edgeId,
            pos.chainId,
            shape.prevEdgeWrap(edgeId), //
            edge.a.equalsPoint(vertex) ? -1 : +1);
      }
      sortCcw(vertex, edges);

      for (int i = 0; i < edges.size(); ++i) {
        EdgeWithInfo curr = edges.get(i);

        // Skip forward to next outgoing edge.
        if (curr.sign > 0) {
          continue;
        }

        // Scan until we find our incoming edge and tally chain counts.
        chainSums.clear(); // Map from chain id to sum of signs of edges.
        int j;
        for (j = 1; j < edges.size(); ++j) {
          EdgeWithInfo edgeInfo = edges.get((i + j) % edges.size());
          if (curr.chainId == edgeInfo.chainId && curr.prevEdgeId == edgeInfo.edgeId) {
            for (int csi = 0; csi < chainSums.size(); ++csi) {
              if (chainSums.getSecond(csi) != 0) {
                error.init(
                    S2Error.Code.OVERLAPPING_GEOMETRY,
                    "Shape %d has one or more chains that cross at a vertex",
                    shapeId);
                return false;
              }
            }
            break;
          }
          // Find the index of the chain sum for this chain, or if it's new, add it with initial
          // value zero.
          int chainSumIndex = findOrAddFirst(chainSums, edgeInfo.chainId);
          // Increment the sum by the sign of the edge.
          int sum = chainSums.getSecond(chainSumIndex) + edgeInfo.sign;
          chainSums.setSecond(chainSumIndex, sum);
        }

        // We went all the way around, where's the incoming edge?
        assert j != edges.size();
      }

      return true;
    }

    /**
     * Does a linear search and returns the lowest index in the given IntPairVector where the first
     * element of the stored pairs equals the given 'firstValue'. If no such index exists, adds a
     * new pair to the end of the vector, setting the first element to the given 'firstValue', and
     * the second element to zero, and returns the new index.
     */
    private static int findOrAddFirst(IntPairVector pairVector, int firstValue) {
      int firstIndex = -1;
      for (int i = 0; i < pairVector.size(); i++) {
        if (pairVector.getFirst(i) == firstValue) {
          firstIndex = i;
          break;
        }
      }

      if (firstIndex == -1) {
        firstIndex = pairVector.size();
        pairVector.add(firstValue, 0);
      }
      return firstIndex;
    }

    /**
     * Sorts a PooledList of EdgeWithInfo in counter-clockwise order around an origin point. By
     * itself, ordering this way is ambiguous in terms of the starting point, so we use the first
     * edge to form a reference point to form a total ordering.
     *
     * <p>Reverse duplicate edges are ordered so that the one with origin as the start point comes
     * before the other.
     *
     * <p>The edges must all contain origin as one of their vertices.
     */
    public static void sortCcw(S2Point origin, PooledList<EdgeWithInfo> edges) {
      EdgeWithInfo first = edges.get(0);
      S2Point firstVertex =
          first.getStart().equalsPoint(origin) ? first.getEnd() : first.getStart();
      Preconditions.checkArgument(!firstVertex.equalsPoint(origin));

      for (EdgeWithInfo edge : edges) {
        Preconditions.checkArgument(edge.hasEndpoint(origin));
      }

      Comparator<EdgeWithInfo> comparator =
          (EdgeWithInfo a, EdgeWithInfo b) -> {
            // orderedCCW will return true if a == b, which will violate the irreflexivity
            // requirement of a strict weak ordering.
            if (a.hasSamePoints(b)) {
              return -1;
            }
            // Order reverse duplicates so that the one with edge.getStart() == origin is first.
            if (a.hasReversePoints(b)) {
              return a.getStart().equalsPoint(origin) ? 1 : -1;
            }
            // If edge a is the first edge, then it is always first in the ordering.
            if (a.hasSamePoints(first)) {
              return 1;
            }
            // If edge b is the first edge, then edge a always comes after it.
            if (b.hasSamePoints(first)) {
              return -1;
            }

            // Otherwise check orientation of vertices.
            S2Point aPoint = a.getStart().equalsPoint(origin) ? a.getEnd() : a.getStart();
            S2Point bPoint = b.getStart().equalsPoint(origin) ? b.getEnd() : b.getStart();
            return orderedCCW(firstVertex, aPoint, bPoint, origin) ? 1 : -1;
          };

      edges.sort(comparator);
    }

    /**
     * Checks that a point is fully outside any polygons in the current cell.
     *
     * <p>Returns false if a point is not interior to any polygons, true otherwise, with error
     * populated.
     */
    private boolean pointContained(int shapeId, S2Point point, S2Error error) {
      S2IndexCellData cell = currentCell();
      for (S2ClippedShape clipped : cell.clippedShapes()) {
        if (clipped.shapeId() == shapeId) {
          continue;
        }

        S2Shape shape = index().getShapes().get(clipped.shapeId());

        // Only check for containment in polygons.
        if (shape.dimension() != 2) {
          continue;
        }

        if (cell.shapeContains(clipped, point)) {
          error.init(
              S2Error.Code.OVERLAPPING_GEOMETRY,
              "Shape %d has one or more edges contained in another shape.",
              shapeId);
          return true;
        }
      }

      return false;
    }

    /**
     * Checks that a given chain of a shape is oriented properly relative to one cell center in the
     * index. We only need to check one center because we assume that the transition between cell
     * indices is consistent, thus we just need to make sure that the interior state isn't flipped.
     *
     * <p>Returns false when the chain isn't properly oriented with S2Error set with the details,
     * true otherwise.
     *
     * <p>REQUIRES: Shape must be two dimensional.
     *
     * <p>REQUIRES: Chain must be closed.
     *
     * <p>REQUIRES: Chain edges must be connected (no gaps).
     *
     * <p>REQUIRES: Chain must have at least three distinct points (cover some area).
     */
    private boolean checkChainOrientation(
        S2Iterator<Cell> iter, S2Shape shape, int shapeId, int chainId, S2Error error) {
      int chainLength = shape.getChainLength(chainId);

      // Given that:
      //  1. Edges in the chain are connected continuously.
      //  2. The chain is closed.
      //  3. The chain has at least two distinct points.
      //
      // Then we can test whether the chain is oriented properly relative to the cell center by
      // testing one edge of the chain for proper orientation.

      S2ContainsVertexQuery query = new S2ContainsVertexQuery();
      MutableEdge edge = new MutableEdge();
      for (int offset = 0; offset < chainLength; ++offset) {
        S2Point vertex = shape.getChainVertex(chainId, offset);
        query.init(vertex);

        // Seek to the cell containing vertex and get the clipped shape.
        if (!iter.locate(vertex)) {
          error.init(S2Error.Code.DATA_LOSS, "Shape vertex was not indexed");
          return false;
        }
        S2Point center = iter.id().toPoint();

        S2ClippedShape clipped = iter.entry().findClipped(shapeId);
        assert clipped != null;

        // Compute winding number and vertex sign at the same time.
        int winding = clipped.containsCenter() ? 1 : 0;
        EdgeCrosser crosser = new EdgeCrosser(center, vertex);

        for (int i = 0; i < clipped.numEdges(); ++i) {
          shape.getEdge(clipped.edge(i), edge);

          // Tally up the total change in winding number from center to vertex.
          winding += crosser.signedEdgeOrVertexCrossing(edge.a, edge.b);

          // Include any edges incident on vertex in the contains vertex query.
          if (vertex.equalsPoint(edge.a)) {
            // Clipped shape edge is from 'vertex' to 'edge.b'.
            query.addOutgoing(edge.b);
          } else if (vertex.equalsPoint(edge.b)) {
            // Clipped shape edge is from 'edge.a' to 'vertex'.
            query.addIncoming(edge.a);
          }
        }

        boolean duplicates = query.duplicateEdges();
        int sign = 0;

        // If we have a sign of zero on the vertex, all the edges incident on it were reverse
        // duplicates and we can't use it to test orientation, continue trying to find another
        // vertex.
        if (!duplicates) {
          sign = query.containsSign();
          if (sign == 0) {
            continue;
          }
        }

        // The sign bit obtained by crossing edges should be consistent with the sign produced by
        // the S2ContainsVertexQuery.
        if (duplicates || winding != (sign < 0 ? 0 : 1)) {
          error.init(
              S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS,
              "Shape %d has one or more edges with interior on the right.",
              shapeId);
          return false;
        }
        return true;
      }
      return true;
    }

    /** Checks if any edge in the current cell is a duplicate (or reverse duplicate). */
    private boolean checkForDuplicateEdges(S2Error error) {
      int dim0 = options().allowDuplicatePolylineEdges() ? 2 : 1;
      int dim1 = 2;

      // This is O(N^2) but cells don't have many edges in them and benchmarks show this to be
      // faster than trying to sort the whole array ahead of time.
      List<EdgeAndIdChain> edges = currentCell().dimRangeEdges(dim0, dim1);
      int numEdges = edges.size();
      for (int i = 0; i < numEdges; ++i) {
        for (int j = i + 1; j < numEdges; ++j) {
          boolean duplicate = edges.get(i).isEqualTo(edges.get(j));
          if (!options().allowReverseDuplicates()) {
            duplicate |= edges.get(i).isReverseOf(edges.get(j));
          }

          if (duplicate) {
            error.init(
                S2Error.Code.OVERLAPPING_GEOMETRY, "One or more duplicate polygon edges detected");
            return false;
          }
        }
      }

      return true;
    }

    /** Checks if any edges in the current cell have an interior crossing. */
    private boolean checkForInteriorCrossings(S2Error error) {
      // Get all the polyline and polygon edges.
      List<EdgeAndIdChain> edges = currentCell().dimRangeEdges(1, 2);

      // If we're allowing polyline edges to cross polyline edges, then we only have to check
      // against polygon edges.
      int checkStart = 0;
      if (options().allowPolylineInteriorCrossings()) {
        checkStart = currentCell().dimEdges(1).size();
      }

      // This can happen when we're allowing polyline crossings and only have polylines, there's
      // nothing to do.
      if (checkStart >= edges.size()) {
        return true;
      }

      int numEdges = edges.size();
      for (int i = 0; i + 1 < numEdges; ++i) {
        // We never have to check against edges at a lower index, because if we intersect them,
        // we'll have already checked that at this point.
        int j = max(checkStart, i + 1);

        // We can skip adjacent (sharing an endpoint) edges.
        if (edges.get(i).end().equalsPoint(edges.get(j).start())) {
          if (++j >= numEdges) {
            break;
          }
        }

        EdgeCrosser crosser = new EdgeCrosser(edges.get(i).start(), edges.get(i).end());
        for (; j < numEdges; ++j) {
          if (crosser.c() == null || !crosser.c().equalsPoint(edges.get(j).start())) {
            crosser.restartAt(edges.get(j).start());
          }

          if (crosser.robustCrossing(edges.get(j).end()) > 0) {
            error.init(
                S2Error.Code.OVERLAPPING_GEOMETRY,
                "Chain %d edge %d crosses chain %d edge %d",
                edges.get(i).chainId(),
                edges.get(i).offset(),
                edges.get(i).chainId(),
                edges.get(i).offset());
            return false;
          }
        }
      }
      return true;
    }

    /**
     * Checks the shapes in the current cell to ensure that any touch points are allowed under the
     * configured semantics.
     *
     * <p>Returns false if any vertices touch at an invalid point with the error value set, true
     * otherwise.
     */
    private boolean checkTouchesAreValid(S2Error error) {
      boolean[] needDim = {true, true, true};
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          TouchTypePair allowed = options().allowedTouches(i, j);
          boolean anyAllowed = allowed.isEqualTo(TouchTypePair.ANY_TO_ANY);
          needDim[i] &= !anyAllowed;
        }
      }

      // If all touches are allowed, then we don't have to check, easy.
      if (!needDim[0] && !needDim[1] && !needDim[2]) {
        return true;
      }

      // Gather unique vertices from each dimension that needs to be checked.
      testVertices.clear();
      for (S2ClippedShape clipped : currentCell().clippedShapes()) {
        int shapeId = clipped.shapeId();
        S2Shape shape = index().getShapes().get(shapeId);
        int dim = shape.dimension();

        // Skip if we're not checking this dimension.
        if (!needDim[dim]) {
          continue;
        }

        for (EdgeAndIdChain edge : currentCell().shapeEdges(shapeId)) {
          // For polylines, we have to handle the start and ending edges specially, since we can
          // have open chains that have boundary points.
          if (dim == 1) {
            boolean onBoundary = polylineVertexIsBoundaryPoint(shape, edge, 0);
            testVertices.add().set(edge.start(), edge.edgeId(), shapeId, dim, onBoundary);

            // Check vertex 1 too only if it's a boundary point, otherwise we would test v1 twice
            // after we grab v0 of the next edge.
            onBoundary = polylineVertexIsBoundaryPoint(shape, edge, 1);
            if (onBoundary) {
              testVertices.add().set(edge.end(), edge.edgeId(), shapeId, dim, true);
            }
          } else {
            // All polygon vertices are on the boundary, all point vertices are not.
            testVertices.add().set(edge.start(), edge.edgeId(), shapeId, dim, dim == 2);
          }
        }
      }

      // For each test vertex, scan over the other shapes and verify that any vertex touches are
      // allowed under the current semantics.
      for (TestVertex testPoint : testVertices) {
        for (S2ClippedShape clipped : currentCell().clippedShapes()) {
          int shapeId = clipped.shapeId();
          S2Shape shape = index().getShapes().get(shapeId);
          int dim = shape.dimension();

          for (EdgeAndIdChain edge : currentCell().shapeEdges(shapeId)) {
            // Don't compare an edge against itself.
            if (testPoint.shapeId == shapeId && testPoint.edgeId == edge.edgeId()) {
              continue;
            }

            // Figure out which (if any) vertex of this edge we touched.
            int vertIndex = -1;
            if (testPoint.vertex.equalsPoint(edge.start())) {
              vertIndex = 0;
            }
            if (testPoint.vertex.equalsPoint(edge.end())) {
              vertIndex = 1;
            }
            if (vertIndex < 0) {
              continue;
            }

            // Closed polylines are always allowed, so if we get a hit on the same shape and it's a
            // polyline, check that it's not from the closure.
            if (testPoint.shapeId == shapeId && dim == 1) {
              if (vertIndex == 0 && shape.prevEdgeWrap(edge.edgeId()) == testPoint.edgeId) {
                continue;
              }
              if (vertIndex == 1 && shape.nextEdgeWrap(edge.edgeId()) == testPoint.edgeId) {
                continue;
              }
            }

            // Points vertices are always interior, Polygon vertices are always boundary, and it may
            // be either for polylines.
            boolean onBoundary = (dim == 2);
            if (dim == 1) {
              onBoundary = polylineVertexIsBoundaryPoint(shape, edge, vertIndex);
            }

            TouchType typeA = TouchType.INTERIOR;
            if (testPoint.onBoundary) {
              typeA = TouchType.BOUNDARY;
            }

            TouchType typeB = TouchType.INTERIOR;
            if (onBoundary) {
              typeB = TouchType.BOUNDARY;
            }

            // We require that touches be symmetric, so that the touch is invalid only if it's
            // invalid from A->B and B->A. This lets us support behavior such as requiring that one
            // or the other point be an interior point.
            TouchTypePair allowed = options().allowedTouches(testPoint.dim, dim);
            if (!permittedTouches(allowed, typeA, typeB)
                && !permittedTouches(allowed, typeB, typeA)) {
              error.init(
                  S2Error.Code.OVERLAPPING_GEOMETRY,
                  "Index has geometry with invalid vertex touches.");
              return false;
            }
          }
        }
      }
      return true;
    }

    /**
     * Returns true if a vertex (0 or 1) of an edge of a polyline is a boundary point or not. This
     * is only true if the vertex is either the start or end point of a chain and the chain is open.
     * Returns false otherwise.
     *
     * <p>REQUIRES: vertex == 0 or vertex == 1
     */
    private static boolean polylineVertexIsBoundaryPoint(
        S2Shape shape, EdgeAndIdChain edge, int vertex) {
      assert vertex == 0 || vertex == 1;

      if (edge.offset() == 0 && vertex == 0) {
        // The first vertex of the first edge of the chain.
        return shape.prevEdgeWrap(edge.edgeId()) == -1;
      } else if (edge.offset() == shape.getChainLength(edge.chainId()) - 1 && vertex == 1) {
        // The last vertex of the last edge of the chain.
        return shape.nextEdgeWrap(edge.edgeId()) == -1;
      }

      return false;
    }
  }

  /**
   * This class represents a semantic class of geometry that's compatible with the requirements of
   * the S2Polygon and S2Polyline isValid() methods. It extends the S2ValidQuery requirements,
   * specifically:
   *
   * <p>General
   *
   * <ul>
   *   <li>Degenerate edges of the form {A,A} are not allowed.
   *   <li>All the shapes in an S2ShapeIndex must be the same dimensionality.
   *   <li>Duplicate vertices in a chain are not allowed. I.e. a chain cannot touch itself even at
   *       one point.
   *   <li>Different chains may touch, but only in such a way they don't create duplicate edges.
   * </ul>
   */
  public static class S2LegacyValidQuery extends S2ValidQuery {

    /** Constructs a new S2LegacyValidQuery, modifying options to match legacy semantics. */
    public S2LegacyValidQuery() {
      options.setAllowDegenerateEdges(false).setAllowReverseDuplicates(false);
    }

    @Override
    protected boolean start(S2Error error) {
      if (!super.start(error)) {
        return false;
      }

      // We can't mix dimensions under legacy semantics.
      int dim = -1;
      for (S2Shape shape : index().getShapes()) {
        if (dim < 0) {
          dim = shape.dimension();
        }

        if (dim != shape.dimension()) {
          error.init(
              S2Error.Code.INVALID_DIMENSION,
              "Mixed dimensional geometry is invalid for legacy semantics.");
          return false;
        }
      }

      return true;
    }

    @Override
    protected boolean checkShape(S2Iterator<Cell> iter, S2Shape shape, int shapeId, S2Error error) {
      // Count the number of empty chains. Non-empty chains must have at least three vertices.
      if (shape.dimension() == 2) {
        boolean hasEmptyLoops = false;
        for (int c = 0; c < shape.numChains(); ++c) {
          int chainLength = shape.getChainLength(c);
          if (chainLength == 0) {
            // Empty or full.
            hasEmptyLoops = true;
          } else if (chainLength < 3) {
            error.init(
                S2Error.Code.LOOP_NOT_ENOUGH_VERTICES,
                "Shape %d has a non-empty chain with less than three edges.",
                shapeId);
            return false;
          }
        }

        if (hasEmptyLoops && shape.numChains() > 1) {
          error.init(
              S2Error.Code.POLYGON_EMPTY_LOOP, "Shape %d has too many empty chains", shapeId);
          return false;
        }
      }

      if (!super.checkShape(iter, shape, shapeId, error)) {
        return false;
      }

      return true;
    }

    @Override
    protected boolean startCell(S2Error error) {
      // Check for duplicate vertices within a chain.
      S2IndexCellData cell = currentCell();
      for (S2ClippedShape clipped : cell.clippedShapes()) {
        List<EdgeAndIdChain> edges = cell.shapeEdges(clipped.shapeId());

        for (int i = 0; i < edges.size(); ++i) {
          for (int j = i + 1; j < edges.size(); ++j) {
            if (edges.get(j).chainId() != edges.get(i).chainId()) {
              continue;
            }

            if (edges.get(j).start().equalsPoint(edges.get(i).start())) {
              error.init(
                  S2Error.Code.DUPLICATE_VERTICES,
                  "Chain %d of shape %d has duplicate vertices",
                  edges.get(i).chainId(),
                  clipped.shapeId());
              return false;
            }
          }
        }
      }

      return super.startCell(error);
    }
  }

  /**
   * Sorts a List of S2Edges in counter-clockwise order around an origin point. By itself, ordering
   * this way is ambiguous in terms of the starting point, so we take an edge to form a reference
   * point to form a total ordering.
   *
   * <p>Reverse duplicate edges are ordered so that the one with origin as v0 comes before the
   * other.
   *
   * <p>The edges and first edge must all contain origin as one of their vertices.
   */
  public static <E extends S2Edge> void sortEdgesCcw(S2Point origin, E first, List<E> data) {
    Preconditions.checkArgument(
        first.getStart().equalsPoint(origin) || first.getEnd().equalsPoint(origin));
    S2Point firstVertex = first.getStart().equalsPoint(origin) ? first.getEnd() : first.getStart();
    Preconditions.checkArgument(!firstVertex.equalsPoint(origin));

    for (S2Edge edge : data) {
      Preconditions.checkArgument(
          edge.getStart().equalsPoint(origin) || edge.getEnd().equalsPoint(origin));
    }

    SortableCollection sortable =
        new SortableCollection() {
          @Override
          public int size() {
            return data.size();
          }

          @Override
          public void truncate(int end) {
            if (end < data.size()) {
              data.subList(end, data.size()).clear();
            }
          }

          @Override
          public boolean less(int leftIndex, int rightIndex) {
            E a = data.get(leftIndex);
            E b = data.get(rightIndex);

            // orderedCCW will return true if a == b, which will violate the irreflexivity
            // requirement of a strict weak ordering.
            if (areEqual(a, b)) {
              return false;
            }
            // Order reverse duplicates so that the one with edge.getStart() == origin is first.
            if (areReversed(a, b)) {
              return a.getStart().equalsPoint(origin);
            }
            // If edge a is the first edge, then it is always first in the ordering.
            if (areEqual(a, first)) {
              return true;
            }
            // If edge b is the first edge, then edge a always comes after it.
            if (areEqual(b, first)) {
              return false;
            }

            // Otherwise check orientation of vertices.
            S2Point aPoint = a.getStart().equalsPoint(origin) ? a.getEnd() : a.getStart();
            S2Point bPoint = b.getStart().equalsPoint(origin) ? b.getEnd() : b.getStart();
            return orderedCCW(firstVertex, aPoint, bPoint, origin);
          }

          @Override
          public void swap(int leftIndex, int rightIndex) {
            E tmp = data.get(leftIndex);
            data.set(leftIndex, data.get(rightIndex));
            data.set(rightIndex, tmp);
          }
        };

    sortable.sort();
  }

  private static boolean areEqual(S2Edge a, S2Edge b) {
    return a.getStart().equalsPoint(b.getStart()) && a.getEnd().equalsPoint(b.getEnd());
  }

  private static boolean areReversed(S2Edge a, S2Edge b) {
    return a.getStart().equalsPoint(b.getEnd()) && a.getEnd().equalsPoint(b.getStart());
  }

  /**
   * Returns true if the given S2Point is valid, meaning none of its components are infinite or NaN.
   */
  private static boolean validPoint(S2Point p) {
    return Double.isFinite(p.getX()) && Double.isFinite(p.getY()) && Double.isFinite(p.getZ());
  }

  private S2ValidationQueries() {}
}
