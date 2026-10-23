/*
 * Copyright 2026 Google Inc.
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

import static com.google.common.base.Preconditions.checkState;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.geometry.S2Builder.SnapFunction;
import com.google.common.geometry.S2BuilderSnapFunctions.IdentitySnapFunction;
import com.google.common.geometry.S2WindingOperation.WindingRule;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.Arrays;
import java.util.List;

/**
 * Computes boolean operations on collections of polygons or {@link S2ShapeIndex} instances
 * efficiently using {@link S2WindingOperation}.
 *
 * <h2>Determinism and Polygon Order Invariance</h2>
 *
 * <p>Under the hood, output geometry is assembled using {@link S2Builder}, and the output's vertex
 * sequence is sensitive to the order in which input shapes and edges are processed:
 *
 * <ul>
 *   <li><b>Static helpers:</b> Convenience methods such as {@link #union(Iterable,
 *       S2WindingOperation.Options)} and {@link #intersection(S2Polygon...)} pre-sort input
 *       polygons by {@link S2Polygon#compareTo}, guaranteeing that the results are {@link
 *       S2Polygon#equals equal} regardless of the input order.
 *   <li><b>Instance and Index APIs:</b> When using {@link #add(S2ShapeIndex)} or adding shapes
 *       sequentially to an {@link S2MultiOperation} instance, inputs are not automatically sorted.
 *       Callers requiring deterministic output must ensure that inputs are provided in a
 *       consistent, deterministic order (both the iteration order of the shapes and the cyclic
 *       vertex ordering of each shape's loops).
 * </ul>
 */
public class S2MultiOperation {

  /** The type of boolean operation to perform. */
  public enum OpType {
    /** Computes the union of all added shapes. Regions covered by at least one shape are kept. */
    UNION,

    /** Computes the intersection of all added shapes. Regions covered by ALL shapes are kept. */
    INTERSECTION,

    /**
     * Computes the symmetric difference of the added shapes.
     *
     * <p>A region belongs to the result if it is covered by an <b>odd number</b> of the input
     * shapes. For two shapes, this is the classic XOR operation (regions in A or B, but not both).
     */
    SYMMETRIC_DIFFERENCE
  }

  /** The type of boolean operation being performed. */
  public final OpType opType;

  /** Returns the type of boolean operation being performed. */
  public OpType opType() {
    return opType;
  }

  private final S2WindingOperation windingOp = new S2WindingOperation();
  private final S2BuilderLayer resultLayer;
  private int refWinding = 0;
  private int num2dShapes = 0;

  /** Default snap function using {@link S2EdgeUtil#INTERSECTION_MERGE_RADIUS}. */
  public static final SnapFunction DEFAULT_SNAP_FUNCTION =
      new IdentitySnapFunction(S2EdgeUtil.INTERSECTION_MERGE_RADIUS);

  /** Returns the default options, configured with {@link #DEFAULT_SNAP_FUNCTION}. */
  public static S2WindingOperation.Options defaultOptions() {
    return new S2WindingOperation.Options(DEFAULT_SNAP_FUNCTION);
  }

  /**
   * Constructs a new S2MultiOperation with default options ({@link #DEFAULT_SNAP_FUNCTION}).
   *
   * @param opType The type of boolean operation to perform (UNION, INTERSECTION, etc.).
   * @param resultLayer The layer where the resulting geometry will be built.
   */
  public S2MultiOperation(OpType opType, S2BuilderLayer resultLayer) {
    this(opType, resultLayer, defaultOptions());
  }

  /**
   * Constructs a new S2MultiOperation with specific options.
   *
   * @param opType The type of boolean operation to perform.
   * @param resultLayer The layer where the resulting geometry will be built.
   * @param options Options for the underlying S2WindingOperation.
   */
  public S2MultiOperation(
      OpType opType, S2BuilderLayer resultLayer, S2WindingOperation.Options options) {
    this.opType = opType;
    this.resultLayer = resultLayer;
    this.windingOp.init(resultLayer, options);
  }

  /** Clears the underlying operation. */
  public void clear() {
    this.refWinding = 0;
    this.num2dShapes = 0;
    this.windingOp.init(resultLayer, windingOp.options());
  }

  /** Sets the options for the underlying {@link S2WindingOperation}. */
  @CanIgnoreReturnValue
  public S2MultiOperation setOptions(S2WindingOperation.Options options) {
    checkState(num2dShapes == 0, "Options cannot be changed after shapes are added.");
    this.windingOp.init(resultLayer, options);
    return this;
  }

  /**
   * Adds shapes from an S2ShapeIndex to the operation.
   *
   * @param index The index containing the shapes to add.
   * @throws IllegalArgumentException if any non-2D shapes are present in the index.
   */
  public void add(S2ShapeIndex index) {
    addShapeIndexInternal(index, false);
  }

  /** Adds a polygon to the operation. */
  public void add(S2Polygon polygon) {
    addPolygonInternal(polygon, false);
  }

  /** Adds the complement of shapes from an S2ShapeIndex to the operation. */
  public void addComplement(S2ShapeIndex index) {
    addShapeIndexInternal(index, true);
  }

  /** Adds the complement of a polygon to the operation. */
  public void addComplement(S2Polygon polygon) {
    addPolygonInternal(polygon, true);
  }

  // Ingests an S2Polygon directly, traversing its loops and XORing loop.containsOrigin(). This
  // intentionally avoids building polygon.index() or allocating a temporary S2Polygon.Shape
  // wrapper.
  private void addPolygonInternal(S2Polygon polygon, boolean complement) {
    num2dShapes++;
    boolean polygonContainsOrigin = false;
    for (S2Loop loop : polygon.getLoops()) {
      polygonContainsOrigin ^= loop.containsOrigin();
      if (loop.isEmptyOrFull()) {
        continue;
      }
      List<S2Point> vertices = loop.orientedVertices();
      windingOp.addLoop(complement ? Lists.reverse(vertices) : vertices);
    }
    if (complement ^ polygonContainsOrigin) {
      refWinding += 1;
    }
  }

  private void addShapeIndexInternal(S2ShapeIndex index, boolean complement) {
    S2ContainsPointQuery query = index.isFresh() ? new S2ContainsPointQuery(index) : null;
    List<S2Shape> shapes = index.getShapes();
    for (int shapeId = 0; shapeId < shapes.size(); ++shapeId) {
      S2Shape shape = shapes.get(shapeId);
      if (shape.dimension() != 2) {
        throw new IllegalArgumentException(
            "Only 2D shapes are supported in S2MultiOperation. Found shape with dimension: "
                + shape.dimension());
      }
      num2dShapes++;
      boolean shapeContainsOrigin =
          (query != null) ? query.shapeContains(shapeId, S2.origin()) : shape.containsOrigin();
      if (complement ^ shapeContainsOrigin) {
        refWinding += 1;
      }
      for (List<S2Point> chain : shape.chains()) {
        windingOp.addLoop(complement ? Lists.reverse(chain) : chain);
      }
    }
  }

  /** Builds the result. */
  public void build() {
    S2Error error = new S2Error();
    WindingRule rule = WindingRule.POSITIVE;
    int effectiveRefWinding = refWinding;

    switch (opType) {
      case UNION -> rule = WindingRule.POSITIVE;
      case INTERSECTION -> {
        effectiveRefWinding = refWinding - num2dShapes + 1;
        rule = WindingRule.POSITIVE;
      }
      case SYMMETRIC_DIFFERENCE -> rule = WindingRule.ODD;
    }

    boolean success = windingOp.build(S2.origin(), effectiveRefWinding, rule, error);
    checkState(success, "Failed to build multi-operation: %s", error.text());
  }

  // Static sugar for common cases

  private static S2Polygon booleanOp(
      OpType opType, S2WindingOperation.Options options, Iterable<S2Polygon> polygons) {
    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(opType, layer, options);
    for (S2Polygon p : polygons) {
      op.add(p);
    }
    op.build();
    return result;
  }

  /**
   * Computes the union of the given polygons with default options ({@link #DEFAULT_SNAP_FUNCTION}).
   *
   * <p>Returns an empty polygon for 0 inputs, or a copy of the single input polygon directly
   * (without running {@link S2Builder} snapping) for 1 input. For 2 or more polygons, inputs are
   * pre-sorted by {@link S2Polygon#compareTo} to guarantee permutation-invariant output.
   *
   * @param polygons The polygons to union.
   * @return A new polygon containing the union of the input polygons.
   */
  public static S2Polygon union(S2Polygon... polygons) {
    return union(Arrays.asList(polygons), defaultOptions());
  }

  /**
   * Computes the union of the given polygons with default options ({@link #DEFAULT_SNAP_FUNCTION}).
   *
   * <p>Returns an empty polygon for 0 inputs, or a copy of the single input polygon directly
   * (without running {@link S2Builder} snapping) for 1 input. For 2 or more polygons, inputs are
   * pre-sorted by {@link S2Polygon#compareTo} to guarantee permutation-invariant output.
   *
   * @param polygons The polygons to union.
   * @return A new polygon containing the union of the input polygons.
   */
  public static S2Polygon union(Iterable<S2Polygon> polygons) {
    return union(polygons, defaultOptions());
  }

  /**
   * Computes the union of the given polygons with the specified options.
   *
   * <p>Returns an empty polygon for 0 inputs, or a copy of the single input polygon directly
   * (without running {@link S2Builder} snapping) for 1 input. For 2 or more polygons, inputs are
   * pre-sorted by {@link S2Polygon#compareTo} to guarantee permutation-invariant output.
   *
   * @param polygons The polygons to union.
   * @param options Options for the underlying {@link S2WindingOperation}.
   * @return A new polygon containing the union of the input polygons.
   */
  public static S2Polygon union(Iterable<S2Polygon> polygons, S2WindingOperation.Options options) {
    ImmutableList<S2Polygon> sorted = ImmutableList.sortedCopyOf(polygons);
    if (sorted.isEmpty()) {
      return new S2Polygon();
    } else if (sorted.size() == 1) {
      return new S2Polygon(sorted.get(0));
    }
    return booleanOp(OpType.UNION, options, sorted);
  }

  /**
   * Computes the intersection of the given polygons with default options ({@link
   * #DEFAULT_SNAP_FUNCTION}).
   *
   * <p>Returns a full polygon for 0 inputs, or a copy of the single input polygon directly (without
   * running {@link S2Builder} snapping) for 1 input. For 2 or more polygons, inputs are pre-sorted
   * by {@link S2Polygon#compareTo} to guarantee permutation-invariant output.
   *
   * @param polygons The polygons to intersect.
   * @return A new polygon containing the intersection of the input polygons.
   */
  public static S2Polygon intersection(S2Polygon... polygons) {
    return intersection(Arrays.asList(polygons));
  }

  /**
   * Computes the intersection of the given polygons with default options ({@link
   * #DEFAULT_SNAP_FUNCTION}).
   *
   * <p>Returns a full polygon for 0 inputs, or a copy of the single input polygon directly (without
   * running {@link S2Builder} snapping) for 1 input. For 2 or more polygons, inputs are pre-sorted
   * by {@link S2Polygon#compareTo} to guarantee permutation-invariant output.
   *
   * @param polygons The polygons to intersect.
   * @return A new polygon containing the intersection of the input polygons.
   */
  public static S2Polygon intersection(Iterable<S2Polygon> polygons) {
    ImmutableList<S2Polygon> sorted = ImmutableList.sortedCopyOf(polygons);
    if (sorted.isEmpty()) {
      return new S2Polygon(S2Loop.full());
    } else if (sorted.size() == 1) {
      return new S2Polygon(sorted.get(0));
    }
    return booleanOp(OpType.INTERSECTION, defaultOptions(), sorted);
  }

  /**
   * Computes the difference: regions covered by the minuend but not by any of the subtrahends,
   * using default options ({@link #DEFAULT_SNAP_FUNCTION}).
   *
   * <p>Returns a copy of {@code minuend} directly (without running {@link S2Builder} snapping) if
   * {@code minuend} is empty or {@code subtrahends} is empty.
   *
   * @param minuend The geometry from which to subtract.
   * @param subtrahends The geometries to subtract from the minuend.
   * @return A new polygon containing the difference.
   */
  public static S2Polygon difference(S2Polygon minuend, List<S2Polygon> subtrahends) {
    if (minuend.isEmpty() || subtrahends.isEmpty()) {
      return new S2Polygon(minuend);
    }
    S2Polygon result = new S2Polygon();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(result);
    S2MultiOperation op = new S2MultiOperation(OpType.INTERSECTION, layer, defaultOptions());
    op.add(minuend);
    ImmutableList.sortedCopyOf(subtrahends).forEach(op::addComplement);
    op.build();
    return result;
  }
}
