/*
 * Copyright 2018 Google Inc.
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
import com.google.common.geometry.S2ContainsPointQuery.S2VertexModel;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import com.google.common.geometry.primitives.IntVector;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.Collection;

/**
 * This class wraps an S2ShapeIndex object with the additional methods needed to implement the
 * S2Region API, in order to allow S2RegionCoverer to compute S2CellId coverings of arbitrary
 * collections of geometry.
 *
 * <p>It also contains a method visitIntersectingShapes() that may be used to efficiently visit all
 * shapes that intersect an arbitrary S2CellId (not limited to cells in the index).
 *
 * <p>Example usage:
 *
 * <pre>
 * S2CellUnion getCovering(S2ShapeIndex index) {
 *   S2RegionCoverer coverer = new S2RegionCoverer();
 *   coverer.setMaxCells(20);
 *   return coverer.getCovering(new S2ShapeIndexRegion(index));
 * }
 * </pre>
 *
 * <p>This class uses a number of temporary mutable objects to keep allocation down, and so is not
 * thread-safe. To use it in parallel, each thread should construct its own instance (this is not
 * expensive).
 */
public class S2ShapeIndexRegion implements S2Region {
  /** The query for contains(S2Point) tests. */
  private final S2ContainsPointQuery containsQuery;

  /** The current index wrapped by this region. */
  private final S2ShapeIndex index;

  /** The iterator over the cells of the underlying S2ShapeIndex. */
  private final S2Iterator<S2ShapeIndex.Cell> it;

  /** Temporary cell union for internal usage. */
  private final S2CellUnion union = new S2CellUnion();

  /** Temporary edge for internal usage. */
  private final MutableEdge edge = new MutableEdge();

  /** Temporary bound for internal usage. */
  private final R2Rect bound = new R2Rect();

  /** Temporary R2 point for internal usage. */
  private final R2Vector p0 = new R2Vector();

  /** Temporary R2 point for internal usage. */
  private final R2Vector p1 = new R2Vector();

  /** Temporary shape IDs used in query processing. */
  private final IntVector ids = new IntVector();

  /**
   * Creates a new region with the given index, and a {@link S2VertexModel#SEMI_OPEN semi-open}
   * vertex model.
   */
  public S2ShapeIndexRegion(S2ShapeIndex index) {
    this(index, S2VertexModel.SEMI_OPEN);
  }

  /** Creates a new region with the given index, and a given {@link S2VertexModel}. */
  public S2ShapeIndexRegion(S2ShapeIndex index, S2VertexModel model) {
    this.index = index;
    this.it = index.iterator();
    this.containsQuery = new S2ContainsPointQuery(index, new S2ContainsPointQuery.Options(model));
  }

  @Override
  public S2Cap getCapBound() {
    getCellUnionBound(union.cellIds());
    return union.getCapBound();
  }

  @Override
  public S2LatLngRect getRectBound() {
    getCellUnionBound(union.cellIds());
    return union.getRectBound();
  }

  /**
   * Clears the given list of cells and adds the cell union of this index. An index of shapes in one
   * face adds up to 4 cells, otherwise up to 6 may be added.
   */
  @Override
  public void getCellUnionBound(Collection<S2CellId> cellIds) {
    // We find the range of S2Cells spanned by the index and choose a level such that the entire
    // index can be covered with just a few cells. There are two cases:
    //
    // - If the index intersects two or more faces, then for each intersected face we add one cell
    //   to the covering. Rather than adding the entire face, instead we add the smallest S2Cell
    //   that covers the S2ShapeIndex cells within that face.
    //
    // - If the index intersects only one face, then we first find the smallest cell S that contains
    //   the index cells (just like the case above). However rather than using the cell S itself,
    //   instead we repeat this process for each of its child cells. In other words, for each child
    //   cell C we add the smallest S2Cell C' that covers the index cells within C. This extra step
    //   is relatively cheap and produces much tighter coverings when the S2ShapeIndex consists of a
    //   small region near the center of a large S2Cell.
    //
    // The following code uses only a single S2Iterator object because creating an S2Iterator may be
    // relatively expensive for S2ShapeIndex instances (e.g. it may involve substantial memory
    // allocation to build a lazily-assembled index).
    cellIds.clear();

    // Find the last S2CellId in the index.
    it.finish();
    if (!it.prev()) {
      // Empty index.
      return;
    }
    S2CellId lastIndexId = it.id();
    it.restart();
    S2CellId currentIndexId = it.id();
    if (!currentIndexId.equals(lastIndexId)) {
      // The index has at least two cells. Choose an S2CellId level such that the entire index can
      // be spanned with at most 6 cells (if the index spans multiple faces) or 4 cells (if the
      // index spans a single face).
      int level = currentIndexId.getCommonAncestorLevel(lastIndexId) + 1;

      // For each cell C at the chosen level, we compute the smallest S2Cell that covers the
      // S2ShapeIndex cells within C.
      S2CellId lastId = lastIndexId.parent(level);
      for (S2CellId id = currentIndexId.parent(level); !id.equals(lastId); id = id.next()) {
        // If the cell C does not contain any index cells, then skip it.
        S2CellId max = id.rangeMax();
        if (max.lessThan(currentIndexId)) {
          continue;
        }

        // Find the range of index cells contained by C and then shrink C so that it just covers
        // those cells.
        it.seek(max.next());
        it.prev();
        coverRange(currentIndexId, it.id(), cellIds);
        it.next();
        currentIndexId = it.id();
      }
    }
    coverRange(currentIndexId, lastIndexId, cellIds);
  }

  /**
   * Computes the smallest S2Cell that covers the S2Cell range (first, last) and adds this cell to
   * "cellIds".
   *
   * @throws IllegalArgumentException "first" and "last" don't have a common ancestor.
   */
  private static void coverRange(S2CellId first, S2CellId last, Collection<S2CellId> cellIds) {
    if (first.equals(last)) {
      // The range consists of a single index cell.
      cellIds.add(first);
    } else {
      // Add the lowest common ancestor of the given range.
      int level = first.getCommonAncestorLevel(last);
      Preconditions.checkArgument(level >= 0, "First and last must have a common ancestor.");
      cellIds.add(first.parent(level));
    }
  }

  /**
   * Returns true if the given point is contained by any two-dimensional shape (i.e., polygon). Zero
   * and one-dimensional shapes are ignored by this method.
   */
  @Override
  public boolean contains(S2Point p) {
    if (it.locate(p)) {
      S2Point center = it.center();
      S2ShapeIndex.Cell cell = it.entry();
      for (int s = 0; s < cell.numShapes(); ++s) {
        if (containsQuery.shapeContains(center, cell.clipped(s), p)) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Returns true if 'target' is contained by any single shape. If the cell is covered by a union of
   * different shapes then it may return false.
   *
   * <p>This implementation is conservative but not exact; if a shape just barely contains the given
   * cell then it may return false. The maximum error is less than 10 * DBL_EPSILON radians (or
   * about 15 nanometers).
   */
  @Override
  public boolean contains(S2Cell target) {
    S2ShapeIndex.CellRelation relation = it.locate(target.id());

    // If the relation is DISJOINT, then "target" is not contained. Similarly if the relation is
    // SUBDIVIDED then "target" is not contained, since index cells are subdivided only if they
    // (nearly) intersect too many edges.
    if (relation != S2ShapeIndex.CellRelation.INDEXED) {
      return false;
    }

    // Otherwise, the iterator points to an index cell containing "target". If any shape contains
    // the target cell, we return true.
    // assert (it.id().contains(target.id()));
    S2ShapeIndex.Cell cell = it.entry();
    S2Point center = it.center();
    for (int s = 0; s < cell.numShapes(); ++s) {
      S2ClippedShape clipped = cell.clipped(s);
      // The shape contains the target cell iff the shape contains the cell center and none of its
      // edges intersects the (padded) cell interior.
      if (it.id().equals(target.id())) {
        if (clipped.numEdges() == 0 && clipped.containsCenter()) {
          return true;
        }
      } else {
        // It is faster to call AnyEdgeIntersects() before Contains().
        S2Shape shape = index.getShapes().get(clipped.shapeId());
        if (shape.hasInterior()
            && !anyEdgeIntersects(clipped, target)
            && containsQuery.shapeContains(center, clipped, target.getCenter())) {
          return true;
        }
      }
    }
    return false;
  }

  /**
   * A ShapeVisitor is a function that is called with shapes that intersects a target S2Cell.
   *
   * <p>{@link ShapeVisitor#test} will be called for each shape, where containsTarget means that the
   * shape fully contains the target S2Cell. The function should return true to continue visiting
   * intersecting shapes, or false to terminate the algorithm early.
   */
  public interface ShapeVisitor {
    boolean test(int shapeId, boolean containsTarget);
  }

  /**
   * Visits all shapes that intersect "target", terminating early if the "visitor" returns false, in
   * which case visitIntersectingShapes returns false as well. Each shape is visited at most once.
   * This method can also be used to visit all shapes that fully contain "target"
   * (visitContainingShapes) by simply having the ShapeVisitor predicate immediately return true
   * when containsTarget is false.
   */
  @CanIgnoreReturnValue
  public boolean visitIntersectingShapes(S2Cell target, ShapeVisitor visitor) {
    S2ShapeIndex.CellRelation relation = it.locate(target.id());
    switch (relation) {
      case DISJOINT:
        return true;

      case SUBDIVIDED:
        // A shape contains the target cell iff (1) it appears in at least one descendent cell,
        // and (2) it contains the center of all descendent cells, and (3) it has no edges in any
        // descendent cell. Add the shapeId<<1|contains bits, sortDistinct the bits, and then since
        // each shape may occur twice, once with each contains bit in that order, and we want to
        // report !contains if any id was not contained, then simply skip any cases where the shape
        // ID is the same as the last shape ID, since the preceding contains bit must be 0.
        ids.clear();
        for (S2CellId max = target.id().rangeMax();
            !it.done() && it.id().lessOrEquals(max);
            it.next()) {
          S2ShapeIndex.Cell cell = it.entry();
          for (int s = 0; s < cell.numShapes(); ++s) {
            S2ClippedShape clipped = cell.clipped(s);
            boolean containsTarget = clipped.numEdges() == 0 && clipped.containsCenter();
            ids.add((clipped.shapeId() << 1) | (containsTarget ? 1 : 0));
          }
        }
        ids.sort();
        ids.unique();

        for (int i = 0, last = -1; i < ids.size(); i++) {
          int id = ids.get(i);
          boolean contains = (id & 1) != 0;
          int shapeId = id >>> 1;
          if (last != shapeId && !visitor.test(shapeId, contains)) {
            return false;
          }
          last = shapeId;
        }
        return true;

      case INDEXED:
        S2ShapeIndex.Cell cell = it.entry();
        for (int s = 0; s < cell.numShapes(); ++s) {
          // The shape contains the target cell iff the shape contains the cell center and none of
          // its edges intersect the (padded) cell interior.
          S2ClippedShape clipped = cell.clipped(s);
          boolean contains = false;
          if (cell.id() == target.id().id()) {
            // The index contains the target exactly.
            contains = clipped.numEdges() == 0 && clipped.containsCenter();
          } else {
            // The index contains a parent of the target.
            if (!anyEdgeIntersects(clipped, target)) {
              if (!containsQuery.shapeContains(it.center(), clipped, target.getCenter())) {
                // The shape and target are disjoint.
                continue;
              }
              // The shape contains the target center, and no edges intersect.
              contains = true;
            }
          }
          if (!visitor.test(clipped.shapeId(), contains)) {
            return false;
          }
        }
        return true;

      default:
        throw new IllegalStateException("Unknown S2ShapeIndex.CellRelation " + relation);
    }
  }

  /**
   * Returns true if any shape intersects "target".
   *
   * <p>This implementation is conservative but not exact; if a shape is just barely disjoint from
   * the given cell then it may return true. The maximum error is less than 10 * DBL_EPSILON radians
   * (or about 15 nanometers).
   */
  @Override
  public boolean mayIntersect(S2Cell target) {
    S2ShapeIndex.CellRelation relation = it.locate(target.id());

    // If "target" does not overlap any index cell, there is no intersection.
    if (relation == S2ShapeIndex.CellRelation.DISJOINT) {
      return false;
    }

    // If "target" is subdivided into one or more index cells, then there is an intersection to
    // within the S2ShapeIndex error bound.
    if (relation == S2ShapeIndex.CellRelation.SUBDIVIDED) {
      return true;
    }

    // Otherwise, the iterator points to an index cell containing "target". If "target" is an index
    // cell itself, there is an intersection because index cells are created only if they have at
    // least one edge or they are entirely contained by the loop.
    // assert (it.id().contains(target.id()));
    if (it.compareTo(target.id()) == 0) {
      return true;
    }

    // Test whether any shape intersects the target cell or contains its center.
    S2ShapeIndex.Cell cell = it.entry();
    S2Point center = it.center();
    for (int s = 0; s < cell.numShapes(); ++s) {
      S2ClippedShape clipped = cell.clipped(s);
      if (anyEdgeIntersects(clipped, target)
          || containsQuery.shapeContains(center, clipped, target.getCenter())) {
        return true;
      }
    }
    return false;
  }

  /**
   * Returns true if any edge of the indexed shape "clipped" intersects the cell "target". It may
   * also return true if an edge is very close to "target"; the maximum error is less than 10 *
   * DBL_EPSILON radians (about 15 nanometers).
   */
  private boolean anyEdgeIntersects(S2ClippedShape clipped, S2Cell target) {
    target.setBoundUV(bound);
    bound.expand(S2EdgeUtil.MAX_CELL_EDGE_ERROR);
    int face = target.face();
    S2Shape shape = index.getShapes().get(clipped.shapeId());
    int numEdges = clipped.numEdges();
    for (int i = 0; i < numEdges; ++i) {
      shape.getEdge(clipped.edge(i), edge);
      if (S2EdgeUtil.clipToPaddedFace(edge.a, edge.b, face, S2EdgeUtil.MAX_CELL_EDGE_ERROR, p0, p1)
          && S2EdgeUtil.intersectsRect(p0, p1, bound)) {
        return true;
      }
    }
    return false;
  }
}
