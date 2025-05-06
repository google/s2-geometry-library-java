/*
 * Copyright 2022 Google Inc.
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

import com.google.common.geometry.S2ClosestEdgeQuery.CellTarget;
import com.google.common.geometry.S2ClosestEdgeQuery.PointTarget;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import jsinterop.annotations.JsMethod;

/**
 * This class provides a way to expand an arbitrary collection of geometry by a fixed radius (an
 * operation variously known as "buffering", "offsetting", or "Minkowski sum with a disc") in order
 * to compute an S2CellId covering (see S2RegionCoverer). The resulting covering contains all points
 * within the given radius of any point in the original geometry.
 *
 * <p>This class does not actually buffer the geometry; instead it implements the S2Region API by
 * computing the distance from candidate S2CellIds to the original geometry. If this distance is
 * below the given radius then the S2CellId intersects the buffered geometry. For example, if the
 * original geometry consists of a single S2Point then the buffered geometry is exactly equivalent
 * to an S2Cap with the given radius. (Note that the region is not approximated as a polygonal
 * loop.)
 *
 * <p>Example usage:
 *
 * {@snippet :
 * S2CellUnion getBufferedCovering(S2ShapeIndex index, S1Angle radius) {
 *   S2ShapeIndexBufferedRegion region = new S2ShapeIndexBufferedRegion(index, radius);
 *   S2CellUnion covering = S2RegionCoverer.DEFAULT.getCovering(region);
 *   return covering;
 * }
 *
 * }
 */
public class S2ShapeIndexBufferedRegion implements S2Region {
  /** The radius to expand the index geometry by. */
  private final S1ChordAngle radius;

  /** Used to determine if the distance to points or cells is within the specified "radius". */
  private final S2ClosestEdgeQuery<S1ChordAngle> query;

  /** This shapeIndexRegion wraps the given S2ShapeIndex with an S2Region API. */
  private S2ShapeIndexRegion shapeIndexRegion = null;

  /**
   * To handle (radius == 0) correctly, we need to test whether distances are less than or equal to
   * "radius". This is done by testing whether distances are less than radius.successor().
   */
  private final S1ChordAngle radiusSuccessor;

  /**
   * Constructs a {@link S2Region} representing all points within the given radius of any point in
   * the given S2ShapeIndex.
   */
  public S2ShapeIndexBufferedRegion(S2ShapeIndex index, S1ChordAngle radius) {
    this.radius = radius;
    radiusSuccessor = radius.successor();
    query = S2ClosestEdgeQuery.builder().setIncludeInteriors(true).build(index);
  }

  /** Convenience constructor that accepts an S1Angle for the radius. */
  public S2ShapeIndexBufferedRegion(S2ShapeIndex index, S1Angle radius) {
    this(index, S1ChordAngle.fromS1Angle(radius));
  }

  /** Returns the underlying S2ShapeIndex that was supplied to init() or the constructor. */
  public S2ShapeIndex index() {
    return query.index();
  }

  /** Returns the buffering radius that was supplied to init() or the constructor. */
  public S1ChordAngle radius() {
    return radius;
  }

  /**
   * Used to defer construction of the S2ShapeIndexRegion until the first time it is used, as doing
   * so will trigger building the index if not already built. Then returns the cached instance.
   */
  private S2ShapeIndexRegion makeS2ShapeIndexRegion() {
    if (shapeIndexRegion == null) {
      shapeIndexRegion = new S2ShapeIndexRegion(index());
    }
    return shapeIndexRegion;
  }

  @Override
  public S2Cap getCapBound() {
    S2Cap origCap = makeS2ShapeIndexRegion().getCapBound();
    if (origCap.isEmpty()) {
      return origCap;
    }
    return S2Cap.fromAxisChord(origCap.axis(), S1ChordAngle.add(origCap.radius(), radius));
  }

  @Override
  public S2LatLngRect getRectBound() {
    S2LatLngRect origRect = makeS2ShapeIndexRegion().getRectBound();
    return origRect.expandedByDistance(radius.toAngle());
  }

  @Override
  public void getCellUnionBound(Collection<S2CellId> results) {
    // Get the max level to add neighbors at, or return the full union if the buffer radius is huge.
    double radians = radius.toAngle().radians();
    int maxLevel = S2Projections.MIN_WIDTH.getMaxLevel(radians) - 1;
    if (maxLevel < 0) {
      S2Cap.full().getCellUnionBound(results);
      return;
    }

    // We start with a covering of the original S2ShapeIndex, and then expand it by replacing each
    // cell with a block of 4 cells whose union contains the original cell buffered by the radius.
    //
    // This increases the number of cells in the covering by a factor of 4 and increases the covered
    // area by a factor of 16, so it is not a very good covering, but it is much better than always
    // returning the 6 face cells.
    List<S2CellId> origCellIds = new ArrayList<>(6);
    new S2ShapeIndexRegion(index()).getCellUnionBound(origCellIds);
    results.clear();
    for (S2CellId id : origCellIds) {
      if (id.isFace()) {
        S2Cap.full().getCellUnionBound(results);
        return;
      }
      id.getVertexNeighbors(min(maxLevel, id.level() - 1), results);
    }
  }

  @Override
  @JsMethod(name = "containsCell")
  public boolean contains(S2Cell cell) {
    // Return true if the buffered region is guaranteed to cover the whole globe.
    if (radiusSuccessor.greaterThan(S1ChordAngle.STRAIGHT)) {
      return true;
    }

    // To implement this method perfectly requires using the directed Hausdorff distance, which is
    // expensive. However the following heuristic is almost as good in practice and much cheaper to
    // compute.

    // Return true if the unbuffered region contains this cell.
    if (makeS2ShapeIndexRegion().contains(cell)) {
      return true;
    }

    // Otherwise approximate the cell by its bounding cap.
    //
    // NOTE(ericv): It would be slightly more accurate to first find the closest point in the
    // indexed geometry to the cell, and then measure the actual maximum distance from that point to
    // the cell (a poor man's Hausdorff distance). But based on actual tests this is not worthwhile.
    S2Cap cap = cell.getCapBound();
    if (radius.lessThan(cap.radius())) {
      return false;
    }

    // Return true if the distance to the cell center plus the radius of the cell's bounding cap is
    // less than or equal to "radius".
    PointTarget<S1ChordAngle> target = new PointTarget<>(cell.getCenter());
    return query.isDistanceLess(target, S1ChordAngle.sub(radiusSuccessor, cap.radius()));
  }

  @Override
  @JsMethod(name = "containsPoint")
  public boolean contains(S2Point p) {
    // Returns true if the distance from the point to the target is less than or equal to "radius".
    return query.isDistanceLess(new PointTarget<>(p), radiusSuccessor);
  }

  @Override
  public boolean mayIntersect(S2Cell cell) {
    // Returns true if the distance from the point to the target is less than or equal to "radius".
    return query.isDistanceLess(new CellTarget<>(cell), radiusSuccessor);
  }
}
