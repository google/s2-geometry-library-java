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

import com.google.common.geometry.S2ShapeUtil.PointVisitor;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import com.google.errorprone.annotations.CheckReturnValue;

/** Targets for finding minimum distances. */
@CheckReturnValue
class S2MinDistanceTargets {
  private S2MinDistanceTargets() {}

  /** PointTarget specializes S2BestEdgesQueryBase.PointTarget for minimum distances. */
  public static class PointTarget<D2 extends S1Distance<D2>>
      extends S2BestEdgesQueryBase.PointTarget<D2> implements S2BestDistanceTarget<D2> {
    public PointTarget(S2Point p) {
      super(p);
    }

    @Override
    public S2Cap getCapBound() {
      return S2Cap.fromAxisChord(point, S1ChordAngle.ZERO);
    }
  }

  /** EdgeTarget specializes S2BestEdgesQueryBase.EdgeTarget for minimum distances. */
  public static class EdgeTarget<D2 extends S1Distance<D2>>
      extends S2BestEdgesQueryBase.EdgeTarget<D2> implements S2BestDistanceTarget<D2> {
    public EdgeTarget(S2Point a, S2Point b) {
      super(a, b);
    }

    @Override
    public S2Cap getCapBound() {
      double r2 = getHalfEdgeLength2();
      return S2Cap.fromAxisChord(a.add(b).normalize(), S1ChordAngle.fromLength2(r2));
    }
  }

  /** CellTarget specializes S2BestEdgesQueryBase.CellTarget for minimum distances. */
  public static class CellTarget<D2 extends S1Distance<D2>>
      extends S2BestEdgesQueryBase.CellTarget<D2> implements S2BestDistanceTarget<D2> {
    public CellTarget(S2Cell c) {
      super(c);
    }

    @Override
    public S2Cap getCapBound() {
      return cell.getCapBound();
    }
  }

  /** S2ShapeIndexTarget is used for finding minimum distances to shape indexes. */
  public static class ShapeIndexTarget<D extends S1Distance<D>> implements S2BestDistanceTarget<D> {
    private final S2ShapeIndex index;
    private final S2BestEdgesQueryBase.Builder<D> queryBuilder;
    private S2BestEdgesQueryBase<D> minDistanceQuery = null;

    /**
     * Clients using S1ChordAngle as their S1Distance type may find it convenient to use
     * {@link S2ClosestEdgeQuery#createShapeIndexTarget(S2ShapeIndex)}.
     *
     * <p>Otherwise, constructing a ShapeIndexTarget for a specific S1Distance type requires
     * providing a S2BestEdgesQueryBase.Builder for closest edges and the templated S1Distance
     * type.
     */
    public ShapeIndexTarget(S2ShapeIndex index, S2BestEdgesQueryBase.Builder<D> queryBuilder) {
      this.index = index;
      this.queryBuilder = queryBuilder;
    }

    /**
     * Note that changing the maxError option after a ShapeIndexTarget has been used requires
     * rebuilding the internal query.
     */
    @Override
    public boolean setMaxError(D maxError) {
      queryBuilder.setMaxError(maxError);
      minDistanceQuery = null;
      return true;
    }

    /**
     * Note that changing the maxError option after a ShapeIndexTarget has been used requires
     * rebuilding the internal query state.
     */
    @Override
    public void setIncludeInteriors(boolean includeInteriors) {
      queryBuilder.setIncludeInteriors(includeInteriors);
      minDistanceQuery = null;
    }

    @Override
    public boolean includeInteriors() {
      return queryBuilder.includeInteriors();
    }

    @Override
    public S2Cap getCapBound() {
      return new S2ShapeIndexRegion(index).getCapBound();
    }

    /**
     * If the distance to the given Target from this ShapeIndexTarget's indexed geometry is less
     * than minDist, update minDist and return true. Otherwise return false.
     */
    private boolean updateMinDistance(
        S2BestDistanceTarget<D> target, DistanceCollector<D> minDist) {
      if (minDistanceQuery == null) {
        minDistanceQuery = queryBuilder.build(index);
      }

      if (minDistanceQuery.atBestLimit(minDist)) {
        return false;  // No better result is possible.
      }
      // Is there a result closer than the current minDist?
      return minDistanceQuery.updateBestDistance(target, minDist);
    }

    @Override
    @CanIgnoreReturnValue
    public boolean updateBestDistance(S2Point p, DistanceCollector<D> minDist) {
      return updateMinDistance(new PointTarget<D>(p), minDist);
    }

    @Override
    @CanIgnoreReturnValue
    public boolean updateBestDistance(S2Point v0, S2Point v1, DistanceCollector<D> minDist) {
      return updateMinDistance(new EdgeTarget<D>(v0, v1), minDist);
    }

    @Override
    @CanIgnoreReturnValue
    public boolean updateBestDistance(S2Cell cell, DistanceCollector<D> minDist) {
      return updateMinDistance(new CellTarget<D>(cell), minDist);
    }

    @CanIgnoreReturnValue
    @Override
    public boolean visitConnectedComponentPoints(PointVisitor visitor) {
      for (S2Shape shape : index.getShapes()) {
        if (shape == null) {
          continue;
        }
        // Shapes that don't have any edges require a special case (below).
        boolean testedPoint = false;
        for (int chain = 0; chain < shape.numChains(); ++chain) {
          if (shape.getChainLength(chain) == 0) {
            continue;
          }
          // Visit the first vertex of the shape chain.
          testedPoint = true;
          if (!visitor.apply(shape.getChainVertex(chain, 0))) {
            return false;
          }
        }

        if (!testedPoint) {
          // Special case to handle full polygons.
          S2Shape.ReferencePoint ref = shape.getReferencePoint();
          if (!ref.contained()) {
            continue;
          }
          // Visit the reference point for the full polygon.
          if (!visitor.apply(ref)) {
            return false;
          }
        }
      }
      return true;
    }
  }

}
