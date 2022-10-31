/*
 * Copyright 2021 Google Inc.
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

/**
 * S2BestDistanceTarget is an interface representing an abstract geometric object to which a "best
 * distance" is being computed. The meaning of "best" is determined by the implementation, but
 * minimum or maximum are typical.
 */
@CheckReturnValue
public interface S2BestDistanceTarget<D extends S1Distance<D>> {
  /**
   * For every connected component of this target, a point on that component is passed to the
   * visitor. For example, if this target consists of points "a" and "b", the visitor will be called
   * twice, once with each point. If the target consists of the edge "ab", the visitor will be
   * called once, with a point on the edge. If the visitor ever returns false, visiting connected
   * components stops and this method returns false. Otherwise, all connected components will be
   * visited and this method returns true.
   */
  @CanIgnoreReturnValue
  public boolean visitConnectedComponentPoints(PointVisitor visitor);

  /**
   * If includeInteriors() is true, distance will be measured to the boundary and interior of
   * polygons in this target, rather than to polygon boundaries only. This is only meaningful for
   * targets with interiors, i.e. polygons, which are only present in S2ShapeIndexTargets. Note that
   * S2 cell target interiors are always included.
   */
  public default void setIncludeInteriors(boolean includeInteriors) {}

  /**
   * The default implementation returns false, as most implementations don't have interiors. Note
   * that S2 cell target interiors are always included.
   */
  public default boolean includeInteriors() {
    return false;
  }

  /**
   * Returns an S2Cap bounding the region of the sphere that has the best possible distance to this
   * target.
   */
  public S2Cap getCapBound();

  /**
   * If the distance to the point p from this target is better than "bestDist", then updates
   * "bestDist" and returns true. Otherwise returns false.
   */
  @CanIgnoreReturnValue
  public boolean updateBestDistance(S2Point p, DistanceCollector<D> bestDist);

  /**
   * If the distance to the edge (v0, v1) from this target is better than "bestDist", then updates
   * "bestDist" and returns true. Otherwise returns false.
   */
  @CanIgnoreReturnValue
  public boolean updateBestDistance(S2Point v0, S2Point v1, DistanceCollector<D> bestDist);

  /**
   * If the distance to the given S2Cell (including its interior) from this target is better than
   * "bestDist", then updates "bestDist" and returns true. Otherwise returns false.
   */
  @CanIgnoreReturnValue
  public boolean updateBestDistance(S2Cell cell, DistanceCollector<D> bestDist);

  /**
   * Specifies that whenever one of the updateBestDistance() methods above returns "true", the
   * returned distance is allowed to be up to "maxError" worse than the true best distances. In
   * other words, it gives this target object permission to terminate its distance calculation as
   * soon as it has determined that:
   *
   * <p>(1) the best distance is better than bestDist and
   *
   * <p>(2) the best possible further improvement is less than "maxError".
   *
   * <p>If the target takes advantage of "maxError" to optimize its distance calculation, this
   * method must return "true". (Most target types can use the default implementation which simply
   * returns false.)
   */
  public default boolean setMaxError(D maxError) {
    return false;
  }

  /**
   * The following method is provided as a convenience for classes that compute distances to a
   * collection of indexed geometry, such as S2ClosestPointQuery, S2ClosestEdgeQuery, and
   * S2ClosestCellQuery. It returns the maximum number of indexed edges for which it is faster to
   * compute the distance by brute force (e.g., by testing every edge) rather than by using an
   * index. The appropriate value is different for each target type and index type and can be
   * estimated for a given (distance target, index type) pair by running benchmarks.
   *
   * <p>By default this method returns -1, indicating that it is not implemented.
   */
  public default int maxBruteForceIndexSize() {
    return -1;
  }
}
