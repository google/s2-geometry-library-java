/*
 * Copyright 2014 Google Inc.
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

import java.io.Serializable;
import java.util.Arrays;
import java.util.Collection;

/**
 * An S2RegionIntersection represents an intersection of overlapping regions. It is convenient for
 * computing a covering of the intersection of a set of regions. The regions are assumed to be
 * immutable. Note: An intersection of no regions covers the entire sphere.
 */
public class S2RegionIntersection implements S2Region, Serializable {
  // Regions is non-private so that it can be accessed from the custom field serializer.
  final S2Region[] regions;
  private transient S2LatLngRect cachedRectBound = null;

  /** Create an intersection from a copy of {@code regions}. */
  public S2RegionIntersection(Collection<S2Region> regions) {
    this.regions = regions.toArray(new S2Region[0]);
  }

  /** Returns true if all the regions fully contain the cell. */
  @Override
  public boolean contains(S2Cell cell) {
    for (S2Region region : regions) {
      if (!region.contains(cell)) {
        return false;
      }
    }

    return true;
  }

  /** Returns true if all the regions fully contain the point. */
  @Override
  public boolean contains(S2Point point) {
    for (S2Region region : regions) {
      if (!region.contains(point)) {
        return false;
      }
    }

    return true;
  }

  @Override
  public S2Cap getCapBound() {
    // This could be optimized to return a tighter bound, but doesn't seem worth it unless
    // profiling shows otherwise.
    return getRectBound().getCapBound();
  }

  @Override
  public S2LatLngRect getRectBound() {
    if (cachedRectBound != null) {
      return cachedRectBound;
    }

    S2LatLngRect.Builder builder = new S2LatLngRect.Builder(S2LatLngRect.full());
    for (S2Region region : regions) {
      builder.intersection(region.getRectBound());
    }
    cachedRectBound = builder.build();
    return cachedRectBound;
  }

  /** Returns true if the cell may intersect all regions in this collection. */
  @Override
  public boolean mayIntersect(S2Cell cell) {
    for (S2Region region : regions) {
      if (!region.mayIntersect(cell)) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns true if this S2RegionIntersection is equal to another S2RegionIntersection, where each
   * region must be equal and in the same order. This method is intended only for testing purposes.
   * NOTE: This should be rewritten to disregard order if such functionality is ever required.
   */
  @Override
  public boolean equals(Object thatObject) {
    if (!(thatObject instanceof S2RegionIntersection)) {
      return false;
    }
    S2RegionIntersection that = (S2RegionIntersection) thatObject;
    return Arrays.equals(regions, that.regions);
  }

  @Override
  public int hashCode() {
    return Arrays.hashCode(regions);
  }
}
