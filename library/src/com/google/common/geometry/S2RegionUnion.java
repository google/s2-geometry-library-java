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
 * An S2RegionUnion represents a union of possibly overlapping regions. It is convenient for
 * computing a covering of a set of regions. The regions are assumed to be immutable.
 *
 * <p>However, note that currently, using S2RegionCoverer to compute coverings of S2RegionUnions may
 * produce coverings with considerably less than the requested number of cells in cases of
 * overlapping or tiling regions. This occurs because the current S2RegionUnion.Contains
 * implementation for S2Cells only returns true if the cell is fully contained by one of the
 * regions. So, cells along internal boundaries in the region union will be subdivided by the
 * coverer even though this is unnecessary, using up the maxSize cell budget. Then, when the coverer
 * normalizes the covering, groups of 4 sibling cells along these internal borders will be replaced
 * by parents, resulting in coverings that may have significantly fewer than maxSize cells, and so
 * are less accurate. This is not a concern for unions of disjoint regions.
 */
public class S2RegionUnion implements S2Region, Serializable {
  // Regions is non-private so that it can be accessed from the custom field serializer.
  final S2Region[] regions;
  private transient S2Cap cachedCapBound = null;
  private transient S2LatLngRect cachedRectBound = null;

  /**
   * Creates a region that intersects or contains cells if any of the given regions does. Makes a
   * copy of 'regions'.
   */
  public S2RegionUnion(Collection<S2Region> regions) {
    this.regions = regions.toArray(new S2Region[0]);
  }

  /** Only returns true if one of the regions in the union fully contains the cell. */
  @Override
  public boolean contains(S2Cell cell) {
    for (S2Region region : regions) {
      if (region.contains(cell)) {
        return true;
      }
    }

    return false;
  }

  /** Only returns true if one of the regions contains the point. */
  @Override
  public boolean contains(S2Point point) {
    for (S2Region region : regions) {
      if (region.contains(point)) {
        return true;
      }
    }

    return false;
  }

  @Override
  public S2Cap getCapBound() {
    if (cachedCapBound != null) {
      return cachedCapBound;
    }

    cachedCapBound = S2Cap.empty();
    for (S2Region region : regions) {
      cachedCapBound = cachedCapBound.addCap(region.getCapBound());
    }
    return cachedCapBound;
  }

  @Override
  public S2LatLngRect getRectBound() {
    if (cachedRectBound != null) {
      return cachedRectBound;
    }

    cachedRectBound = S2LatLngRect.empty();
    for (S2Region region : regions) {
      cachedRectBound = cachedRectBound.union(region.getRectBound());
    }
    return cachedRectBound;
  }

  /** Returns true if the cell may intersect any region in this collection. */
  @Override
  public boolean mayIntersect(S2Cell cell) {
    for (S2Region region : regions) {
      if (region.mayIntersect(cell)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Returns true if this S2RegionUnion is equal to another S2RegionUnion, where each region must be
   * equal and in the same order. This method is intended only for testing purposes. NOTE: This
   * should be rewritten to disregard order if such functionality is ever required.
   */
  @Override
  public boolean equals(Object thatObject) {
    if (!(thatObject instanceof S2RegionUnion)) {
      return false;
    }
    S2RegionUnion that = (S2RegionUnion) thatObject;
    return Arrays.equals(regions, that.regions);
  }

  @Override
  public int hashCode() {
    return Arrays.hashCode(regions);
  }
}
