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

import com.google.common.collect.Iterables;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A sharding function that provides shard IDs whose boundaries intersect an {@link S2Region}. This
 * class is especially suited to testing regions against shards that are usually very different in
 * size from the regions; in that case, {@link S2CellUnion coverings} of the region tend to cover
 * too much area (any simple covering of the Pacific Ocean for example) or too complex (small
 * regions are often contained by a single cell of the shard's covering).
 */
public class S2RegionSharder {
  private final S2CellIndex index = new S2CellIndex();

  /**
   * Creates a new sharder.
   * @param boundaries the boundaries of each shard indexed by shard ID
   */
  public S2RegionSharder(List<S2CellUnion> boundaries) {
    for (int i = 0; i < boundaries.size(); i++) {
      index.add(boundaries.get(i), i);
    }
    index.build();
  }

  /**
   * Returns an index into the original list of {@code S2CellUnion} given to the constructor, which
   * indicates the shard whose covering has the most overlap with {@code region}, or returns
   * {@code defaultShard} if no shard overlaps the region.
   */
  public int getMostIntersectingShard(S2Region region, int defaultShard) {
    // Return the best shard by intersection area.
    Map<Integer, S2CellUnion> shardCoverings = intersections(region);
    int bestShard = defaultShard;
    long bestSum = 0;
    for (int shardId : shardCoverings.keySet()) {
      S2CellUnion shardCovering = shardCoverings.get(shardId);
      long sum = 0;
      for (S2CellId id : shardCovering.cellIds()) {
        sum += id.lowestOnBit();
      }
      if (sum > bestSum) {
        bestShard = shardId;
        bestSum = sum;
      }
    }
    return bestShard;
  }

  /**
   * Returns a list of shard numbers which intersect with {@code region}. Shard numbers are not
   * guaranteed to be sorted in any particular order. If no shards overlap, returns an empty list.
   */
  public Iterable<Integer> getIntersectingShards(S2Region region) {
    return intersections(region).keySet();
  }

  private Map<Integer, S2CellUnion> intersections(S2Region region) {
    // Compute the intersection between the region covering and each shard covering.
    S2CellUnion regionCovering = new S2CellUnion();
    // TODO(eengle): use GetCellUnionBound when available.
    S2RegionCoverer.DEFAULT.getFastCovering(region.getCapBound(), regionCovering.cellIds());
    Map<Integer, S2CellUnion> shardCoverings = new HashMap<>();
    index.visitIntersectingCells(regionCovering, (cell, shardId) -> {
      S2CellUnion shard = shardCoverings.get(shardId);
      if (shard == null) {
        shard = new S2CellUnion();
        shardCoverings.put(shardId, shard);
      }
      shard.cellIds().add(cell);
      return true;
    });

    // The fast covering is very loose, but typically it only intersects one shard.
    if (shardCoverings.size() == 1) {
      return shardCoverings;
    }

    // Clip each shard to the region.
    S2CellUnion tempIntersection = new S2CellUnion();
    for (S2CellUnion shardCovering : shardCoverings.values()) {
      // Get the intersection since the shard's covering may have smaller cells than the region's.
      // We know the region covering is normalized but must normalize each shard covering first.
      shardCovering.normalize();
      tempIntersection.getIntersection(shardCovering, regionCovering);
      // Remove cells that don't intersect the region. Since the fast covering tends to intersect
      // multiple shards for regions that are near a boundary between shards, there is a high chance
      // of testing the region against the same cell more than once. These tests can be expensive,
      // but rather than make this algorithm more complicated, we choose to push clients toward
      // regions that are fast, such as S2ShapeIndexRegion.
      Iterables.removeIf(tempIntersection, id -> !region.mayIntersect(new S2Cell(id)));
      // Save the clipped shard covering.
      shardCovering.cellIds().clear();
      shardCovering.cellIds().addAll(tempIntersection.cellIds());
    }
    Iterables.removeIf(shardCoverings.values(), cells -> cells.size() == 0);
    return shardCoverings;
  }
}
