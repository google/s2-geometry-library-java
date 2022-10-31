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

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Iterables;
import java.util.List;
import junit.framework.TestCase;

public class S2RegionSharderTest extends TestCase {
  private static final int DEFAULT_SHARD = 42;

  public void testGetMostIntersectingShard() {
    ImmutableList<S2CellUnion> coverings = ImmutableList.of(
        S2CellUnion.copyFrom(ImmutableList.of(
            S2CellId.fromFacePosLevel(0, 0, 10))),
        S2CellUnion.copyFrom(ImmutableList.of(
            S2CellId.fromFacePosLevel(1, 1, 9),
            S2CellId.fromFacePosLevel(3, 0, 8))),
        S2CellUnion.copyFrom(ImmutableList.of(
            S2CellId.fromFacePosLevel(5, 0, 10))));

    S2RegionSharder sharder = new S2RegionSharder(coverings);

    // Overlap with only 1 shard.
    assertEquals(0, sharder.getMostIntersectingShard(
        S2CellUnion.copyFrom(ImmutableList.of(
            S2CellId.fromFacePosLevel(0, 0, 11))), DEFAULT_SHARD));

    // Overlap with multiple shards, picks the shard with more overlap.
    assertEquals(1, sharder.getMostIntersectingShard(
        S2CellUnion.copyFrom(ImmutableList.of(
            S2CellId.fromFacePosLevel(0, 0, 10),
            S2CellId.fromFacePosLevel(3, 0, 9),
            S2CellId.fromFacePosLevel(3, 1, 9))),
        DEFAULT_SHARD));

    // Overlap with no shards.
    assertEquals(DEFAULT_SHARD, sharder.getMostIntersectingShard(
        S2CellUnion.copyFrom(ImmutableList.of(S2CellId.fromFacePosLevel(4, 0, 10))),
        DEFAULT_SHARD));
  }

  public void testGetIntersectingShards() {
    List<S2CellUnion> coverings = ImmutableList.of(
        S2CellUnion.copyFrom(ImmutableList.of(S2CellId.fromFacePosLevel(0, 0, 10))),
        S2CellUnion.copyFrom(ImmutableList.of(
            S2CellId.fromFacePosLevel(1, 1, 9),
            S2CellId.fromFacePosLevel(3, 0, 8))),
        S2CellUnion.copyFrom(ImmutableList.of(S2CellId.fromFacePosLevel(5, 0, 10))));

    S2RegionSharder sharder = new S2RegionSharder(coverings);

    // Overlap with only 1 shard
    assertEquals(
        ImmutableSet.of(0),
        ImmutableSet.copyOf(sharder.getIntersectingShards(
            S2CellUnion.copyFrom(ImmutableList.of(S2CellId.fromFacePosLevel(0, 0, 11))))));

    // Overlap with multiple shards, picks the shard with more overlap.
    assertEquals(
        ImmutableSet.of(0, 1),
        ImmutableSet.copyOf(sharder.getIntersectingShards(
            S2CellUnion.copyFrom(ImmutableList.of(
                S2CellId.fromFacePosLevel(0, 0, 10),
                S2CellId.fromFacePosLevel(3, 0, 9),
                S2CellId.fromFacePosLevel(3, 1, 9))))));

    // Overlap with no shards.
    assertTrue(Iterables.isEmpty(sharder.getIntersectingShards(
        S2CellUnion.copyFrom(ImmutableList.of(S2CellId.fromFacePosLevel(4, 0, 10))))));
  }
}
