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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import java.util.List;
import java.util.Set;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link S2RegionSharder}. */
@RunWith(JUnit4.class)
public class S2RegionSharderTest {
  private static final int DEFAULT_SHARD = 42;

  @Test
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

  @Test
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
    assertTrue(
        sharder
            .getIntersectingShards(
                S2CellUnion.copyFrom(ImmutableList.of(S2CellId.fromFacePosLevel(4, 0, 10))))
            .isEmpty());
  }

  @Test
  public void testGetIntersectingShards_sortedCellIds() {
    // Test a situation where S2Region.getCellUnionBound produces an unsorted list of cell IDs,
    // triggering a binary search failure in S2CellIndex.visitIntersectingCells. This example was
    // observed to fail with an unnormalized covering.

    // Two points near the corner of Faces 0 (+X), 1 (+Y), and 2 (+Z). p0 is on F0, p1 is on F1.
    S2Point p0 = new S2Point(1.0, 0.99, 0.99).normalize();
    S2Point p1 = new S2Point(0.99, 1.0, 0.99).normalize();

    S2ShapeIndex index = new S2ShapeIndex();
    index.add(S2Point.Shape.singleton(p0));
    index.add(S2Point.Shape.singleton(p1));

    S2ShapeIndexBufferedRegion bufferedRegion =
        new S2ShapeIndexBufferedRegion(index, S1Angle.degrees(0.1));

    // Create a single shard on the corner of Face 1 and pass it to a new S2RegionSharder.
    S2CellId shardCell = S2CellId.fromPoint(new S2Point(0.999, 1.0, 0.999).normalize()).parent(10);
    ImmutableList<S2CellUnion> coverings =
        ImmutableList.of(S2CellUnion.copyFrom(ImmutableList.of(shardCell)));
    S2RegionSharder sharder = new S2RegionSharder(coverings);

    Set<Integer> intersectingShards = sharder.getIntersectingShards(bufferedRegion);
    assertTrue("Shard should intersect the buffered region", intersectingShards.contains(0));
  }

  @Test
  public void testGetIntersectingShardsTieBreaking() {
    S2CellId c0 = S2CellId.fromFace(0).child(0);
    S2CellId c1 = S2CellId.fromFace(1).child(0);
    assertEquals(c0.lowestOnBit(), c1.lowestOnBit());

    ImmutableList<S2CellUnion> coverings =
        ImmutableList.of(
            S2CellUnion.copyFrom(ImmutableList.of(c0)), S2CellUnion.copyFrom(ImmutableList.of(c1)));

    for (int i = 0; i < 2; ++i) {
      // Tie-breaking: if two shards have same intersection area, pick lowest index.
      //
      // Shards 0 and 1 cover c1 and c0 on the first iteration and c0 and c1 on the second
      // iteration.
      // With region c0+c1, shard 0 intersects c0, shard 1 intersects c1.
      // Intersection sums are equal, so the shard with smaller index must be chosen.
      S2RegionSharder sharder = new S2RegionSharder(coverings);
      assertEquals(
          0, sharder.getMostIntersectingShard(S2CellUnion.copyFrom(ImmutableList.of(c0, c1)), 42));

      // Swap the order of the shards and repeat.
      coverings =
          ImmutableList.of(
              S2CellUnion.copyFrom(ImmutableList.of(c1)),
              S2CellUnion.copyFrom(ImmutableList.of(c0)));
    }
  }
}
