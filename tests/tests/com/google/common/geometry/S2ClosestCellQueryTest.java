// Copyright 2026 Google LLC
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
package com.google.common.geometry;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.geometry.S2ClosestCellQuery.Result;
import java.util.List;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public class S2ClosestCellQueryTest extends GeometryTestCase {

  private S2CellIndex index;
  private S2ClosestCellQuery query;

  @Before
  public void setupIndex() {
    index = new S2CellIndex();
  }

  @Test
  public void testEmptyIndex() {
    index.build();
    query = new S2ClosestCellQuery(index);
    Result result = query.findClosestCell(new S2ClosestCellQuery.PointTarget(new S2Point(1, 0, 0)));
    assertTrue(result.isEmpty());
    assertEquals(S1ChordAngle.INFINITY, result.distance());
  }

  @Test
  public void testFindClosestCell() {
    S2CellId cellId1 = S2CellId.fromToken("1");
    S2CellId cellId2 = S2CellId.fromToken("3");

    index.add(cellId1, 1);
    index.add(cellId2, 2);
    index.build();

    query = new S2ClosestCellQuery(index);

    // Target is cellId1
    S2ClosestCellQuery.Target target = new S2ClosestCellQuery.PointTarget(cellId1.toPoint());
    Result result = query.findClosestCell(target);
    assertFalse(result.isEmpty());
    assertEquals(1, result.label());
    assertTrue(result.distance().isZero());
  }

  @Test
  public void testMaxDistance() {
    S2CellId cellId1 = S2CellId.fromFace(0);
    S2CellId cellId2 = S2CellId.fromFace(2); // Far away

    index.add(cellId1, 1);
    index.add(cellId2, 2);
    index.build();

    query = new S2ClosestCellQuery(index);
    query.setMaxDistance(S1Angle.degrees(10));

    S2ClosestCellQuery.Target target = new S2ClosestCellQuery.PointTarget(cellId1.toPoint());
    List<Result> results = query.findClosestCells(target);
    assertEquals(1, results.size());
    assertEquals(1, results.get(0).label());
  }

  @Test
  public void testMaxResults() {
    for (int i = 0; i < 10; i++) {
      index.add(S2CellId.fromFace(0).childBegin(10).advance(i), i);
    }
    index.build();

    query = new S2ClosestCellQuery(index);
    query.setMaxResults(5);

    S2ClosestCellQuery.Target target = new S2ClosestCellQuery.PointTarget(S2CellId.fromFace(0).toPoint());
    List<Result> results = query.findClosestCells(target);
    assertEquals(5, results.size());
  }

  @Test
  public void testIsDistanceLess() {
    S2CellId cellId = S2CellId.fromFace(0);
    index.add(cellId, 1);
    index.build();

    query = new S2ClosestCellQuery(index);
    S2ClosestCellQuery.Target target = new S2ClosestCellQuery.PointTarget(cellId.toPoint());

    assertTrue(query.isDistanceLess(target, S1ChordAngle.fromDegrees(10)));
    assertFalse(query.isDistanceLess(target, S1ChordAngle.ZERO));
    assertTrue(query.isDistanceLessOrEqual(target, S1ChordAngle.ZERO));
  }
}
