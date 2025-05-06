/*
 * Copyright 2024 Google Inc.
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
import static java.util.Comparator.naturalOrder;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.errorprone.annotations.Keep;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

@RunWith(JUnit4.class)
public final class S2CellIteratorJoinTest extends GeometryTestCase {

  // Cell IDs covering Central Park in New York City
  private static final String[] kCentralParkATokens =
      new String[] {
        "89c2589", "89c258a1", "89c258a3", "89c258bc", "89c258c1", "89c258ec", "89c258f4"
      };

  // Cell Ids also covering Central Park but a subset of CentralParkA.
  private static final String[] kCentralParkBTokens =
      new String[] {"89c2589", "89c258a03", "89c258a1c", "89c258a3", "89c258bd", "89c258be1"};

  private static S2Iterator<IntEntry> makeListS2CellIterator(String[] tokens) {
    return S2Iterator.fromList(tokenList(tokens));
  }

  // Test that we can join two different types of iterators.
  @Test
  public void hetrogenousJoinCompiles() {
    new S2CellIteratorJoin<>(new S2ShapeIndex().iterator(), new S2PointIndex<Integer>().iterator())
        .join((iterA, iterB) -> true);
  }

  @Test
  public void exactJoinWorks() {
    S2Iterator<IntEntry> listIterA = makeListS2CellIterator(kCentralParkATokens);
    S2Iterator<IntEntry> listIterB = makeListS2CellIterator(kCentralParkBTokens);

    List<Pair<S2CellId, S2CellId>> rows = new ArrayList<>();
    new S2CellIteratorJoin<>(listIterA, listIterB)
        .join(
            (iterA, iterB) -> {
              S2CellId a = iterA.id();
              S2CellId b = iterB.id();
              assertTrue(a.contains(b));
              rows.add(Pair.of(a, b));
              return true;
            });

    assertEquals(min(kCentralParkATokens.length, kCentralParkBTokens.length), rows.size());

    ImmutableList<Pair<S2CellId, S2CellId>> truth =
        ImmutableList.of(
            pairFromTokens("89c2589", "89c2589"),
            pairFromTokens("89c258a1", "89c258a03"),
            pairFromTokens("89c258a1", "89c258a1c"),
            pairFromTokens("89c258a3", "89c258a3"),
            pairFromTokens("89c258bc", "89c258bd"),
            pairFromTokens("89c258bc", "89c258be1"));

    assertEquals(rows.size(), truth.size());
    for (int i = 0; i < truth.size(); ++i) {
      assertEquals(truth.get(i), rows.get(i));
    }
  }

  @Test
  public void exactFalseJoinReturnsImmediately() {
    S2Iterator<IntEntry> listIterA = makeListS2CellIterator(kCentralParkATokens);
    S2Iterator<IntEntry> listIterB = makeListS2CellIterator(kCentralParkBTokens);

    List<Pair<S2CellId, S2CellId>> rows = new ArrayList<>();
    boolean cancelled =
        new S2CellIteratorJoin<>(listIterA, listIterB)
            .join(
                (iterA, iterB) -> {
                  rows.add(Pair.of(iterA.id(), iterB.id()));
                  return false;
                });
    assertFalse(cancelled);
    assertEquals(1, rows.size());
  }

  @Test
  public void tolerantFalseJoinReturnsImmediately() {
    S2Iterator<IntEntry> listIterA = makeListS2CellIterator(kCentralParkATokens);
    S2Iterator<IntEntry> listIterB = makeListS2CellIterator(kCentralParkBTokens);

    // Just need a non-zero tolerance to trigger tolerant join logic.
    S1ChordAngle dist = S1ChordAngle.fromDegrees(.001);

    List<Pair<S2CellId, S2CellId>> rows = new ArrayList<>();
    boolean cancelled =
        new S2CellIteratorJoin<>(listIterA, listIterB, dist)
            .join(
                (iterA, iterB) -> {
                  rows.add(Pair.of(iterA.id(), iterB.id()));
                  return false;
                });

    assertFalse(cancelled);
    assertEquals(1, rows.size());
  }

  @Test
  public void exactJoinSeekingWorks() {
    // Sometimes we have to seek more than once to find a pair of cells that overlap. 2d5e3 below
    // doesn't overlap anything in listIterB, and so shouldn't be reported. Instead the join should
    // seek twice and skip it, make sure it does.
    S2Iterator<IntEntry> listIterA =
        makeListS2CellIterator(new String[] {"2d5dd7", "2d5ddc", "2d5e3", "2d5e801", "2d5e803"});
    S2Iterator<IntEntry> listIterB = makeListS2CellIterator(new String[] {"2d5d", "2d5e84"});

    ImmutableList<Pair<S2CellId, S2CellId>> truth =
        ImmutableList.of(
            pairFromTokens("2d5dd7", "2d5d"), //
            pairFromTokens("2d5ddc", "2d5d"), //
            pairFromTokens("2d5e801", "2d5e84"), //
            pairFromTokens("2d5e803", "2d5e84") //
            );

    List<Pair<S2CellId, S2CellId>> rows = new ArrayList<>();
    new S2CellIteratorJoin<>(listIterA, listIterB)
        .join(
            (iterA, iterB) -> {
              rows.add(Pair.of(iterA.id(), iterB.id()));
              return true;
            });

    assertEquals(truth.size(), rows.size());
    for (int i = 0; i < truth.size(); ++i) {
      assertEquals(truth.get(i), rows.get(i));
    }
  }

  @Test
  public void nearJoinWorks() {
    S2Iterator<IntEntry> listIterA = makeListS2CellIterator(kCentralParkATokens);
    S2Iterator<IntEntry> listIterB = makeListS2CellIterator(kCentralParkBTokens);

    S1ChordAngle kTolerance = S1ChordAngle.fromDegrees(1);

    Set<Pair<S2CellId, S2CellId>> rows = new HashSet<>();
    new S2CellIteratorJoin<>(listIterA, listIterB, kTolerance)
        .join(
            (iterA, iterB) -> {
              rows.add(Pair.of(iterA.id(), iterB.id()));
              return true;
            });

    ImmutableList<Pair<S2CellId, S2CellId>> truth =
        ImmutableList.of(
            pairFromTokens("89c2589", "89c2589"),
            pairFromTokens("89c258a1", "89c258a03"),
            pairFromTokens("89c258a1", "89c258a1c"),
            pairFromTokens("89c258a3", "89c258a3"),
            pairFromTokens("89c258bc", "89c258bd"),
            pairFromTokens("89c258bc", "89c258be1"));

    // The contents of the exact join should be a subset of the rough join.
    for (Pair<S2CellId, S2CellId> row : truth) {
      assertTrue(rows.contains(row));
    }

    // All the pairs should be within tolerance of each other.
    for (Pair<S2CellId, S2CellId> row : truth) {
      S1ChordAngle distance = new S2Cell(row.first).getDistance(new S2Cell(row.second));
      assertTrue(distance.lessOrEquals(kTolerance));
    }

    // These cells are more than 0 degrees apart (i.e. not touching) and not part of the exact
    // results, but should be included in the tolerant result.
    ImmutableList<Pair<S2CellId, S2CellId>> tolerantTruth =
        ImmutableList.of(
            pairFromTokens("89c258a1", "89c258bd"),
            pairFromTokens("89c258a1", "89c258be1"),
            pairFromTokens("89c258a3", "89c258a03"),
            pairFromTokens("89c258a3", "89c258be1"),
            pairFromTokens("89c258bc", "89c258a03"),
            pairFromTokens("89c258bc", "89c258a1c"),
            pairFromTokens("89c258c1", "89c258a03"),
            pairFromTokens("89c258c1", "89c258a1c"),
            pairFromTokens("89c258c1", "89c258a3"),
            pairFromTokens("89c258c1", "89c258bd"),
            pairFromTokens("89c258c1", "89c258be1"),
            pairFromTokens("89c258ec", "89c258a03"),
            pairFromTokens("89c258ec", "89c258a1c"),
            pairFromTokens("89c258ec", "89c258a3"),
            pairFromTokens("89c258ec", "89c258bd"),
            pairFromTokens("89c258ec", "89c258be1"),
            pairFromTokens("89c258f4", "89c258a03"),
            pairFromTokens("89c258f4", "89c258a1c"),
            pairFromTokens("89c258f4", "89c258a3"),
            pairFromTokens("89c258f4", "89c258bd"),
            pairFromTokens("89c258f4", "89c258be1"));

    for (Pair<S2CellId, S2CellId> row : tolerantTruth) {
      assertTrue(rows.contains(row));
    }
  }

  // Verify that the join returns all the A cells in a contiguous manner.
  @Test
  public void tolerantJoinIsLeftDriven() {
    // Coordinate frame to use to generate fractal. We intentionally build a fractal that spans the
    // boundary between faces.
    Matrix kFrame = upFrameAt(0, -45);

    S2FractalBuilder fractal = new S2FractalBuilder(data.rand);
    fractal.setLevelForApproxMaxEdges(100);

    S2Polygon polygon = new S2Polygon(fractal.makeLoop(kFrame, S1Angle.degrees(10)));
    S2ShapeIndex index = polygon.index();
    index.applyUpdates();

    Set<S2CellId> cellsSeen = new HashSet<>();
    S2CellId[] currCell = new S2CellId[] {S2CellId.sentinel()};

    new S2CellIteratorJoin<>(index.iterator(), index.iterator(), S1ChordAngle.fromDegrees(2))
        .join(
            (iterA, iterB) -> {
              if (!iterA.id().equals(currCell[0])) {
                boolean aIsNew = !cellsSeen.contains(iterA.id());
                assertTrue(aIsNew);
                currCell[0] = iterA.id();
                cellsSeen.add(currCell[0]);
                return aIsNew;
              }
              return true;
            });
  }

  // Verify that the join does indeed return _all_ pairs of cells that are within a tolerance of
  // each other.
  @Test
  public void allPairsSeen() {
    // Coordinate frame to use to generate fractal. We intentionally build a
    // fractal that spans the boundary between faces.
    Matrix kFrame = upFrameAt(0, -45);

    S2FractalBuilder fractal = new S2FractalBuilder(data.rand);
    fractal.setLevelForApproxMaxEdges(1000);

    S2Polygon polygon = new S2Polygon(fractal.makeLoop(kFrame, S1Angle.degrees(10)));
    S2ShapeIndex index = polygon.index();
    index.applyUpdates();

    // Pull cells out of the index.
    List<S2Cell> cells = new ArrayList<>();
    S2Iterator<S2ShapeIndex.Cell> iter = index.iterator();
    while (!iter.done()) {
      cells.add(new S2Cell(iter.id()));
      iter.next();
    }

    // Build all pairs that are closer than the tolerance by brute force.
    S1ChordAngle kTolerance = S1ChordAngle.fromDegrees(2);
    Set<Pair<S2CellId, S2CellId>> brutePairs = new HashSet<>();
    for (S2Cell cell0 : cells) {
      for (S2Cell cell1 : cells) {
        if (cell0.getDistance(cell1).lessThan(kTolerance)) {
          brutePairs.add(Pair.of(cell0.id(), cell1.id()));
        }
      }
    }

    Set<Pair<S2CellId, S2CellId>> joinPairs = new HashSet<>();
    new S2CellIteratorJoin<>(index.iterator(), index.iterator(), kTolerance)
        .join(
            (iterA, iterB) -> {
              joinPairs.add(Pair.of(iterA.id(), iterB.id()));
              return true;
            });

    assertEquals(brutePairs, joinPairs);
  }

  // Verifies that simple join produces expected results.
  @Test
  public void simple() {
    S2CellId a = S2CellId.fromFace(0);
    S2CellId b = S2CellId.fromFace(1);
    S2CellUnion left = new S2CellUnion();
    left.cellIds().add(a);
    left.cellIds().add(b);
    S2CellUnion right = new S2CellUnion();
    right.cellIds().add(b.child(0));
    right.cellIds().add(b.child(2));
    List<S2CellId> results = new ArrayList<>();
    assertTrue(new S2CellIteratorJoin<>(left.s2Iterator(), right.s2Iterator())
        .simpleJoin(c -> results.add(c.entry())));
    assertEquals(ImmutableList.of(b), results);
  }

  @SuppressWarnings("FloatingPointLiteralPrecision") // to have the same constants as the C++ test.
  @Test
  public void b299938257Regression() {
    // This used to trigger a bug where the join wasn't checking for the end of the iterator before
    // de-referencing it to check for the cell id.
    S2PointIndex<Integer> pointIndex = new S2PointIndex<>();
    S2Point[] points =
        new S2Point[] {
          new S2Point(0.998782953991165789, -0.034851647907011431, -0.034899476426537568),
          new S2Point(1.000000000000000000, -0.000000000000005489, -0.000000000000005494),
          new S2Point(0.998782953991165789, -0.034851647907011431, 0.034899476426537568),
          new S2Point(1.000000000000000000, -0.000000000000005489, 0.000000000000005494)
        };

    for (S2Point point : points) {
      pointIndex.add(point, 0);
    }

    Matrix kFrame = upFrameAt(0, 0);

    S2FractalBuilder fractal = new S2FractalBuilder(data.rand);
    fractal.setLevelForApproxMaxEdges(100);

    S2Polygon polygon = new S2Polygon(fractal.makeLoop(kFrame, S1Angle.degrees(1)));
    S2ShapeIndex index = polygon.index();
    index.applyUpdates();

    int[] count = {0};
    new S2CellIteratorJoin<>(index.iterator(), pointIndex.iterator(), S1ChordAngle.fromDegrees(0.5))
        .join(
            (iterA, iterB) -> {
              ++count[0];
              return true;
            });

    // This is arbitrary, but for the given fractal and radius we should see eight cells within
    // range of the point index.
    assertEquals(8, count[0]);
  }

  private Pair<S2CellId, S2CellId> pairFromTokens(String tokenA, String tokenB) {
    return Pair.of(S2CellId.fromToken(tokenA), S2CellId.fromToken(tokenB));
  }

  // Returns a frame in the "up" (positive z direction at a given point). Use for fractals on the
  // equator.
  private Matrix upFrameAt(double lat, double lng) {
    return TestDataGenerator.upFrameAt(S2LatLng.fromDegrees(lat, lng));
  }

  /**
   * Given an array of strings which are valid S2CellId tokens, returns a List of IntEntries with
   * consecutive integers.
   */
  private static List<IntEntry> tokenList(String[] tokens) {
    List<S2CellId> cellIds = new ArrayList<>();
    for (String token : tokens) {
      cellIds.add(S2CellId.fromToken(token));
    }
    cellIds.sort(naturalOrder());

    ArrayList<IntEntry> entries = new ArrayList<>();
    int count = 0;
    for (S2CellId cellId : cellIds) {
      entries.add(new IntEntry(cellId, ++count));
    }
    return entries;
  }

  private static class IntEntry implements S2Iterator.Entry {
    private final S2CellId id;
    private final Integer value;

    public IntEntry(S2CellId id, Integer value) {
      this.id = id;
      this.value = value;
    }

    @Override
    public long id() {
      return id.id();
    }

    @Keep
    public Integer value() {
      return value;
    }

    @Override
    public boolean equals(Object o) {
      if (!(o instanceof IntEntry)) {
        return false;
      }
      IntEntry other = (IntEntry) o;
      return id.equals(other.id) && value.equals(other.value);
    }

    @Override
    public int hashCode() {
      return id.hashCode() + value.hashCode();
    }
  }
}
