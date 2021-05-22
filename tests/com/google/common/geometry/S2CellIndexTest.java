/*
 * Copyright 2019 Google Inc.
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

import com.google.common.annotations.GwtCompatible;
import com.google.common.collect.ImmutableSet;
import com.google.common.geometry.S2CellIndex.CellIterator;
import com.google.common.geometry.S2CellIndex.ContentsIterator;
import com.google.common.geometry.S2CellIndex.Labels;
import com.google.common.geometry.S2CellIndex.NonEmptyRangeIterator;
import com.google.common.geometry.S2CellIndex.RangeIterator;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

@GwtCompatible
public class S2CellIndexTest extends GeometryTestCase {
  private final S2CellIndex index = new S2CellIndex();
  private final List<LabelledCell> contents = new ArrayList<>();

  @Override
  public void setUp() {
    super.setUp();
    index.clear();
    contents.clear();
  }

  public void testEmpty() {
    quadraticValidate();
  }

  public void testOneFaceCell() {
    add("0/", 0);
    quadraticValidate();
  }

  public void testOneLeafCell() {
    add("1/012301230123012301230123012301", 12);
    quadraticValidate();
  }

  public void testDuplicateValues() {
    add("0/", 0);
    add("0/", 0);
    add("0/", 1);
    add("0/", 17);
    quadraticValidate();
  }

  public void testDisjointCells() {
    add("0/", 0);
    add("3/", 0);
    quadraticValidate();
  }

  public void testNestedCells() {
    // Tests nested cells, including cases where several cells have the same
    // range_min() or range_max() and with randomly ordered labels.
    add("1/", 3);
    add("1/0", 15);
    add("1/000", 9);
    add("1/00000", 11);
    add("1/012", 6);
    add("1/01212", 5);
    add("1/312", 17);
    add("1/31200", 4);
    add("1/3120000", 10);
    add("1/333", 20);
    add("1/333333", 18);
    add("5/", 3);
    add("5/3", 31);
    add("5/3333", 27);
    quadraticValidate();
  }

  public void testRandomCellUnions() {
    // Construct cell unions from random S2CellIds at random levels.  Note that
    // because the cell level is chosen uniformly, there is a very high
    // likelihood that the cell unions will overlap.
    for (int i = 0; i < 100; ++i) {
      add(getRandomCellUnion(), i);
    }
    quadraticValidate();
  }

  public void testContentsIteratorSuppressesDuplicates() {
    // Checks that ContentsIterator stops reporting values once it reaches a
    // node of the cell tree that was visited by the previous call to begin().
    add("2/1", 1);
    add("2/1", 2);
    add("2/10", 3);
    add("2/100", 4);
    add("2/102", 5);
    add("2/1023", 6);
    add("2/31", 7);
    add("2/313", 8);
    add("2/3132", 9);
    add("3/1", 10);
    add("3/12", 11);
    add("3/13", 12);
    quadraticValidate();

    ContentsIterator contents = index.contents();
    expectContents("1/123", contents);
    expectContents("2/100123", contents, "2/1", 1, "2/1", 2, "2/10", 3, "2/100", 4);

    // Check that a second call with the same key yields no additional results.
    expectContents("2/100123", contents);

    // Check that seeking to a different branch yields only the new values.
    expectContents("2/10232", contents, "2/102", 5, "2/1023", 6);

    // Seek to a node with a different root.
    expectContents("2/313", contents, "2/31", 7, "2/313", 8);

    // Seek to a descendant of the previous node.
    expectContents("2/3132333", contents, "2/3132", 9);

    // Seek to an ancestor of the previous node.
    expectContents("2/213", contents);

    // A few more tests of incremental reporting.
    expectContents("3/1232", contents, "3/1", 10, "3/12", 11);
    expectContents("3/133210", contents, "3/13", 12);
    expectContents("3/133210", contents);
    expectContents("5/0", contents);

    // Now try moving backwards, which is expected to yield values that were already reported above.
    expectContents("3/13221", contents, "3/1", 10, "3/13", 12);
    expectContents("2/31112", contents, "2/31", 7);
  }

  // Tests various corner cases for the binary search optimization in visitIntersectingCells.
  public void testIntersectionOptimization() {
    add("1/001", 1);
    add("1/333", 2);
    add("2/00", 3);
    add("2/0232", 4);
    build();
    checkIntersection(makeCellUnion("1/010", "1/3"));
    checkIntersection(makeCellUnion("2/010", "2/011", "2/02"));
  }

  public void testIntersectionRandomUnions() {
    // Construct cell unions from random S2CellIds at random levels.  Note that because the cell
    // level is chosen uniformly, there is a very high likelihood that the cell unions will overlap.
    for (int i = 0; i < 100; ++i) {
      add(getRandomCellUnion(), i);
    }
    build();
    // Now repeatedly query a cell union constructed in the same way.
    for (int i = 0; i < 200; ++i) {
      checkIntersection(getRandomCellUnion());
    }
  }

  public void testIntersectionSemiRandomUnions() {
    // This test also uses random S2CellUnions, but the unions are specially constructed so that
    // interesting cases are more likely to arise.
    for (int iter = 0; iter < 200; iter++) {
      index.clear();
      S2CellId id = S2CellId.fromDebugString("1/0123012301230123");
      S2CellUnion union = new S2CellUnion();
      for (int i = 0; i < 100; i++) {
        if (oneIn(10)) {
          add(id, i);
        }
        if (oneIn(4)) {
          union.cellIds().add(id);
        }
        if (oneIn(2)) {
          id = id.nextWrap();
        }
        if (oneIn(6) && !id.isFace()) {
          id = id.parent();
        }
        if (oneIn(6) && !id.isLeaf()) {
          id = id.childBegin();
        }
      }
      union.normalize();
      build();
      checkIntersection(union);
    }
  }

  public void testLabelsEmptyNormalize() {
    Labels labels = new Labels();
    labels.normalize();
    assertEquals(0, labels.size());
  }

  /** Returns a cell union from a small number of random cells at random levels. */
  private S2CellUnion getRandomCellUnion() {
    S2CellUnion union = new S2CellUnion();
    for (int j = 0; j < 10; j++) {
      union.cellIds().add(getRandomCellId());
    }
    union.normalize();
    return union;
  }

  /** Returns a cell union of the given debug strings. */
  private static S2CellUnion makeCellUnion(String... strs) {
    S2CellUnion union = new S2CellUnion();
    for (String str : strs) {
      union.cellIds().add(S2CellId.fromDebugString(str));
    }
    return union;
  }

  /**
   * Adds the (cellId, label) pair to index and also contents (which is used for independent
   * validation).
   */
  private void add(S2CellId cellId, int label) {
    index.add(cellId, label);
    contents.add(new LabelledCell(cellId, label));
  }

  /** As {@link #add(S2CellId, int)} but accepts a debug string. */
  private void add(String cellStr, int label) {
    add(S2CellId.fromDebugString(cellStr), label);
  }

  private void add(S2CellUnion cellUnion, int label) {
    index.add(cellUnion, label);
    for (S2CellId cellId : cellUnion) {
      contents.add(new LabelledCell(cellId, label));
    }
  }

  private void build() {
    index.build();
  }

  private void quadraticValidate() {
    // Verifies that the index computes the correct set of (cellId, label) pairs for every possible
    // leaf cell.  The running time of this function is quadratic in the size of the index.
    build();
    verifyCellIterator();
    verifyIndexContents();
    verifyRangeIterators();
  }

  private void verifyCellIterator() {
    // Verifies that CellIterator visits each (cellId, label) pair exactly once.
    List<LabelledCell> actual = new ArrayList<>();
    for (CellIterator it = index.cells(); !it.done(); it.next()) {
      actual.add(new LabelledCell(it.cellId(), it.label()));
    }
    expectEqual(contents, actual);
  }

  private void verifyRangeIterators() {
    // Test finish(), which is not otherwise tested below.
    RangeIterator it = index.ranges();
    it.begin();
    it.finish();
    assertTrue(it.done());

    // And also for non-empty ranges.
    NonEmptyRangeIterator nonEmpty = index.nonEmptyRanges();
    nonEmpty.begin();
    nonEmpty.finish();
    assertTrue(nonEmpty.done());

    // Iterate through all the ranges in the index.  We simultaneously iterate through the non-empty
    // ranges and check that the correct ranges are found.
    S2CellId prevStart = S2CellId.none();
    S2CellId nonEmptyPrevStart = S2CellId.none();
    for (it.begin(), nonEmpty.begin(); !it.done(); it.next()) {
      // Check that seeking in the current range takes us to this range.
      RangeIterator it2 = index.ranges();
      S2CellId start = it.startId();
      it2.seek(it.startId());
      assertEquals(start, it2.startId());
      it2.seek(it.limitId().prev());
      assertEquals(start, it2.startId());

      // And also for non-empty ranges.
      NonEmptyRangeIterator nonEmpty2 = index.nonEmptyRanges();
      S2CellId nonEmptyStart = nonEmpty.startId();
      nonEmpty2.seek(it.startId());
      assertEquals(nonEmptyStart, nonEmpty2.startId());
      nonEmpty2.seek(it.limitId().prev());
      assertEquals(nonEmptyStart, nonEmpty2.startId());

      // Test prev() and next().
      if (it2.prev()) {
        assertEquals(prevStart, it2.startId());
        it2.next();
        assertEquals(start, it2.startId());
      } else {
        assertEquals(start, it2.startId());
        assertEquals(S2CellId.none(), prevStart);
      }

      // And also for non-empty ranges.
      if (nonEmpty2.prev()) {
        assertEquals(nonEmptyPrevStart, nonEmpty2.startId());
        nonEmpty2.next();
        assertEquals(nonEmptyStart, nonEmpty2.startId());
      } else {
        assertEquals(nonEmptyStart, nonEmpty2.startId());
        assertEquals(S2CellId.none(), nonEmptyPrevStart);
      }

      // Keep the non-empty iterator synchronized with the regular one.
      if (!it.isEmpty()) {
        assertEquals(it.startId(), nonEmpty.startId());
        assertEquals(it.limitId(), nonEmpty.limitId());
        assertFalse(nonEmpty.done());
        nonEmptyPrevStart = nonEmptyStart;
        nonEmpty.next();
      }
      prevStart = start;
    }
    // Verify that the NonEmptyRangeIterator is also finished.
    assertTrue(nonEmpty.done());
  }

  private void verifyIndexContents() {
    // Verifies that RangeIterator and ContentsIterator can be used to determine the exact set of
    // (s2cell_id, label) pairs that contain any leaf cell, where "minCellId" is the first S2CellId
    // that has not been validated yet.
    S2CellId minCellId = S2CellId.begin(S2CellId.MAX_LEVEL);
    RangeIterator range = index.ranges();
    for (range.begin(); !range.done(); range.next()) {
      assertEquals(minCellId, range.startId());
      assertEquals(-1, minCellId.compareTo(range.limitId()));
      assertTrue(range.limitId().isLeaf());
      minCellId = range.limitId();

      // Build a list of expected (cellId, label) pairs for this range.
      List<LabelledCell> expected = new ArrayList<>();
      for (LabelledCell x : contents) {
        S2CellId xMin = x.cellId.rangeMin();
        S2CellId xMax = x.cellId.rangeMax();
        if (xMin.compareTo(range.startId()) <= 0 && xMax.next().compareTo(range.limitId()) >= 0) {
          // The cell contains the entire range.
          expected.add(x);
        } else {
          // Verify that the cell does not intersect the range.
          assertFalse(
              xMin.compareTo(range.limitId().prev()) <= 0
                  && xMax.compareTo(range.startId()) >= 0);
        }
      }
      List<LabelledCell> actual = new ArrayList<>();
      ContentsIterator contents = index.contents();
      for (contents.startUnion(range); !contents.done(); contents.next()) {
        actual.add(new LabelledCell(contents.cellId(), contents.label()));
      }
      expectEqual(expected, actual);
    }
    assertEquals(S2CellId.end(S2CellId.MAX_LEVEL), minCellId);
  }

  // Tests that visitIntersectingCells() and getIntersectingLabels() return
  // correct results for the given target.
  private void checkIntersection(S2CellUnion target) {
    List<LabelledCell> expected = new ArrayList<>();
    List<LabelledCell> actual = new ArrayList<>();
    Set<Integer> expectedLabels = new HashSet<>();
    for (CellIterator it = index.cells(); !it.done(); it.next()) {
      if (target.intersects(it.cellId())) {
        expected.add(new LabelledCell(it.cellId(), it.label()));
        expectedLabels.add(it.label());
      }
    }
    index.visitIntersectingCells(target, (cell_id, label) -> {
      actual.add(new LabelledCell(cell_id, label));
      return true;
    });
    expectEqual(expected, actual);
    assertEquals(expectedLabels, ImmutableSet.copyOf(index.getIntersectingLabels(target)));
  }

  private void expectContents(String targetStr, ContentsIterator contents, Object... strLabel) {
    // Given an S2CellId "targetStr" in human-readable form, expects that the first leaf cell
    // contained by this target will intersect the exact set of (cellId, label) pairs expected by
    // "strLabel".
    RangeIterator range = index.ranges();
    range.seek(S2CellId.fromDebugString(targetStr).rangeMin());
    List<LabelledCell> expected = new ArrayList<>();
    List<LabelledCell> actual = new ArrayList<>();
    for (int i = 0; i < strLabel.length; i += 2) {
      expected.add(
          new LabelledCell(
              S2CellId.fromDebugString((String) strLabel[i]), (Integer) strLabel[i + 1]));
    }
    for (contents.startUnion(range); !contents.done(); contents.next()) {
      actual.add(new LabelledCell(contents.cellId(), contents.label()));
    }
    expectEqual(expected, actual);
  }

  private static void expectEqual(List<LabelledCell> expected, List<LabelledCell> actual) {
    // Verifies that "expected" and "actual" have the same contents.  Note duplicate are allowed.
    Comparator<LabelledCell> order = Comparator.comparing(x -> x.cellId);
    order = order.thenComparingInt(x -> x.label);
    Collections.sort(expected, order);
    Collections.sort(actual, order);
    assertEquals(expected, actual);
  }

  private static final class LabelledCell {
    private final S2CellId cellId;
    private final int label;

    LabelledCell(S2CellId cell, int label) {
      this.cellId = cell;
      this.label = label;
    }

    @Override
    public boolean equals(Object o) {
      if (o instanceof LabelledCell) {
        LabelledCell l = (LabelledCell) o;
        return cellId.equals(l.cellId) && label == l.label;
      } else {
        return false;
      }
    }

    @Override
    public int hashCode() {
      return cellId.hashCode();
    }

    @Override
    public String toString() {
      return cellId.toToken() + "=" + label;
    }
  }
}
