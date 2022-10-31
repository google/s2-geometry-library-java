/*
 * Copyright 2015 Google Inc.
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

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import com.google.common.geometry.S2PointIndex.Entry;

/** Verifies {@link S2PointIndex}. */
public final class S2PointIndexTest extends GeometryTestCase {
  private final S2PointIndex<Integer> index = new S2PointIndex<>();
  private Multiset<S2PointIndex.Entry<Integer>> contents;

  @Override
  public void setUp() {
    super.setUp();
    contents = HashMultiset.create();
    index.reset();
  }

  /** Adds an entry to the index and the contents. */
  private void add(S2Point point, int data) {
    Entry<Integer> entry = S2PointIndex.createEntry(point, data);
    index.add(entry);
    contents.add(entry);
  }

  /** Verifies the index contents against {@link #contents}. */
  private void verify() {
    Multiset<S2PointIndex.Entry<Integer>> entries = HashMultiset.create();
    for (S2Iterator<Entry<Integer>> it = index.iterator(); !it.done(); it.next()) {
      entries.add(it.entry());
    }
    assertEquals(contents, entries);
  }

  /** Verifies the iterator methods. This is somewhat redundant with S2ShapeIndexTest. */
  public void checkIteratorMethods() {
    S2Iterator<Entry<Integer>> it = index.iterator();
    assertTrue(it.atBegin());
    it.finish();
    assertTrue(it.done());

    // Iterate through all the cells in the index.
    S2CellId prev = S2CellId.none();
    S2CellId min = S2CellId.begin(S2CellId.MAX_LEVEL);
    for (it.restart(); !it.done(); it.next()) {
      S2CellId id = it.id();
      assertEquals(id, S2CellId.fromPoint(it.entry().point()));

      S2Iterator<Entry<Integer>> it2 = index.iterator();
      if (id.equals(prev)) {
        it2.seek(id);
      }

      // Generate a cellunion that covers the range of empty leaf cells between the last cell and
      // this one. Then make sure that seeking to any of those cells takes us to the immediately
      // following cell.
      S2CellUnion skipped = new S2CellUnion();
      skipped.initFromBeginEnd(min, id.rangeMin());
      for (int i = 0; i < skipped.size(); ++i) {
        it2.seek(skipped.cellId(i));
        assertEquals(id, it2.id());
      }
      // Test prev(), next(), seek(), and seekForward().
      if (prev.isValid()) {
        assertFalse(it.atBegin());
        it2 = index.iterator();
        it2.position(it);
        it2.prev();
        assertEquals(prev, it2.id());
        it2.next();
        assertEquals(id, it2.id());
        it2.seek(prev);
        assertEquals(prev, it2.id());
        it2.seekForward(id);
        assertEquals(id, it2.id());
        it2.seekForward(prev);
        assertEquals(id, it2.id());
      }
      prev = id;
      min = id.rangeMax().next();
    }
  }

  public void testNoPoints() {
    checkIteratorMethods();
  }

  public void testRandomPoints() {
    for (int i = 0; i < 1000; ++i) {
      add(data.getRandomPoint(), (int) data.uniform(0, 100));
    }
    verify();
    checkIteratorMethods();
  }

  public void testEntryEquality() {
    S2Point point = data.getRandomPoint();
    int data = 10;
    Entry<Integer> entry = S2PointIndex.createEntry(point, data);
    Entry<Integer> entry2 = S2PointIndex.createEntry(point, data);
    assertEquals(entry, entry2);
    entry2 = S2PointIndex.createEntry(point, 20);
    assertFalse(entry.equals(entry2));
    entry2 = S2PointIndex.createEntry(point, null);
    assertFalse(entry.equals(entry2));
    assertFalse(entry2.equals(entry));
    entry = S2PointIndex.createEntry(point, null);
    assertEquals(entry, entry2);
  }
}
