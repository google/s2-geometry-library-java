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

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertNotEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;
import static org.junit.Assert.assertSame;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import com.google.common.geometry.S2PointIndex.Entry;
import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;
import java.util.Comparator;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Verifies {@link S2PointIndex}. */
@RunWith(JUnit4.class)
public final class S2PointIndexTest extends GeometryTestCase {
  private final S2PointIndex<Integer> index = new S2PointIndex<>();
  private Multiset<S2PointIndex.Entry<Integer>> contents;

  @Before
  public void reset() {
    contents = HashMultiset.create();
    index.reset();
  }

  /** Adds an entry to the index and the contents. */
  private void add(S2Point point, int data) {
    Entry<Integer> entry = S2PointIndex.createEntry(point, data);
    System.out.println("Adding entry: " + entry.point().toDegreesString() + ", " + entry.data());
    index.add(entry);
    contents.add(entry);
  }

  /** Verifies the index contents against {@link #contents}. */
  private void verify() {
    System.out.println("Verifying index contents");
    Multiset<S2PointIndex.Entry<Integer>> entries = HashMultiset.create();
    for (S2Iterator<Entry<Integer>> it = index.iterator(); !it.done(); it.next()) {
      Entry<Integer> entry = it.entry();
      System.out.println("  Got entry: " + entry.point().toDegreesString() + ", " + entry.data());
      entries.add(entry);
    }
    assertEquals(contents, entries);
  }

  @Test
  public void testEntryOrdering() {
    S2Point point = S2Point.X_POS;
    Entry<Integer> xposEntry1 = S2PointIndex.createEntry(point, 1);
    Entry<Integer> xposEntry2 = S2PointIndex.createEntry(point, 2);
    Entry<Integer> xposEntry2b = S2PointIndex.createEntry(point, 2);
    Entry<Integer> xposEntryNull = S2PointIndex.createEntry(point, null);
    Entry<Integer> xposEntryNullb = S2PointIndex.createEntry(point, null);

    // Ordering by different data values is stable with stableOrder().
    assertEquals(-1, Entry.<Integer>stableOrder().compare(xposEntry1, xposEntry2));
    assertEquals(1, Entry.<Integer>stableOrder().compare(xposEntry2, xposEntry1));
    // Ordering by different data values is not stable with order().
    assertEquals(-1, Entry.<Integer>order().compare(xposEntry1, xposEntry2));
    assertEquals(-1, Entry.<Integer>order().compare(xposEntry2, xposEntry1));

    // Null is always first for both order() and stableOrder().
    assertEquals(-1, Entry.<Integer>stableOrder().compare(xposEntryNull, xposEntry1));
    assertEquals(1, Entry.<Integer>stableOrder().compare(xposEntry1, xposEntryNull));
    assertEquals(-1, Entry.<Integer>order().compare(xposEntryNull, xposEntry1));
    assertEquals(1, Entry.<Integer>order().compare(xposEntry1, xposEntryNull));

    // Identical instances for both order() and stableOrder().
    assertEquals(0, Entry.<Integer>stableOrder().compare(xposEntry2, xposEntry2));
    assertEquals(0, Entry.<Integer>stableOrder().compare(xposEntryNull, xposEntryNull));
    assertEquals(0, Entry.<Integer>order().compare(xposEntry2, xposEntry2));
    assertEquals(0, Entry.<Integer>order().compare(xposEntryNull, xposEntryNull));

    // Equal instances for both order() and stableOrder().
    assertEquals(0, Entry.<Integer>stableOrder().compare(xposEntry2, xposEntry2b));
    assertEquals(0, Entry.<Integer>stableOrder().compare(xposEntryNull, xposEntryNullb));
    assertEquals(0, Entry.<Integer>order().compare(xposEntry2, xposEntry2b));
    assertEquals(0, Entry.<Integer>order().compare(xposEntryNull, xposEntryNullb));
  }

  @Test
  public void testPositioningIteratorEmpty() {
    reset();

    S2Iterator<Entry<Integer>> it = index.iterator();
    assertTrue(it.atBegin());
    assertTrue(it.done());

    assertFalse(it.next());
    assertFalse(it.prev());
  }

  @Test
  public void testPositioningIteratorSingleElement() {
    reset();
    add(data.getRandomPoint(), 0);
    verify();

    S2Iterator<Entry<Integer>> it = index.iterator();
    assertTrue(it.atBegin());
    assertFalse(it.done());
    Entry<Integer> e0 = it.entry();
    assertNotNull(e0);

    // prev() when atBegin() should return false, stay at element 0.
    assertFalse(it.prev());
    assertTrue(it.atBegin());
    assertFalse(it.done());
    assertSame(e0, it.entry());

    // next() from element 0 returns false, but makes done() true and entry() null.
    assertFalse(it.next());
    assertFalse(it.atBegin());
    assertTrue(it.done());
    assertNull(it.entry());

    // prev() from done() should move back to the single element.
    assertTrue(it.prev());
    assertTrue(it.atBegin());
    assertFalse(it.done());
    assertSame(it.entry(), e0);

    // prev() when atBegin() should return false, stay at element 0.
    assertFalse(it.prev());
    assertTrue(it.atBegin());
    assertFalse(it.done());
    assertSame(it.entry(), e0);
  }

  @Test
  public void testPositioningIteratorThreeElements() {
    reset();
    add(data.getRandomPoint(), 2);
    add(data.getRandomPoint(), 0);
    add(data.getRandomPoint(), 1);

    verify();

    S2Iterator<Entry<Integer>> it = index.iterator();
    assertTrue(it.atBegin());
    assertFalse(it.done());
    Entry<Integer> e0 = it.entry();
    assertNotNull(e0);

    // prev() when atBegin() should return false, stay at element 0.
    assertFalse(it.prev());
    assertTrue(it.atBegin());
    assertFalse(it.done());
    assertSame(it.entry(), e0);

    // next() from element 0 moves to element 1.
    assertTrue(it.next());
    assertFalse(it.atBegin());
    assertFalse(it.done());
    Entry<Integer> e1 = it.entry();
    assertSame(it.entry(), e1);

    // prev() from element 1 should move back to element 0.
    assertTrue(it.prev());
    assertTrue(it.atBegin());
    assertFalse(it.done());
    assertSame(it.entry(), e0);

    // next() twice takes us to the last element, but we are not done() yet.
    assertTrue(it.next());
    assertSame(it.entry(), e1);
    assertTrue(it.next());
    assertFalse(it.atBegin());
    assertFalse(it.done());
    Entry<Integer> e2 = it.entry();
    assertNotNull(e2);
    assertNotEquals(e2, e0);
    assertNotEquals(e2, e1);

    // prev() from element 2 should move back to element 1.
    assertTrue(it.prev());
    assertFalse(it.atBegin());
    assertFalse(it.done());
    assertSame(it.entry(), e1);

    // next to element 2
    assertTrue(it.next());
    assertFalse(it.atBegin());
    assertFalse(it.done());
    assertSame(it.entry(), e2);

    // next() when on the last element returns false, but makes done() true and entry() null.
    assertFalse(it.next());
    assertFalse(it.atBegin());
    assertTrue(it.done());
    assertNull(it.entry());

    // prev() from done() should move back to the last element.
    assertTrue(it.prev());
    assertFalse(it.atBegin());
    assertFalse(it.done());
    assertSame(it.entry(), e2);

    // finish() should move position "past the last" and make entry() return null.
    it.finish();
    assertFalse(it.atBegin());
    assertTrue(it.done());
    assertNull(it.entry());

    // prev() from finish() should move back to the last element.
    it.prev();
    assertFalse(it.atBegin());
    assertFalse(it.done());
    assertNotNull(it.entry());
    assertSame(it.entry(), e2);
  }

  /** Verifies the iterator methods. This is somewhat redundant with S2ShapeIndexTest. */
  public void checkIteratorMethods() {
    S2Iterator<Entry<Integer>> it = index.iterator();
    assertTrue(it.atBegin());
    if (index.numPoints() > 0) {
      assertNotNull(it.entry());
      assertFalse(it.done());
    } else {
      assertNull(it.entry());
      assertTrue(it.done());
    }

    // Try to go to the position after the first.
    boolean moved = it.next();
    if (index.numPoints() > 1) {
      assertTrue(moved);
      assertFalse(it.atBegin());
    } else {
      assertFalse(moved);
      assertTrue(it.atBegin());
    }
    if (index.numPoints() > 0) {
      assertNotNull(it.entry());
    } else {
      assertNull(it.entry());
    }

    // Go to the end of the index.
    it.finish();
    assertTrue(it.done());
    assertNull(it.entry());

    // Try to go to the position previous to the end.
    moved = it.prev();
    if (index.numPoints() < 2) {
      assertFalse(moved);
    } else {
      assertTrue(moved);
    }
    if (index.numPoints() > 0) {
      Entry<Integer> entry = it.entry();
      assertNotNull(entry);
    } else {
      assertNull(it.entry());
    }

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
        it2 = it.copy();
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

  @Test
  public void testNoPoints() {
    checkIteratorMethods();
  }

  @Test
  public void testEntryComparatorAndTreeSetContains() {
    Entry<String> e1 = S2PointIndex.createEntry(data.getRandomPoint(), "1");
    Entry<String> e2 = S2PointIndex.createEntry(data.getRandomPoint(), "2");
    Entry<String> e3 = S2PointIndex.createEntry(data.getRandomPoint(), "3");

    Comparator<Entry<String>> order = S2PointIndex.Entry.order();
    assertEquals(0, order.compare(e1, e1));
    assertEquals(0, order.compare(e2, e2));
    assertNotEquals(0, order.compare(e1, e2));
    assertNotEquals(0, order.compare(e2, e1));

    ObjectAVLTreeSet<Entry<String>> entries = new ObjectAVLTreeSet<>(order);
    entries.add(e1);
    entries.add(e2);

    assertTrue(entries.contains(e1));
    assertTrue(entries.contains(e2));
    assertFalse(entries.contains(e3));

    Entry<String> first = entries.first();
    assertTrue(entries.contains(first));
    Entry<String> last = entries.last();
    assertTrue(entries.contains(last));
  }

  @Test
  public void testRandomPoints() {
    for (int i = 0; i < 3; ++i) {
      add(data.getRandomPoint(), (int) data.uniform(0, 100));
    }
    verify();

    for (int i = 0; i < 1000; ++i) {
      add(data.getRandomPoint(), (int) data.uniform(0, 100));
    }
    verify();
    checkIteratorMethods();
  }

  @Test
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
