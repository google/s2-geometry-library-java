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

import static java.util.Comparator.naturalOrder;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import com.google.common.geometry.MutableS2Point.MutableS2PointList;
import com.google.common.geometry.primitives.Pullable.PullIterator;
import java.util.ArrayList;
import java.util.List;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/**
 * Unit tests for {@link MutableS2Point} and {@link MutableS2PointList}. Also tests the
 * Pullable interface provided by MutableS2PointList.
 */
@RunWith(JUnit4.class)
public final class MutableS2PointTest extends GeometryTestCase {

  @Test
  public void testMutableS2PointImpl() {
    S2Point p123 = new S2Point(1, 2, 3);
    S2Point p234 = new S2Point(2, 3, 4);
    MutableS2Point.MutableS2PointImpl mutablePoint = new MutableS2Point.MutableS2PointImpl();
    assertFalse(mutablePoint.isEqualTo(p123));
    mutablePoint.set(p123);
    assertTrue(mutablePoint.isEqualTo(p123));
    assertFalse(mutablePoint.isEqualTo(p234));
    mutablePoint.set(2, 3, 4);
    assertTrue(mutablePoint.isEqualTo(p234));
    assertFalse(mutablePoint.isEqualTo(p123));
    assertTrue(mutablePoint.getX() == 2);
    assertTrue(mutablePoint.getY() == 3);
    assertTrue(mutablePoint.getZ() == 4);
  }

  @Test
  public void testMutableS2PointCompare() {
    MutableS2PointList pair = MutableS2PointList.pair();

    // Test where point 1 is greater than point 0, comparison only needs to consider X.
    pair.set(0, new S2Point(1, 2, 3));
    pair.set(1, new S2Point(2, 3, 4));

    assertTrue(pair.get(0).compareTo(pair.get(1)) < 0);
    assertTrue(pair.get(1).compareTo(pair.get(0)) > 0);
    assertEquals(0, pair.get(0).compareTo(pair.get(0)));

    // Same comparisons, but on the list directly
    assertTrue(pair.compare(0, 1) < 0);
    assertTrue(pair.compare(1, 0) > 0);
    assertEquals(0, pair.compare(0, 0));

    // Test where point 1 is less than point 0, with equal X, tie-breaking on Y.
    pair.set(1, new S2Point(1, 1, 1));
    assertTrue(pair.get(0).compareTo(pair.get(1)) > 0);
    assertTrue(pair.get(1).compareTo(pair.get(0)) < 0);
    assertEquals(0, pair.get(1).compareTo(pair.get(1)));

    // Same comparisons, but on the list directly
    assertTrue(pair.compare(0, 1) > 0);
    assertTrue(pair.compare(1, 0) < 0);
    assertEquals(0, pair.compare(1, 1));

    // Test tie-break on Z.
    pair.set(0, new S2Point(1, 1, 0));
    assertTrue(pair.get(0).compareTo(pair.get(1)) < 0);
    assertTrue(pair.get(1).compareTo(pair.get(0)) > 0);

    // Same comparisons, but on the list directly
    assertTrue(pair.compare(0, 1) < 0);
    assertTrue(pair.compare(1, 0) > 0);
  }

  @Test
  public void testListGetAndSet() {
    MutableS2PointList list = MutableS2PointList.pair();
    assertEquals(2, list.size());
    list.set(0, new S2Point(1, 2, 3));
    list.set(1, new S2Point(2, 3, 4));

    MutableS2Point mutablePoint = list.newElement();
    list.get(0, mutablePoint);
    assertTrue(mutablePoint.isEqualTo(new S2Point(1, 2, 3)));
    list.get(1, mutablePoint);
    assertTrue(mutablePoint.isEqualTo(new S2Point(2, 3, 4)));

    // Test the MutableS2Point interface provided by get().
    MutableS2Point inList = list.get(0);
    assertTrue(inList.isEqualTo(new S2Point(1, 2, 3)));
    inList.set(2, 3, 4);
    assertTrue(inList.isEqualTo(list.get(1)));
    inList.set(new S2Point(1, 0, 1));
    assertTrue(list.get(0).isEqualTo(new S2Point(1, 0, 1)));
    inList.set(new S2Point(2, 3, 4));
    assertTrue(list.isEqualTo(0, new S2Point(2, 3, 4)));
    list.set(1, inList);
    assertTrue(list.isEqualTo(1, new S2Point(2, 3, 4)));
  }

  @Test
  public void testMoreListOperations() {
    MutableS2PointList list = new MutableS2PointList();
    assertTrue(list.isEmpty());
    list.add(1, 2, 3);
    list.add(2, 3, 4);
    assertFalse(list.isEmpty());

    // Get a list entry as an S2Point and check it.
    S2Point immutablePoint = list.getImmutable(1);
    S2Point expectedPoint = new S2Point(2, 3, 4);
    assertTrue(immutablePoint.equalsPoint(expectedPoint));

    MutableS2Point mutablePoint = list.newElement();

    // Access to entries outside of the range throws an exception.
    assertThrows(IndexOutOfBoundsException.class, () -> list.get(2));
    assertThrows(IndexOutOfBoundsException.class, () -> list.getImmutable(2));
    assertThrows(IndexOutOfBoundsException.class, () -> list.set(2, new S2Point(1, 2, 3)));
    assertThrows(IndexOutOfBoundsException.class, () -> list.get(2, mutablePoint));
    assertThrows(IndexOutOfBoundsException.class, () -> list.getImmutable(2));

    list.add(new S2Point(3, 4, 5));
    assertEquals(3, list.size());
    assertTrue(list.get(2).isEqualTo(new S2Point(3, 4, 5)));

    // Test the various isEqualTo methods
    MutableS2Point p0 = list.newElement();
    list.get(0, p0);
    assertTrue(list.isEqualTo(0, p0));
    assertFalse(list.isEqualTo(1, p0));

    list.get(0, mutablePoint);
    assertTrue(mutablePoint.isEqualTo(p0));

    assertTrue(list.isEqualTo(0, 0));
    assertFalse(list.isEqualTo(1, 0));

    assertFalse(list.isEqualTo(2, p0));
    assertFalse(list.isEqualTo(0, 2));
    list.copy(0, 2);
    assertTrue(list.isEqualTo(2, p0));
    assertTrue(list.isEqualTo(0, 2));
  }

  @Test
  public void testAccessInvalidMutableS2Point() {
    MutableS2PointList list = new MutableS2PointList(3);
    list.set(0, new S2Point(0, 0, 1));
    list.set(1, new S2Point(0, 0, 1));
    list.set(2, new S2Point(0, 0, 2));

    MutableS2Point inList = list.get(0);
    assertTrue(inList.isEqualTo(new S2Point(0, 0, 1)));

    // After clear(), the MutableS2Point 'inList' is no longer valid, although that's not checked.
    list.clear();
    assertEquals(0, list.size());
    // inList is not actually "in" the list any more, but the array is still there.
    assertTrue(inList.isEqualTo(new S2Point(0, 0, 1)));
    // Add a point (as coordinates) at index 0, which is what inList is pointing at.
    list.add(1, 1, 1);
    // Now inList is valid again, but has different contents.
    assertTrue(inList.isEqualTo(new S2Point(1, 1, 1)));
  }

  @Test
  public void testCopyConstructor() {
    // Size zero, default capacity.
    MutableS2PointList l1 = new MutableS2PointList();
    l1.add(new S2Point(1, 2, 3));
    l1.add(new S2Point(2, 3, 4));
    assertEquals(2, l1.size());
    assertEquals(MutableS2PointList.DEFAULT_CAPACITY, l1.capacity());

    // Copy constructor.
    MutableS2PointList l2 = new MutableS2PointList(l1);

    // Size and capacity of the copied list should both be 2.
    assertEquals(2, l2.size());
    assertEquals(2, l2.capacity());
    // And contents should be the same.
    assertTrue(l2.get(0).isEqualTo(new S2Point(1, 2, 3)));
    assertTrue(l2.get(1).isEqualTo(new S2Point(2, 3, 4)));

    // Enlarge the copied list, and verify that the capacity increases to the next power of two.
    l2.enlarge(5);
    assertEquals(5, l2.size());
    assertEquals(8, l2.capacity());
  }

  @Test
  public void testAsMutablePointList() {
    MutableS2PointList mutablePointList = new MutableS2PointList();
    mutablePointList.add(1, 2, 3);
    mutablePointList.add(2, 3, 4);
    mutablePointList.add(3, 4, 5);

    List<MutableS2Point> asList = mutablePointList.asList();
    assertEquals(3, asList.size());
    assertTrue(asList.get(0).isEqualTo(new S2Point(1, 2, 3)));
    assertTrue(asList.get(1).isEqualTo(new S2Point(2, 3, 4)));
    assertTrue(asList.get(2).isEqualTo(new S2Point(3, 4, 5)));

    // Verify that set() returns the previous value.
    MutableS2Point newValue = mutablePointList.newElement();
    newValue.set(4, 5, 6);
    assertTrue(asList.set(1, newValue).isEqualTo(new S2Point(2, 3, 4)));
    assertTrue(mutablePointList.get(1).isEqualTo(new S2Point(4, 5, 6)));
  }

  @Test
  public void testAsPointList() {
    MutableS2PointList mutablePointList = new MutableS2PointList();
    mutablePointList.add(1, 2, 3);
    mutablePointList.add(2, 3, 4);
    mutablePointList.add(3, 4, 5);

    List<S2Point> asPointList = mutablePointList.asPointList();
    assertEquals(3, asPointList.size());
    assertTrue(asPointList.get(0).equalsPoint(new S2Point(1, 2, 3)));
    assertTrue(asPointList.get(1).equalsPoint(new S2Point(2, 3, 4)));
    assertTrue(asPointList.get(2).equalsPoint(new S2Point(3, 4, 5)));

    // Verify that set() returns the previous value.
    assertTrue(asPointList.set(1, new S2Point(4, 5, 6)).equalsPoint(new S2Point(2, 3, 4)));
    assertTrue(mutablePointList.get(1).isEqualTo(new S2Point(4, 5, 6)));
  }

  @Test
  public void testPullIterator() {
    MutableS2PointList mutablePointList = new MutableS2PointList();
    mutablePointList.add(1, 2, 3);
    mutablePointList.add(2, 3, 4);
    mutablePointList.add(3, 4, 5);

    MutableS2Point mutablePoint = mutablePointList.newElement();
    PullIterator it = mutablePointList.iterator(mutablePoint);

    // The mutablePoint isn't set until pull() is called, so has uninitialized contents.
    assertTrue(mutablePoint.isEqualTo(new S2Point(0, 0, 0)));
    assertTrue(it.pull());
    assertTrue(mutablePoint.isEqualTo(new S2Point(1, 2, 3)));
    assertTrue(it.pull());
    assertTrue(mutablePoint.isEqualTo(new S2Point(2, 3, 4)));
    assertTrue(it.pull());
    assertTrue(mutablePoint.isEqualTo(new S2Point(3, 4, 5)));
    assertFalse(it.pull());
  }

  @Test
  public void testMutablePointListForEach() {
    // Set up a list of 10 points, (0,0,0), (1,0,0), ... (9,0,0).
    MutableS2PointList list = new MutableS2PointList();
    for (int i = 0; i < 10; i++) {
      list.add(i, 0, 0);
    }

    int[] counter = { 0 };

    // Visit all the points as coordinates.
    assertTrue(list.forEach((x, y, z) -> {
      // Verify that entries are visited in order.
      int xi = (int) Math.round(x);
      assertEquals(counter[0], xi);
      counter[0]++;
      return true; // Don't stop visiting.
    }));
    assertEquals(10, counter[0]);

    // Reset the counter and visit some of the points as coordinates, stopping early.
    counter[0] = 0;
    assertFalse(list.forEach((x, y, z) -> {
      // Verify that entries are visited in order.
      int xi = (int) Math.round(x);
      assertEquals(counter[0], xi);
      counter[0]++;
      return (xi != 8); // Stop visiting when we see the ninth point, with x=8.
    }));
    assertEquals(9, counter[0]);

    // Visit all the points as index and coordinates.
    counter[0] = 0;
    assertTrue(list.forEach((i, x, y, z) -> {
      // Verify that entries are visited in order.
      int xi = (int) Math.round(x);
      assertEquals(counter[0], xi);
      assertEquals(xi, i);
      counter[0]++;
      return true; // Don't stop visiting.
    }));
    assertEquals(10, counter[0]);

    // Same again but stop early.
    counter[0] = 0;
    assertFalse(list.forEach((i, x, y, z) -> {
      // Verify that entries are visited in order.
      int xi = (int) Math.round(x);
      assertEquals(counter[0], xi);
      assertEquals(xi, i);
      counter[0]++;
      return (xi != 8); // Stop visiting when we see the ninth point, with x=8.
    }));
    assertEquals(9, counter[0]);
  }

  /**
   * Tests construction of 100 random sequences of random points, with random length and various
   * initial capacities, followed by a truncation and modification of the contents.
   */
  @Test
  public void testRandomSequences() {
    for (int i = 0; i < 100; i++) {
      int length = data.nextInt(500);
      ArrayList<S2Point> pointArrayList = new ArrayList<>(length);

      // Test with various initial capacities, but always with initial size zero.
      MutableS2PointList mutablePointList;
      int capacityCase = data.nextInt(4);
      switch (capacityCase) {
        case 0:
          mutablePointList = new MutableS2PointList(); // Default capacity.
          break;
        case 1:
          mutablePointList = new MutableS2PointList(0); // Capacity zero.
          break;
        case 2:
          mutablePointList = MutableS2PointList.ofCapacity(length); // Capacity == length.
          break;
        case 3:
          mutablePointList = MutableS2PointList.ofCapacity(length + 16); // Capacity > length.
          break;
        default:
          throw new AssertionError("Unexpected case: " + data.nextInt(3));
      }

      // Add 'length' random points to the list.
      for (int j = 0; j < length; j++) {
        S2Point point = data.getRandomPoint();
        pointArrayList.add(point);
        mutablePointList.add(point);
      }
      assertEqualLists(pointArrayList, mutablePointList);

      // Truncate both lists. Might be a no-op if they're already shorter.
      if (pointArrayList.size() > 100) {
        pointArrayList.subList(100, pointArrayList.size()).clear();
      }
      mutablePointList.truncate(100);
      assertEqualLists(pointArrayList, mutablePointList);

      // Modify the contents of both lists.
      for (int j = 0; j < pointArrayList.size(); j++) {
        S2Point point = data.getRandomPoint();
        pointArrayList.set(j, point);
        mutablePointList.set(j, point);
      }
      assertEqualLists(pointArrayList, mutablePointList);

      // Sort both lists.
      pointArrayList.sort(naturalOrder());
      mutablePointList.sort();
      assertEqualLists(pointArrayList, mutablePointList);
    }
  }

  private static void assertEqualLists(
      ArrayList<S2Point> expected, MutableS2PointList actual) {
    assertEquals(expected.size(), actual.size());
    for (int i = 0; i < expected.size(); i++) {
      assertTrue(actual.isEqualTo(i, expected.get(i)));
    }
  }
}