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
package com.google.common.geometry.primitives;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertThrows;
import static org.junit.Assert.assertTrue;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.primitives.Pullable.PullIterator;
import com.google.common.geometry.primitives.Pullable.PullList;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link PullCollections}. */
@RunWith(JUnit4.class)
public final class PullCollectionsTest {
  private static final int SEED = 123455;
  private final Random rand = new Random(SEED);

  /** A mutable long for testing. */
  private static class MutableLong {
    long val;

    public static MutableLong of(long value) {
      MutableLong m = new MutableLong();
      m.val = value;
      return m;
    }

    public long val() {
      return val;
    }
  }

  /**
   * A fixed-capacity PullList of MutableLong for testing the default implementations of PullList
   * methods.
   */
  private static class LongPullList implements PullList<MutableLong> {
    final long[] values;
    int size;

    private LongPullList() {
      values = new long[0];
      size = 0;
    }

    /** Constructs an empty LongPullList with the given capacity. */
    public LongPullList(int capacity) {
      values = new long[capacity];
      size = 0;
    }

    /** Constructs an LongPullList that is a copy of the given List of Longs. */
    public LongPullList(List<Long> list) {
      values = new long[list.size()];
      for (int i = 0; i < list.size(); i++) {
        values[i] = list.get(i);
      }
      size = list.size();
    }

    public static LongPullList of(long... values) {
      LongPullList result = new LongPullList(values.length);
      for (long value : values) {
        result.add(MutableLong.of(value));
      }
      return result;
    }

    private void rangeCheck(int index) {
      if (index < 0 || index >= size) {
        throw new IndexOutOfBoundsException("Index " + index + " out of bounds [0," + size + ")");
      }
    }

    @Override
    public MutableLong newElement() {
      return new MutableLong();
    }

    @Override
    public void get(int index, MutableLong element) {
      rangeCheck(index);
      element.val = values[index];
    }

    @Override
    public int size() {
      return size;
    }

    @Override
    public void clear() {
      size = 0;
    }

    @Override
    public void set(int index, MutableLong value) {
      rangeCheck(index);
      values[index] = value.val;
    }

    public void set(int index, long value) {
      rangeCheck(index);
      values[index] = value;
    }

    @Override
    public void add(MutableLong value) {
      ensureCapacity(size + 1);
      set(size++, value.val);
    }

    @Override
    @CanIgnoreReturnValue
    public boolean ensureCapacity(int capacity) {
      if (capacity <= values.length) {
        return false;
      }
      throw new UnsupportedOperationException("Increasing capacity not implemented");
    }

    @Override
    public void resize(int newSize) {
      ensureCapacity(newSize);
      size = newSize;
    }

    @Override
    public void truncate(int newSize) {
      if (newSize < size) {
        size = newSize;
      }
    }

    @Override
    public void enlarge(int newSize) {
      ensureCapacity(newSize);
      if (newSize > size) {
        size = newSize;
      }
    }

    /** A stable sort implementation built on asList(). Very inefficient. */
    public void stableSort(Comparator<MutableLong> comparator) {
      Collections.sort(asList(), comparator);
    }
  }

  private void fillConsecutive(PullList<MutableLong> pullList) {
    MutableLong element = new MutableLong();
    for (int i = 0; i < pullList.size(); i++) {
      element.val = i;
      pullList.set(i, element);
    }
  }

  /**
   * Verifies that the contents and size of the PullList of MutableLongs are the same as the given
   * "expected" ImmutableList of Longs by iterating the PullList.
   */
  public void assertPullList(ImmutableList<Long> expected, PullList<MutableLong> actual) {
    assertEquals(expected.size(), actual.size());
    int index = 0;

    MutableLong element = new MutableLong();
    for (Pullable.PullIterator it = actual.iterator(element); it.pull(); ) {
      long actualVal = element.val;
      long expectedVal = expected.get(index);
      assertEquals(
          "\nExpected:"
              + expectedVal
              + "\nActualVal  :"
              + actual
              + " at index "
              + index,
          expectedVal,
          actualVal);
      index++;
    }
  }

  /**
   * Verifies that the contents and size of the PullList of MutableLongs are the same as the given
   * "expected" array of longs by calling get() for each value.
   */
  public void assertPullList(long[] expected, PullList<MutableLong> actual) {
    assertEquals(expected.length, actual.size());

    MutableLong element = new MutableLong();
    for (int i = 0; i < expected.length; i++) {
      actual.get(i, element);
      long expectedVal = element.val;
      long actualVal = expected[i];
      assertEquals(
          "\nExpected:" + expectedVal + "\nActual  :" + actualVal + "\nMismatch at position " + i,
          expectedVal,
          actualVal);
    }
  }

  // Mostly tests the test class and test methods.
  @Test
  public void testConstruction() {
    // Construct a PullList<MutableLong> from an ImmutableList<Long>.
    ImmutableList<Long> list = ImmutableList.of(Long.valueOf(6), Long.valueOf(7), Long.valueOf(8));
    PullList<MutableLong> fromList = new LongPullList(list);
    assertPullList(list, fromList);

    assertPullList(new long[] {6, 7, 8}, fromList);
    assertFalse(fromList.isEmpty());

    // Construct a PullList<MutableLong> from given longs.
    PullList<MutableLong> fromArray = LongPullList.of(6, 7, 8);
    assertPullList(ImmutableList.of(Long.valueOf(6), Long.valueOf(7), Long.valueOf(8)), fromList);
    assertPullList(new long[] {6, 7, 8}, fromList);
    assertFalse(fromArray.isEmpty());

    // Construct a PullList<MutableLong> of a given capacity and then resize it to capacity.
    PullList<MutableLong> oneHundred = new LongPullList(100);
    assertEquals(0, oneHundred.size());
    assertPullList(new long[] {}, oneHundred);
    assertPullList(ImmutableList.of(), oneHundred);
    oneHundred.resize(100);
    assertEquals(100, oneHundred.size());
    assertFalse(oneHundred.isEmpty());
  }

  // Add, get, set, clear.
  @Test
  public void testBasicOperations() {
    PullList<MutableLong> pullList = new LongPullList(10);
    assertEquals(0, pullList.size());
    pullList.add(MutableLong.of(4));
    assertEquals(1, pullList.size());
    pullList.add(MutableLong.of(8));
    assertEquals(2, pullList.size());

    MutableLong element = new MutableLong();
    pullList.get(0, element);
    assertEquals(4, element.val);
    pullList.get(1, element);
    assertEquals(8, element.val);

    pullList.set(1, MutableLong.of(7));
    pullList.get(1, element);
    assertEquals(7, element.val);

    pullList.clear();
    assertEquals(0, pullList.size());
  }

  @Test
  public void testSwap() {
    PullList<MutableLong> pullList = LongPullList.of(new long[] {100, 101, 102, 103, 104});
    pullList.swap(0, 1);
    assertPullList(new long[] {101, 100, 102, 103, 104}, pullList);
    pullList.swap(3, 4);
    assertPullList(new long[] {101, 100, 102, 104, 103}, pullList);
    pullList.swap(0, 4);
    assertPullList(new long[] {103, 100, 102, 104, 101}, pullList);
  }

  @Test
  public void testResize() {
    // Capacity is 1000, but size is zero.
    PullList<MutableLong> pullList = new LongPullList(1000);
    assertEquals(0, pullList.size());

    pullList.resize(660);
    assertEquals(660, pullList.size());

    // Entries that haven't been set have default values.
    MutableLong element = new MutableLong();
    pullList.get(10, element);
    assertEquals(0, element.val);

    // Set every 25th element in the pullList to equal its index.
    for (int i = 0; i < 600; i += 25) {
      element.val = i;
      pullList.set(i, element);
    }

    // Increasing size shouldn't affect existing pullList contents.
    pullList.resize(800);
    assertEquals(800, pullList.size());
    for (int i = 0; i < 600; i += 25) {
      pullList.get(i, element);
      assertEquals(i, element.val);
    }

    // ensureCapacity only changes capacity when required, and doesn't change size.
    assertFalse(pullList.ensureCapacity(200));
    assertEquals(800, pullList.size());
    pullList.ensureCapacity(400); // Increases capacity to the next power of two.
    assertEquals(800, pullList.size());
    pullList.ensureCapacity(1000);
    assertEquals(800, pullList.size());

    // LongPullList doesn't support increasing capacity.
    assertThrows(UnsupportedOperationException.class, () -> pullList.ensureCapacity(1001));

    // enlarge() to a smaller value does nothing.
    pullList.enlarge(600);
    assertEquals(800, pullList.size());
    // The contents are unchanged.
    for (int i = 0; i < 600; i += 25) {
      pullList.get(i, element);
      assertEquals(i, element.val);
    }

    // truncate() to a larger value does nothing.
    pullList.truncate(2000);
    assertEquals(800, pullList.size());

    // Add another element, making the array size 801.
    element.val = 43;
    pullList.add(element);
    assertEquals(801, pullList.size());
    pullList.get(800, element);
    assertEquals(43, element.val);
  }

  @Test
  public void testOutOfBoundsExceptions() {
    PullList<MutableLong> pullList = new LongPullList(100);
    MutableLong element = new MutableLong();

    // get() on an empty pullList
    assertThrows(IndexOutOfBoundsException.class, () -> pullList.get(0, element));

    // set() on an empty pullList
    assertThrows(IndexOutOfBoundsException.class, () -> pullList.set(0, MutableLong.of(1)));

    // get() and set() past the last element.
    pullList.resize(20);
    assertEquals(20, pullList.size());

    assertThrows(IndexOutOfBoundsException.class, () -> pullList.set(20, MutableLong.of(1)));
    assertThrows(IndexOutOfBoundsException.class, () -> pullList.get(20, element));
  }

  public void assertSorted(PullList<MutableLong> pullList) {
    assertSorted(pullList, Comparator.comparingLong(MutableLong::val));
  }

  public void assertSorted(PullList<MutableLong> pullList, Comparator<MutableLong> comparator) {
    MutableLong curr = new MutableLong();
    MutableLong next = new MutableLong();
    for (int i = 0; i < pullList.size() - 1; i++) {
      pullList.get(i, curr);
      pullList.get(i + 1, next);
      assertTrue(
          "Expected value at "
              + i
              + " to compare be less than or equal to value at "
              + (i + 1)
              + " but values were "
              + curr.val
              + " and "
              + next.val
              + ".  PullList = "
              + pullList,
          comparator.compare(curr, next) <= 0);
    }
  }

  /** Compare MutableLongs by their last digit only. Useful for testing that sort is stable. */
  private static final Comparator<MutableLong> LAST_DIGIT_COMPARATOR =
      new Comparator<>() {
        @Override
        public int compare(MutableLong a, MutableLong b) {
          return Long.compare(a.val % 10, b.val % 10);
        }
      };

  @Test
  public void testSort() {
    for (int i = 0; i < 30; i++) {
      int listSize = 10 + rand.nextInt(10);
      PullList<MutableLong> pullList = new LongPullList(listSize);
      pullList.resize(listSize);
      fillConsecutive(pullList);
      Collections.shuffle(pullList.asList());
      pullList.sort(Comparator.comparingLong(MutableLong::val));
      assertSorted(pullList);
    }
  }

  @Test
  public void testSortByLastDigit() {
    for (int i = 0; i < 30; i++) {
      int listSize = 10 + rand.nextInt(10);
      PullList<MutableLong> pullList = new LongPullList(listSize);
      pullList.resize(listSize);
      fillConsecutive(pullList);
      Collections.shuffle(pullList.asList());
      pullList.sort(LAST_DIGIT_COMPARATOR);
      assertSorted(pullList, LAST_DIGIT_COMPARATOR);
    }
  }

  /**
   * Test asList() by verifying that stableSort is stable, i.e. does not change the relative order
   * of elements that compare as equal according to the comparator.
   */
  @Test
  public void testSortIsStable() {
    // The numbers ending in 2 should retain their relative ordering.
    LongPullList pullList = LongPullList.of(4, 3, 2, 1, 32, 22, 12);
    pullList.stableSort(LAST_DIGIT_COMPARATOR);
    assertPullList(new long[] {1, 2, 32, 22, 12, 3, 4}, pullList);

    // Numbers should be sorted into subranges by last digit while retaining their relative order.
    pullList = LongPullList.of(51, 53, 41, 43, 31, 33, 21, 23, 11, 13, 1, 3);
    pullList.stableSort(LAST_DIGIT_COMPARATOR);
    assertPullList(new long[] {51, 41, 31, 21, 11, 1, 53, 43, 33, 23, 13, 3}, pullList);
  }

  /** Verifies that the iterator visits all elements in the pullList. */
  @Test
  public void testIterator() {
    PullList<MutableLong> pullList = new LongPullList(10);
    pullList.resize(10);
    fillConsecutive(pullList);

    ArrayList<Long> visited = new ArrayList<>();
    MutableLong element = new MutableLong();
    for (PullIterator it = pullList.iterator(element); it.pull(); ) {
      visited.add(Long.valueOf(element.val));
    }

    assertEquals(10, visited.size());
    for (int i = 0; i < 10; i++) {
      assertEquals((long) i, visited.get(i).longValue());
    }
  }

  @Test
  public void testForEach() {
    PullList<MutableLong> pullList = new LongPullList(10);
    pullList.resize(3);
    fillConsecutive(pullList);

    ArrayList<Long> visited = new ArrayList<>();
    MutableLong element = new MutableLong();
    pullList.forEach(element, () -> visited.add(Long.valueOf(element.val)));

    assertEquals(3, visited.size());
    for (int i = 0; i < 3; i++) {
      assertEquals((long) i, visited.get(i).longValue());
    }

    // Single element list
    pullList.resize(1);
    visited.clear();
    pullList.forEach(element, () -> visited.add(Long.valueOf(element.val)));
    assertEquals(1, visited.size());
    assertEquals(0, visited.get(0).longValue());

    // Empty list
    pullList.clear();
    visited.clear();
    pullList.forEach(element, () -> visited.add(Long.valueOf(element.val)));
    assertEquals(0, visited.size());
  }
}
