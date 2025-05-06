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
package com.google.common.geometry.primitives;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import com.google.common.collect.ImmutableList;
import com.google.common.geometry.primitives.Ints.ImmutableIntSequence;
import com.google.common.geometry.primitives.Ints.IntComparator;
import com.google.common.geometry.primitives.Ints.IntConsumer;
import com.google.common.geometry.primitives.Ints.IntSequence;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Random;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for IntVector. */
@RunWith(JUnit4.class)
public final class IntVectorTest {
  private static final int SEED = 123455;
  private final Random rand = new Random(SEED);

  private static final IntComparator LEXICOGRAPHIC_COMPARATOR =
      new IntComparator() {
        @Override
        public int compare(int a, int b) {
          return Integer.toString(a).compareTo(Integer.toString(b));
        }
      };

  /**
   * Verifies that the contents and size of the IntVector are the same as the given "expected"
   * ImmutableList of Integers by iterating the IntVector with forEach().
   */
  public void assertIntVector(ImmutableList<Integer> expected, IntVector actual) {
    assertEquals(expected.size(), actual.size());
    Iterator<Integer> iter = expected.iterator();
    actual.forEach((IntConsumer) element ->
        assertEquals(
            "\nExpected:"
                + IntVector.copyOf(expected)
                + "\nActual  :"
                + actual,
            iter.next().intValue(),
            element));
  }

  /**
   * Verifies that the contents and size of the IntVector are the same as the given "expected" array
   * of ints by calling get() for each value.
   */
  public void assertIntVector(int[] expected, IntVector actual) {
    assertEquals(expected.length, actual.size());
    for (int i = 0; i < expected.length; i++) {
      assertEquals(
          "\nExpected:"
              + IntVector.copyOf(expected)
              + "\nActual  :"
              + actual
              + "\nMismatch at position "
              + i,
          expected[i],
          actual.get(i));
    }
  }

  public ImmutableList<Integer> generateRandomImmutableList(int length) {
    ImmutableList.Builder<Integer> b = new ImmutableList.Builder<>();
    for (int i = 0; i < length; i++) {
      b.add(rand.nextInt());
    }
    return b.build();
  }

  /**
   * Tests construction of 100 random sequences of random length up to 9000, from copyOf(List) and
   * with the copy constructor.
   */
  @Test
  public void testRandomSequences() {
    for (int i = 0; i < 100; i++) {
      ImmutableList<Integer> list = generateRandomImmutableList(rand.nextInt(9000));
      // Verify construction from a List<Integer>.
      IntVector vector = IntVector.copyOf(list);
      assertIntVector(list, vector);
      // Verify the copy constructor.
      IntVector copied = new IntVector(vector);
      vector.clear();
      assertIntVector(list, copied);
    }
  }

  @Test
  public void testConstructors() {
    IntVector empty = new IntVector();
    assertIntVector(ImmutableList.of(), empty);
    assertIntVector(new int[] {}, empty);
    assertTrue(empty.isEmpty());

    // Construct an IntVector from a List<Integer>.
    IntVector fromList = IntVector.copyOf(ImmutableList.of(6, 7, 8));
    assertIntVector(ImmutableList.of(6, 7, 8), fromList);
    assertIntVector(new int[] {6, 7, 8}, fromList);
    assertFalse(fromList.isEmpty());

    // Construct an IntVector of a given size.
    IntVector oneHundred = IntVector.ofSize(100);
    assertEquals(100, oneHundred.size());
    assertEquals(128, oneHundred.capacity());
    assertFalse(oneHundred.isEmpty());

    // Exactly the requested size is allocated for a power of two.
    IntVector oneK = IntVector.ofSize(1024);
    assertEquals(1024, oneK.size());
    assertEquals(1024, oneK.capacity());

    // 40K of data will allocate 64K.
    IntVector fortyK = IntVector.ofSize(1024 * 40);
    assertEquals(40960, fortyK.size());
    assertEquals(65536, fortyK.capacity());

    // Constructing an IntVector with a specified capacity and zero size.
    IntVector oneKCapacity = IntVector.ofCapacity(1024);
    assertEquals(0, oneKCapacity.size());
    assertEquals(1024, oneKCapacity.capacity());

    // Capacity rounds up to a power of two.
    IntVector twoThousandCapacity = IntVector.ofCapacity(2000);
    assertEquals(0, twoThousandCapacity.size());
    assertEquals(2048, twoThousandCapacity.capacity());
  }

  @Test
  public void testBasicOperations() {
    IntVector vector = IntVector.empty();
    assertEquals(0, vector.size());
    // Push increases the vector size by one.
    vector.push(4);
    assertEquals(1, vector.size());
    vector.push(8);
    assertEquals(2, vector.size());

    // Mixing get/set with push/pop is probably bad practice but is useful for testing. Pushed
    // elements are at the end of the vector.
    assertEquals(4, vector.get(0));
    assertEquals(8, vector.get(1));
    assertEquals(8, vector.peek());

    vector.set(1, 7);
    assertEquals(7, vector.peek());
    assertEquals(2, vector.size());
    assertEquals(7, vector.pop());
    assertEquals(1, vector.size());

    // Clear resets the vector size to zero.
    vector.clear();
    assertEquals(0, vector.size());
  }

  @Test
  public void testFill() {
    IntVector smallVector = IntVector.ofSize(12);
    assertEquals(12, smallVector.size());
    smallVector.fill(99);
    assertEquals(99, smallVector.get(0));
    assertEquals(smallVector.toString(), 99, smallVector.get(11));

    IntVector bigVector = IntVector.ofSize(12000);
    bigVector.fill(99);
    assertEquals(99, bigVector.get(0));
    assertEquals(99, bigVector.get(2047));
    assertEquals(99, bigVector.get(2048));
    assertEquals(99, bigVector.get(11999));
  }

  @Test
  public void testFillConsecutive() {
    IntVector smallVector = IntVector.ofSize(12);
    smallVector.fillConsecutive();
    assertEquals(0, smallVector.get(0));
    assertEquals(11, smallVector.get(11));
    assertEquals(12, smallVector.size());

    IntVector bigVector = IntVector.ofSize(12000);
    bigVector.fillConsecutive();
    assertEquals(0, bigVector.get(0));
    assertEquals(2047, bigVector.get(2047));
    assertEquals(2048, bigVector.get(2048));
    assertEquals(11999, bigVector.get(11999));
    assertEquals(12000, bigVector.size());
  }

  @Test
  public void testReverse() {
    IntVector primes = IntVector.of(1, 3, 5, 7, 11, 13);
    primes.reverse();
    assertIntVector(new int[] {13, 11, 7, 5, 3, 1}, primes);
  }

  @Test
  public void testSwap() {
    IntVector vector = IntVector.copyOf(new int[] {100, 101, 102, 103, 104});
    vector.swap(0, 1);
    assertIntVector(new int[] {101, 100, 102, 103, 104}, vector);
    vector.swap(3, 4);
    assertIntVector(new int[] {101, 100, 102, 104, 103}, vector);
    vector.swap(0, 4);
    assertIntVector(new int[] {103, 100, 102, 104, 101}, vector);
  }

  @Test
  public void testResizeAndReallocate() {
    IntVector vector = new IntVector();
    assertEquals(0, vector.size());

    // The default is for 16 blocks of size 1024 each, so increasing size to 16,000 won't require
    // reallocating the base array, just allocating blocks.
    vector.resize(16000);
    assertEquals(16000, vector.size());
    assertEquals(16384, vector.capacity());

    // Vector entries that haven't been set have value zero.
    assertEquals(0, vector.get(100));
    assertEquals(0, vector.get(10000));

    // Set every 250th element in the vector to equal its index.
    for (int i = 0; i < vector.size(); i += 250) {
      vector.set(i, i);
    }

    // Increasing size shouldn't affect existing vector contents.
    vector.resize(18000);
    assertEquals(18000, vector.size());
    assertEquals(32768, vector.capacity());
    for (int i = 0; i < 16000; i += 250) {
      assertEquals(i, vector.get(i));
    }

    // The contents of the added elements are zero.
    for (int i = 16000; i < 18000; i += 250) {
      assertEquals(0, vector.get(i));
    }

    // ensureCapacity only changes capacity when required, and will not decrease capacity.
    vector.ensureCapacity(20000);
    assertEquals(32768, vector.capacity());
    assertEquals(18000, vector.size());
    vector.ensureCapacity(30000); // Does not change capacity.
    assertEquals(32768, vector.capacity());
    assertEquals(18000, vector.size());
    vector.shrink(); // Does not change capacity.
    assertEquals(32768, vector.capacity());
    assertEquals(18000, vector.size());
    vector.ensureCapacity(40000); // Increases capacity to the next power of two.
    assertEquals(65536, vector.capacity());
    assertEquals(18000, vector.size());

    vector.ensureCapacity(30000); // Does not decrease capacity.
    assertEquals(65536, vector.capacity());
    assertEquals(18000, vector.size());
    vector.shrink(); // Does reduce capacity.
    assertEquals(32768, vector.capacity());
    assertEquals(18000, vector.size());
    vector.ensureCapacity(16000); // Does not decrease capacity because 16384 < 18000.
    assertEquals(32768, vector.capacity());
    assertEquals(18000, vector.size());


    // enlarge() to a smaller value does nothing.
    vector.enlarge(8000);
    // The size does not change.
    assertEquals(18000, vector.size());
    // The existing contents are unchanged and are still accessible.
    for (int i = 0; i < 16000; i += 1000) {
      assertEquals(i, vector.get(i));
    }

    // truncate() to a larger value does nothing.
    vector.truncate(20000);
    // The size does not change.
    assertEquals(18000, vector.size());
    // The existing contents are unchanged and are still accessible.
    for (int i = 0; i < 16000; i += 1000) {
      assertEquals(i, vector.get(i));
    }

    // Make the vector so large that reallocating the array will be required.
    vector.resize(1_000_000);
    assertEquals(1_000_000, vector.size());
    vector.set(999999, 42);

    // Check that the existing vector contents were not changed or lost during realloc.
    for (int i = 0; i < 16000; i += 250) {
      assertEquals(i, vector.get(i));
    }

    // Push another element, making the array size a million and one.
    vector.push(43);
    assertEquals(1_000_001, vector.size());
    assertEquals(1_048_576, vector.capacity());
    assertEquals(42, vector.get(999999));
    assertEquals(43, vector.get(1000000));

    // truncate() followed by enlarge() clears the content, although reuses storage.
    vector.truncate(8000);
    vector.enlarge(16000);
    // The not-truncated part is unchanged,
    for (int i = 0; i < 8000; i += 1000) {
      assertEquals(i, vector.get(i));
    }
    // The truncated & enlarged part is zeroed out.
    for (int i = 8000; i < 16000; i += 1000) {
      assertEquals(0, vector.get(i));
    }

    // Similarly, clear() followed by enlarge() zeros out content.
    vector.clear();
    vector.enlarge(8000);
    for (int i = 0; i < 8000; i += 1000) {
      assertEquals(0, vector.get(i));
    }
  }

  @Test
  public void testTruncateAndShrink() {
    IntVector vector = new IntVector();
    vector.resize(80000);
    assertEquals(80000, vector.size());
    assertEquals(131072, vector.capacity());
    // shrink() does nothing here, it cannot reduce the data size to a smaller power of two.
    vector.shrink();
    assertEquals(80000, vector.size());
    assertEquals(131072, vector.capacity());

    vector.set(19999, 42);
    vector.truncate(20000);
    assertEquals(20000, vector.size());
    // truncate() does not shrink the underlying data array.
    assertEquals(131072, vector.capacity());

    // But now, calling shrink() reduces the length of the data array to a smaller power of two.
    vector.shrink();
    assertEquals(20000, vector.size());
    assertEquals(32768, vector.capacity());
    assertEquals(42, vector.get(19999));
  }

  @Test
  public void testSubList() {
    IntVector vector = IntVector.ofSize(20000);
    vector.fillConsecutive();

    IntSequence startRange = vector.subList(0, 3);
    assertEquals(3, startRange.size());
    OfInt iter = startRange.intIterator();
    assertTrue(iter.hasNext());
    assertEquals(0, iter.nextInt());
    assertTrue(iter.hasNext());
    assertEquals(1, iter.nextInt());
    assertTrue(iter.hasNext());
    assertEquals(2, iter.nextInt());
    assertFalse(iter.hasNext());

    IntSequence emptyRange = vector.subList(1, 1);
    assertEquals(0, emptyRange.size());
    iter = emptyRange.intIterator();
    assertFalse(iter.hasNext());

    IntSequence endRange = vector.subList(19995, 20000);
    assertEquals(5, endRange.size());
    iter = endRange.intIterator();
    assertTrue(iter.hasNext());
    assertEquals(19995, iter.nextInt());
    assertTrue(iter.hasNext());
    assertEquals(19996, iter.nextInt());
    assertTrue(iter.hasNext());
    assertEquals(19997, iter.nextInt());
    assertTrue(iter.hasNext());
    assertEquals(19998, iter.nextInt());
    assertTrue(iter.hasNext());
    assertEquals(19999, iter.nextInt());
    assertFalse(iter.hasNext());
  }

  @Test
  public void testToArray() {
    IntVector vector = IntVector.ofSize(20);
    vector.fillConsecutive();
    int[] a = vector.toArray();
    assertEquals(20, a.length);
    assertEquals(0, a[0]);
    assertEquals(19, a[19]);
  }

  @Test
  public void testIncrementAndDecrementAt() {
    IntVector vector = IntVector.ofSize(20);
    vector.fillConsecutive();

    assertEquals(0, vector.get(0));
    vector.incrementAt(0);
    assertEquals(1, vector.get(0));
    vector.decrementAt(0);
    assertEquals(0, vector.get(0));

    assertEquals(19, vector.get(19));
    vector.incrementAt(19);
    assertEquals(20, vector.get(19));
    vector.decrementAt(19);
    assertEquals(19, vector.get(19));
  }

  @Test
  public void testOutOfBoundsExceptions() {
    IntVector vector = IntVector.empty();

    // get() on an empty vector
    try {
      int unused = vector.get(0);
      fail();
    } catch (IndexOutOfBoundsException expectedException) {
      // Expected.
    }

    // set() on an empty vector
    try {
      vector.set(0, 1);
      fail();
    } catch (IndexOutOfBoundsException expectedException) {
      // Expected.
    }

    // get() and set() past the last element.
    vector.resize(2000);
    assertEquals(2000, vector.size());
    vector.set(1999, 53);
    assertEquals(53, vector.get(1999));
    try {
      vector.set(2000, 1);
      fail();
    } catch (ArrayIndexOutOfBoundsException expectedException) {
      // Expected.
    }
    try {
      int unused = vector.get(2000);
      fail();
    } catch (ArrayIndexOutOfBoundsException expectedException) {
      // Expected.
    }

    // subList() with fromIndex below 0.
    try {
      IntSequence unused = vector.subList(-1, 10);
      fail();
    } catch (IndexOutOfBoundsException expectedException) {
      // Expected.
    }

    // subList() with fromIndex above vector.size().
    try {
      IntSequence unused = vector.subList(2000, 2010);
      fail();
    } catch (IndexOutOfBoundsException expectedException) {
      // Expected.
    }

    // subList() with fromIndex above toIndex.
    try {
      IntSequence unused = vector.subList(100, 99);
      fail();
    } catch (IllegalArgumentException expectedException) {
      // Expected.
    }

    // subList() with toIndex past length.
    try {
      IntSequence unused = vector.subList(100, 2001);
      fail();
    } catch (IndexOutOfBoundsException expectedException) {
      // Expected.
    }

    // incrementAt() with index below zero.
    try {
      vector.incrementAt(-1);
      fail();
    } catch (IndexOutOfBoundsException expectedException) {
      // Expected.
    }
    // decrementAt() with index below zero.
    try {
      vector.decrementAt(-1);
      fail();
    } catch (IndexOutOfBoundsException expectedException) {
      // Expected.
    }
    // incrementAt() with index past end.
    try {
      vector.incrementAt(2000);
      fail();
    } catch (IndexOutOfBoundsException expectedException) {
      // Expected.
    }
    // incrementAt() with index below zero.
    try {
      vector.decrementAt(2000);
      fail();
    } catch (IndexOutOfBoundsException expectedException) {
      // Expected.
    }
  }

  public void assertSorted(IntVector vector) {
    for (int i = 0; i < vector.size() - 1; i++) {
      int curr = vector.get(i);
      int next = vector.get(i + 1);
      assertTrue(
          "Expected value at " + i + " to be less than or equal to value at " + (i + 1)
              + " but values were " + curr + " and " + next + ".  Vector = " + vector,
          curr <= next);
    }
  }

  /**
   * This is not very sophisticated test coverage, but the sort method is a one-liner with the
   * actual logic tested elsewhere. This method also tests the assertSorted() call above and {@link
   * IntVector#isSorted()} methods.
   */
  @Test
  public void testSort() {
    for (int i = 0; i < 30; i++) {
      IntVector vector = IntVector.ofSize(10 + rand.nextInt(10));
      vector.fillConsecutive();
      vector.shuffle(rand);
      assertFalse("After shuffle: " + vector, vector.isSorted());
      vector.sort();
      assertTrue("After sort: " + vector, vector.isSorted());
      assertSorted(vector);
    }
  }

  /** Compare numbers by their last digit only. Useful for testing that sort is stable. */
  private static final IntComparator LAST_DIGIT_COMPARATOR =
      new IntComparator() {
        @Override
        public int compare(int a, int b) {
          return Integer.compare(a % 10, b % 10);
        }
      };

  /**
   * Test that sort is stable, i.e. does not change the relative order of elements that compare as
   * equal according to the comparator.
   */
  @Test
  public void testSortIsStable() {
    // The numbers ending in 2 should retain their relative ordering.
    IntVector vector = IntVector.of(4, 3, 2, 1, 32, 22, 12);
    vector.sort(LAST_DIGIT_COMPARATOR);
    assertIntVector(new int[] {1, 2, 32, 22, 12, 3, 4}, vector);

    // Numbers should be sorted into subranges by last digit while retaining their relative order.
    vector = IntVector.of(51, 53, 41, 43, 31, 33, 21, 23, 11, 13, 1, 3);
    vector.sort(LAST_DIGIT_COMPARATOR);
    assertIntVector(new int[] {51, 41, 31, 21, 11, 1, 53, 43, 33, 23, 13, 3}, vector);
  }

  @SuppressWarnings("JUnitIncompatibleType")
  @Test
  public void testUnique() {
    IntVector empty = new IntVector();
    empty.unique();
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(), empty);

    IntVector single = IntVector.of(1);
    single.unique();
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1), single);

    IntVector duplicate = IntVector.of(1, 1);
    duplicate.unique();
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1), duplicate);

    IntVector unique = IntVector.of(1, 2);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2), unique);
    unique.unique();
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2), unique);

    IntVector dupAtEnd = IntVector.of(1, 2, 3, 4, 4);
    dupAtEnd.unique();
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2, 3, 4), dupAtEnd);

    IntVector dupAtStart = IntVector.of(3, 3, 4, 5, 6);
    dupAtStart.unique();
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(3, 4, 5, 6), dupAtStart);

    IntVector dupInMiddle = IntVector.of(3, 4, 5, 5, 6);
    dupInMiddle.unique();
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(3, 4, 5, 6), dupInMiddle);
  }

  /**
   * Tests the upperBound() and lowerBound() methods taking a custom IntComparator. upperBound()
   * returns the lowest index containing an element strictly greater than the given key, while
   * lowerBound() returns the lowest index containing an element greater than or equal to the given
   * key (according to the provided comparator).
   */
  @Test
  public void testUpperBoundWithIntComparator() {
    IntVector empty = new IntVector();
    assertEquals(0, empty.lowerBound(1, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(0, empty.upperBound(1, LEXICOGRAPHIC_COMPARATOR));

    IntVector single = IntVector.of(4);
    assertEquals(0, single.lowerBound(1, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(0, single.upperBound(1, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(0, single.lowerBound(4, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(1, single.upperBound(4, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(1, single.lowerBound(5, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(1, single.upperBound(5, LEXICOGRAPHIC_COMPARATOR));

    IntVector duplicate = IntVector.of(3, 3);
    assertEquals(0, duplicate.lowerBound(1, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(0, duplicate.upperBound(1, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(0, duplicate.lowerBound(3, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(2, duplicate.upperBound(3, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(2, duplicate.lowerBound(5, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(2, duplicate.upperBound(5, LEXICOGRAPHIC_COMPARATOR));

    IntVector unique = IntVector.of(3, 4);
    assertEquals(0, unique.lowerBound(1, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(0, unique.upperBound(1, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(0, unique.lowerBound(3, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(1, unique.upperBound(3, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(1, unique.lowerBound(4, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(2, unique.upperBound(4, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(2, unique.lowerBound(5, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(2, unique.upperBound(5, LEXICOGRAPHIC_COMPARATOR));

    // Even numbers in lexicographic order with a run of equal numbers.
    IntVector evens = IntVector.of(10, 2, 4, 4, 4, 6, 8);
    assertEquals(0, evens.lowerBound(1, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(0, evens.upperBound(1, LEXICOGRAPHIC_COMPARATOR));

    assertEquals(0, evens.lowerBound(10, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(1, evens.upperBound(10, LEXICOGRAPHIC_COMPARATOR));

    assertEquals(1, evens.lowerBound(11, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(1, evens.upperBound(11, LEXICOGRAPHIC_COMPARATOR));

    assertEquals(1, evens.lowerBound(2, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(2, evens.upperBound(2, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(2, evens.lowerBound(3, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(2, evens.upperBound(3, LEXICOGRAPHIC_COMPARATOR));

    assertEquals(2, evens.lowerBound(4, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(5, evens.upperBound(4, LEXICOGRAPHIC_COMPARATOR));

    assertEquals(5, evens.upperBound(5, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(5, evens.lowerBound(5, LEXICOGRAPHIC_COMPARATOR));

    assertEquals(6, evens.lowerBound(8, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(7, evens.upperBound(8, LEXICOGRAPHIC_COMPARATOR));

    assertEquals(7, evens.upperBound(9, LEXICOGRAPHIC_COMPARATOR));
    assertEquals(7, evens.lowerBound(9, LEXICOGRAPHIC_COMPARATOR));
  }

  /**
   * The upper and lower bound methods are one-liners with the logic tested elsewhere, so this test
   * is simple.
   */
  @Test
  public void testUpperAndLowerBoundOfSubrange() {
    // A test vector that is sorted in the subrange [1 .. 13)
    //                  index       0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11, 12, 13, 14
    IntVector vector = IntVector.of(9, 0, 1, 1, 2, 3, 4, 5, 5, 5, 7, 9, 11, 13, 0);

    // Upper and lower bounds for a unique number.
    assertEquals(5, vector.lowerBound(3, 1, 13));
    assertEquals(6, vector.upperBound(3, 1, 13));

    // Upper and lower bounds for a missing number in the range.
    assertEquals(10, vector.lowerBound(6, 1, 13));
    assertEquals(10, vector.upperBound(6, 1, 13));

    // Upper and lower bounds for a missing number outside the range.
    assertEquals(13, vector.lowerBound(20, 1, 13));
    assertEquals(13, vector.upperBound(20, 1, 13));

    // Upper and lower bounds for a repeated number.
    assertEquals(2, vector.lowerBound(1, 1, 13));
    assertEquals(4, vector.upperBound(1, 1, 13));
  }

  @SuppressWarnings("JUnitIncompatibleType")
  @Test
  public void testRotate() {
    // Rotation of an empty vector does nothing.
    IntVector empty = new IntVector();
    empty.rotate(1);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(), empty);

    // Rotation of a single-element vector does nothing.
    IntVector single = IntVector.of(1);
    single.rotate(1);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1), single);

    // Rotation of a two-element vector by one place, forward and backward.
    IntVector v12 = IntVector.of(1, 2);
    v12.rotate(1);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(2, 1), v12);
    v12.rotate(-1);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2), v12);

    // Rotation of a three-element vector.
    IntVector v123 = IntVector.of(1, 2, 3);
    v123.rotate(1); // Rotate right one place.
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(3, 1, 2), v123);
    v123.rotate(-2); // Rotate left two places.
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(2, 3, 1), v123);
  }

  @SuppressWarnings("JUnitIncompatibleType")
  @Test
  public void testInsertValue() {
    // Insert into an empty vector.
    IntVector empty = new IntVector();
    empty.insert(0, 3);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(3), empty);

    // Insert at the beginning of a vector.
    IntVector v123 = IntVector.of(1, 2, 3);
    v123.insert(0, 0);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(0, 1, 2, 3), v123);

    // Insert at the end of a vector.
    v123 = IntVector.of(1, 2, 3);
    v123.insert(3, 4);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2, 3, 4), v123);

    // Insert into the middle of a vector.
    IntVector v123456 = IntVector.of(1, 2, 3, 4, 5, 6);
    v123456.insert(2, 100);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2, 100, 3, 4, 5, 6), v123456);

    // Insert which requires increasing capacity.
    IntVector ofSize32 = IntVector.ofSize(32);
    assertEquals(32, ofSize32.size());
    assertEquals(32, ofSize32.capacity());
    ofSize32.insert(5, 5);
    assertEquals(33, ofSize32.size());
    assertEquals(64, ofSize32.capacity());
  }

  @SuppressWarnings("JUnitIncompatibleType")
  @Test
  public void testInsertSequence() {
    IntSequence toInsert = ImmutableIntSequence.of(100, 101, 102);

    // Insert into an empty vector.
    IntVector empty = new IntVector();
    empty.insert(0, toInsert);
    assertEquals(toInsert, empty);

    // Inserting an empty sequence anywhere does nothing.
    IntVector v12345 = IntVector.of(1, 2, 3, 4, 5);
    v12345.insert(0, ImmutableIntSequence.of());
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2, 3, 4, 5), v12345);
    v12345.insert(3, ImmutableIntSequence.of());
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2, 3, 4, 5), v12345);
    v12345.insert(5, ImmutableIntSequence.of());
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2, 3, 4, 5), v12345);

    // Insert at the beginning of a vector.
    IntVector v123 = IntVector.of(1, 2, 3);
    v123.insert(0, toInsert);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(100, 101, 102, 1, 2, 3), v123);

    // Insert at the end of a vector so there is nothing to shift right.
    IntVector v124 = IntVector.of(1, 2, 4);
    v124.insert(3, toInsert);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2, 4, 100, 101, 102), v124);

    // Test a case where the shifted-right end of the vector overlaps its original position.
    IntVector v123456 = IntVector.of(1, 2, 3, 4, 5, 6);
    v123456.insert(2, ImmutableIntSequence.of(100));
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2, 100, 3, 4, 5, 6), v123456);

    // Test a case where the shifted-right end of the vector does not overlap its original position.
    v123456 = IntVector.of(1, 2, 3, 4, 5, 6);
    v123456.insert(3, toInsert);
    // The types of this assertion are mismatched: type `IntVector` is not compatible with
    // `ImmutableIntSequence`. Please investigate and fix!
    assertEquals(ImmutableIntSequence.of(1, 2, 3, 100, 101, 102, 4, 5, 6), v123456);
  }

  @Test
  public void testForEach() {
    IntVector vector = IntVector.empty();
    vector.resize(10);
    vector.fillConsecutive();

    ArrayList<Integer> visited = new ArrayList<>();
    vector.forEach((IntConsumer) visited::add);
    for (int i = 0; i < 10; i++) {
      assertEquals(i, (int) visited.get(i));
    }
  }

  @Test
  public void testConstructAndAppendLists() {
    ImmutableList<Integer> fibbonaci = ImmutableList.of(1, 1, 2, 3, 5, 8, 13, 21, 34);
    // Start with a copy of the list.
    IntVector vector = IntVector.copyOf(fibbonaci);
    // Then push and append.
    vector.push(-1);
    vector.addAll(fibbonaci.reverse());

    // Validate the results.
    assertEquals(19, vector.size());

    assertEquals(1, vector.get(0));
    assertEquals(1, vector.get(1));
    assertEquals(2, vector.get(2));
    assertEquals(34, vector.get(8));
    assertEquals(-1, vector.get(9));
    assertEquals(34, vector.get(10));
    assertEquals(21, vector.get(11));
    assertEquals(1, vector.get(18));
  }

}
