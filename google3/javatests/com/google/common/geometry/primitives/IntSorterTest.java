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

import static com.google.common.geometry.primitives.Ints.INTEGER_COMPARATOR;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import com.google.common.geometry.primitives.Ints.IntComparator;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.Random;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link IntSorter}. */
@RunWith(JUnit4.class)
public final class IntSorterTest {

  private static final int SEED = 123455;
  private final Random rand = new Random(SEED);

  public void assertSorted(IntVector vector, IntComparator comparator) {
    OfInt iter = vector.intIterator();
    int prev = iter.nextInt();
    while (iter.hasNext()) {
      int next = iter.nextInt();
      assertTrue(comparator.compare(prev, next) <= 0);
      prev = next;
    }
  }

  /** Compare numbers lexicographically. "100" < "2". */
  private static final IntComparator LEXICOGRAPHIC_COMPARATOR =
      new IntComparator() {
        @Override
        public int compare(int a, int b) {
          return Integer.toString(a).compareTo(Integer.toString(b));
        }
      };

  /** Basic test check of IntSorter and two comparators. */
  @Test
  public void testComparators() {
    IntVector vector = IntVector.copyOf(100, 110, 11, 12, 120, 4, 32);
    assertEquals(7, vector.size());

    // Sort numerically, verify that the order is correct, and that assertSorted agrees.
    IntSorter.sort(vector);
    assertTrue(IntVector.copyOf(4, 11, 12, 32, 100, 110, 120).isEqualTo(vector));
    assertSorted(vector, INTEGER_COMPARATOR);

    // Sort lexicographically, verify the order is correct, and that assertSorted agrees.
    IntSorter.sort(LEXICOGRAPHIC_COMPARATOR, vector);
    assertTrue(IntVector.copyOf(100, 11, 110, 12, 120, 32, 4).isEqualTo(vector));
    assertSorted(vector, LEXICOGRAPHIC_COMPARATOR);
  }

  /** Test the sort() method that takes a range. */
  @Test
  public void testSortRange() {
    IntVector vector = IntVector.copyOf(new int[] {100, 110, 11, 12, 120, 4, 32});
    IntSorter.sort(vector, 1, 5);
    assertTrue(IntVector.copyOf(100, 4, 11, 12, 110, 120, 32).isEqualTo(vector));
  }

  /** Test sorting 1000 random-length sequences of random integers with different comparators. */
  @Test
  public void testNumericSortRandomSequences() {
    for (int test = 0; test < 1000; test++) {
      // Create a random-content, random-length vector.
      int size = rand.nextInt(8192);
      IntVector vector = IntVector.ofSize(size);
      for (int i = 0; i < size; i++) {
        vector.set(i, rand.nextInt());
      }
      // Sort both lexicographically and numerically. Half the time, sort lexicographically first.
      if (rand.nextBoolean()) {
        IntSorter.sort(INTEGER_COMPARATOR, vector, 0, size - 1);
        assertSorted(vector, INTEGER_COMPARATOR);
        IntSorter.sort(LEXICOGRAPHIC_COMPARATOR, vector, 0, size - 1);
        assertSorted(vector, LEXICOGRAPHIC_COMPARATOR);
      } else {
        IntSorter.sort(LEXICOGRAPHIC_COMPARATOR, vector, 0, size - 1);
        assertSorted(vector, LEXICOGRAPHIC_COMPARATOR);
        IntSorter.sort(INTEGER_COMPARATOR, vector, 0, size - 1);
        assertSorted(vector, INTEGER_COMPARATOR);
      }
    }
  }
}
