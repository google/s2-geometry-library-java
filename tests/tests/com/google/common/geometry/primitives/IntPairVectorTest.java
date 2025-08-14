/*
 * Copyright 2023 Google Inc.
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

import com.google.common.geometry.primitives.Ints.IntPairComparator;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for IntPairVector. */
@RunWith(JUnit4.class)
public final class IntPairVectorTest {

  private void assertPairEquals(IntPairVector v, int index, int first, int second) {
    assertEquals(first, v.getFirst(index));
    assertEquals(second, v.getSecond(index));
  }

  @Test
  public void testBasics() {
    IntPairVector vector = new IntPairVector();
    vector.add(1, 10);
    vector.add(2, 20);
    vector.setFirst(0, 11);
    vector.setSecond(1, 12);

    assertPairEquals(vector, 0, 11, 10);
    assertPairEquals(vector, 1, 2, 12);

    vector.swap(0, 1);
    assertPairEquals(vector, 0, 2, 12);
    assertPairEquals(vector, 1, 11, 10);

    IntPairVector another = new IntPairVector();
    another.add(1, 10);
    another.add(2, 20);
    assertPairEquals(another, 0, 1, 10);
    assertPairEquals(another, 1, 2, 20);
    another.copy(vector);
    assertPairEquals(another, 0, 2, 12);
    assertPairEquals(another, 1, 11, 10);
  }

  @Test
  public void testSortAndContains() {
    IntPairVector vector = new IntPairVector();
    vector.add(3, 30);
    vector.add(1, 13);
    vector.add(1, 12);
    vector.add(1, 11);
    vector.add(2, 20);
    vector.add(2, 2);

    vector.sort();

    // After sorting, the pairs should be:
    // (1,11), (1,12), (1,13), (2,2), (2,20), (3,30)
    assertPairEquals(vector, 0, 1, 11);
    assertPairEquals(vector, 1, 1, 12);
    assertPairEquals(vector, 2, 1, 13);
    assertPairEquals(vector, 3, 2, 2);
    assertPairEquals(vector, 4, 2, 20);
    assertPairEquals(vector, 5, 3, 30);

    // The contains method requires the vector to be sorted.
    assertTrue(vector.contains(1, 11));
    assertTrue(vector.contains(3, 30));
    assertTrue(vector.contains(2, 2));
    assertFalse(vector.contains(2, 3));
    assertFalse(vector.contains(1, 10));
    assertFalse(vector.contains(0, 0));
  }

  @Test
  public void testSortWithComparator() {
    // A comparator of int pairs that compares by the second element, breaking ties with the first.
    IntPairComparator bySecond = new IntPairComparator() {
      @Override
      public int comparePair(int a0, int a1, int b0, int b1) {
        int c = Integer.compare(a1, b1);
        if (c != 0) {
          return c;
        }
        return Integer.compare(a0, b0);
      }
    };

    IntPairVector vector = new IntPairVector();
    vector.add(2, 5);
    vector.add(2, 6);
    vector.add(2, 7);
    vector.add(4, 5);
    vector.add(4, 6);
    vector.add(4, 7);

    vector.sort(bySecond);

    // After sorting, the pairs should be:
    // (2,5), (4,5), (2,6), (4,6), (2,7), (4,7)
    assertPairEquals(vector, 0, 2, 5);
    assertPairEquals(vector, 1, 4, 5);
    assertPairEquals(vector, 2, 2, 6);
    assertPairEquals(vector, 3, 4, 6);
    assertPairEquals(vector, 4, 2, 7);
    assertPairEquals(vector, 5, 4, 7);
  }

  @Test
  public void testLowerBound() {
    IntPairVector vector = new IntPairVector();
    vector.add(2, 5);
    vector.add(2, 6);
    vector.add(2, 7);
    vector.add(4, 5);
    vector.add(4, 6);
    vector.add(4, 7);

    // Case where the target first element is lower than any vector entry.
    assertEquals(0, vector.lowerBound(1, 1));
    // Case where the target second element is lower than any vector entry
    assertEquals(0, vector.lowerBound(2, 1));
    // Find the first three elements.
    assertEquals(0, vector.lowerBound(2, 5));
    assertEquals(1, vector.lowerBound(2, 6));
    assertEquals(2, vector.lowerBound(2, 7));
    // Case where the target first element is not in the vector, between existing entries
    assertEquals(3, vector.lowerBound(3, 1));
    // Case where the target second element is lower than any in the vector.
    assertEquals(3, vector.lowerBound(4, 1));
    // Find the last three elements.
    assertEquals(3, vector.lowerBound(4, 5));
    assertEquals(4, vector.lowerBound(4, 6));
    assertEquals(5, vector.lowerBound(4, 7));
    // Case where the second target element is greater than the last vector element.
    assertEquals(6, vector.lowerBound(4, 8));
    // Case where the target first element is greater than any in the vector.
    assertEquals(6, vector.lowerBound(5, 1));

    // Insert new elements at the beginning, middle, and end of the vector with unique first values.
    vector.add(1, 3);
    vector.add(3, 3);
    vector.add(5, 3);
    vector.sort();

    // Vector is now:
    //  0      1      2      3      4      5      6      7      8    size == 9
    // (1,3), (2,5), (2,6), (2,7), (3,3), (4,5), (4,6), (4,7), (5,3)

    // Find the new elements
    assertEquals(0, vector.lowerBound(1, 3));
    assertEquals(4, vector.lowerBound(3, 3));
    assertEquals(8, vector.lowerBound(5, 3));

    // Look for elements just below the existing new ones, with lower second target values
    assertEquals(0, vector.lowerBound(1, 2));
    assertEquals(4, vector.lowerBound(3, 2));
    assertEquals(8, vector.lowerBound(5, 2));

    // Look for elements just above the existing new ones, with higher second target values
    assertEquals(1, vector.lowerBound(1, 4));
    assertEquals(5, vector.lowerBound(3, 4));
    assertEquals(9, vector.lowerBound(5, 4));
  }
}
