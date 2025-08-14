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

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.util.Random;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Tests for {@link Sorter}. */
@RunWith(JUnit4.class)
public final class SorterTest {

  private static final int SEED = 123455;
  private final Random rand = new Random(SEED);

  // A simple collection of pairs of integers for testing sorting.
  private static class SortablePairs implements Sorter.SortableCollection {
    final IntVector firsts;
    final IntVector seconds;

    public SortablePairs() {
      firsts = new IntVector();
      seconds = new IntVector();
    }

    public void add(int first, int second) {
      firsts.add(first);
      seconds.add(second);
    }

    @Override
    public int size() {
      return firsts.size();
    }

    @Override
    public void truncate(int start) {
      firsts.truncate(start);
      seconds.truncate(start);
    }

    @Override
    public void swap(int a, int b) {
      firsts.swap(a, b);
      seconds.swap(a, b);
    }

    // Order by first, breaking ties with second.
    @Override
    public boolean less(int a, int b) {
      if (firsts.get(a) == firsts.get(b)) {
        return seconds.get(a) < seconds.get(b);
      }
      return firsts.get(a) < firsts.get(b);
    }
  }

  // Checks that the SortablePairs are ordered according to the less() method it provides.
  public void assertSorted(SortablePairs collection) {
    for (int i = 1; i < collection.size(); i++) {
      // There may be ties, so rather than asserting that the pair at position (i-1) is less than
      // the pair at position (i), we assert that the pair at position (i) is not less than the pair
      // at position (i-1).
      assertFalse(collection.less(i, i - 1));
    }
  }

  // Checks that the SortablePairs are unique if ordered by the less() method it provides.
  public void assertUnique(SortablePairs collection) {
    for (int i = 1; i < collection.size(); i++) {
      assertTrue(collection.less(i - 1, i));
    }
  }

  /** Test sorting 1000 random-length sequences of random pairs. */
  @Test
  public void testSortRandomPairSequences() {
    for (int test = 0; test < 1000; test++) {
      // Create a random-content, random-length collection of pairs where each first value is in the
      // range [0, 10] and each second value is in the range [0, 100].
      int size = rand.nextInt(8192);
      SortablePairs pairCollection = new SortablePairs();
      for (int i = 0; i < size; i++) {
        pairCollection.add(rand.nextInt(10), rand.nextInt(100));
      }

      // Sort and verify that it is sorted.
      pairCollection.sort();
      assertSorted(pairCollection);
      pairCollection.unique();
      assertUnique(pairCollection);
    }
  }
}
