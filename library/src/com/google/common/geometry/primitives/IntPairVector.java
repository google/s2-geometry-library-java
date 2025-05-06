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

import com.google.common.base.Preconditions;
import com.google.common.geometry.primitives.Ints.IntPairComparator;
import com.google.common.geometry.primitives.Sorter.SortableCollection;

/**
 * An IntPairVector is an automatically-growing list of pairs of integers. It may be used much like
 * a {@code List<Pair<Integer, Integer>} or {@code Stack<Pair<Integer, Integer>>} but it is more
 * efficient as it avoids boxing Integers and Pairs. It also provides several convenience methods
 * making it useful as a replacement for a C++ {@code vector<pair<int, int>>} when porting to Java.
 */
public class IntPairVector {
  private final IntVector firsts;
  private final IntVector seconds;

  /** Constructs an empty IntPairVector. */
  public IntPairVector() {
    firsts = IntVector.empty();
    seconds = IntVector.empty();
  }

  /** Copy constructor. */
  public IntPairVector(IntPairVector other) {
    firsts = new IntVector(other.firsts);
    seconds = new IntVector(other.seconds);
  }

  /** Returns a new IntVector with the given ints as pairs. There must be an even number of them. */
  public static IntPairVector of(int... contents) {
    Preconditions.checkArgument(contents.length % 2 == 0);
    IntPairVector result = new IntPairVector();
    for (int i = 0; i < contents.length; i += 2) {
      result.add(contents[i], contents[i + 1]);
    }
    return result;
  }

  // This toString implementation is intended to support debugging.
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("IntPairVector of ").append(firsts.size()).append(" pairs: ");
    for (int i = 0; i < firsts.size() && i < 10; i++) {
      sb.append("(").append(firsts.get(i)).append(",").append(seconds.get(i)).append("), ");
    }
    if (firsts.size() > 10) {
      sb.append("...");
    }
    return sb.toString();
  }

  /** Adds a pair of integers to this IntPairVector. */
  public void add(int first, int second) {
    firsts.add(first);
    seconds.add(second);
  }

  /** Gets the first element of the pair at position 'index' in this IntPairVector. */
  public int getFirst(int index) {
    return firsts.get(index);
  }

  /** Gets the second element of the pair at position 'index' in this IntPairVector. */
  public int getSecond(int index) {
    return seconds.get(index);
  }

  /** Sets the pair at position 'index' to ('valFirst', 'valSecond'). */
  public void setPair(int index, int valFirst, int valSecond) {
    firsts.set(index, valFirst);
    seconds.set(index, valSecond);
  }

  /** Sets the first element of the pair at position 'index' to 'value'. */
  public void setFirst(int index, int value) {
    firsts.set(index, value);
  }

  /** Gets the second element of the pair at position 'index' to 'value'. */
  public void setSecond(int index, int value) {
    seconds.set(index, value);
  }

  /** Returns the number of pairs in this IntPairVector. */
  public int size() {
    return firsts.size();
  }

  /** Calls the given consumer with each pair in this IntPairVector. */
  public void forEach(Ints.IntBiConsumer consumer) {
    for (int i = 0; i < size(); i++) {
      consumer.accept(firsts.get(i), seconds.get(i));
    }
  }

  // TODO(torrey): A full-featured and public implementation of IntPairVector should have lowerBound
  // and upperBound implementations that take an IntPairComparator.

  /**
   * Performs binary searches in this vector and returns the lowest index containing a pair (A, B)
   * greater than or equal to the given pair ('first', 'second'), using lexicographic ordering of
   * pair elements (first element before second), and natural ordering of integers. If all pairs are
   * less than ('first', 'second'), returns 'size()'.
   *
   * <P>NOTE: This vector must be pre-sorted. Use {@link #sort()} if required.
   */
  public int lowerBound(int first, int second) {
    // The lowest index where the first existing element is greater than or equal to target.first.
    int low = firsts.lowerBound(first);
    if (low == size() || firsts.get(low) > first) {
      return low;
    }
    Preconditions.checkState(firsts.get(low) == first);
    int high = firsts.upperBound(first);
    return seconds.lowerBound(second, low, high);
  }

  /**
   * Returns true if this IntPairVector contains the given pair.
   *
   * <p>Requires that this IntPairVector is sorted by first elements, breaking ties with second
   * elements. Call {@link #sort()} to ensure this, if needed.
   */
  public boolean contains(int first, int second) {
    // Lowest index of an element >= 'first'.
    int firstLower = firsts.lowerBound(first);
    if (firstLower == firsts.size() || firsts.get(firstLower) != first) {
      return false; // 'firsts' vector doesn't contain 'first'.
    }

    // If the number of pairs beginning with 'first' is small enough, do a linear search instead of
    // a binary search.
    if (firstLower + 8 >= firsts.size() || firsts.get(firstLower + 8) != first) {
      for (int i = firstLower; i < firsts.size() && first == firsts.get(i); i++) {
        if (second == seconds.get(i)) {
          return true;
        }
      }
      return false;
    }

    // Lowest index of an element > 'first'.
    int firstUpper = firsts.upperBound(first);

    // All the pairs from firstLower (inclusive) to firstUpper (exclusive) have 'first' as their
    // first element. Search that subrange of 'seconds' for 'second'.
    int secondLower = seconds.lowerBound(second, firstLower, firstUpper);
    return secondLower < seconds.size() && seconds.get(secondLower) == second;
  }

  /**
   * Exchanges the content of this IntPairVector with the content of the 'other' IntPairVector. The
   * size of the IntPairVector may differ. This is a constant-time operation, regardless of the size
   * of the IntPairVectors.
   */
  public void swap(IntPairVector other) {
    firsts.swap(other.firsts);
    seconds.swap(other.seconds);
  }

  /** Swaps the pairs at index 'a' and index 'b', without bounds checks. */
  public void swap(int indexA, int indexB) {
    firsts.swap(indexA, indexB);
    seconds.swap(indexA, indexB);
  }

  /** Sorts the pairs using the provided IntPairComparator. */
  public void sort(IntPairComparator cmp) {
    new SortableCollection() {
      @Override
      public int size() {
        return firsts.size();
      }
      @Override
      public void truncate(int end) {
        firsts.truncate(end);
        seconds.truncate(end);
      }
      @Override
      public void swap(int a, int b) {
        firsts.swap(a, b);
        seconds.swap(a, b);
      }
      @Override
      public boolean less(int a, int b) {
        return cmp.comparePair(firsts.get(a), seconds.get(a), firsts.get(b), seconds.get(b)) < 0;
      }
    }.sort();
  }

  /**
   * Sorts the pairs by the numeric value of the first element of the pair, breaking ties with the
   * value of the second element of the pair.
   */
  public void sort() {
    sort((firstA, secondA, firstB, secondB) -> Integer.compare(firstA, firstB) != 0
        ? Integer.compare(firstA, firstB)
        : Integer.compare(secondA, secondB));
  }

  /**
   * Ensure that at least 'newCapacity' pairs can be stored in the IntPairVector, increasing or
   * possibly decreasing the sizes of the underlying arrays as required. Does not change the number
   * or content of pairs already stored.
   */
  public void ensureCapacity(int newCapacity) {
    firsts.ensureCapacity(newCapacity);
    seconds.ensureCapacity(newCapacity);
  }

  /**
   * Replace the current contents of this IntPairVector with a copy of the given IntPairVector
   * 'other'.
   */
  public void copy(IntPairVector other) {
    firsts.copy(other.firsts);
    seconds.copy(other.seconds);
  }

  /** Clear the contents of this IntPairVector. */
  public void clear() {
    firsts.clear();
    seconds.clear();
  }

  /** Returns true if this IntPairStack has no elements. */
  public boolean isEmpty() {
    return firsts.isEmpty();
  }

  /**
   * Returns true if this IntPairVector currently contains the same pairs of ints as 'other', in the
   * same order. Does not override Object.equals() because IntPairVector is mutable.
   */
  public boolean isEqualTo(IntPairVector other) {
    return firsts.isEqualTo(other.firsts) && seconds.isEqualTo(other.seconds);
  }

  // Stack methods.

  /** Adds a new pair to the top of the stack, which is the end of the vector. */
  public void pushPair(int first, int second) {
    firsts.push(first);
    seconds.push(second);
  }

  /** Returns the first element at the top of the stack (the end of the vector) of pairs. */
  public int peekFirst() {
    return firsts.peek();
  }

  /** Returns the second element at the top of the stack (the end of the vector) of pairs. */
  public int peekSecond() {
    return seconds.peek();
  }

  /** Removes the top pair from the stack, which is the pair at the end of the vector. */
  public void pop() {
    firsts.pop();
    seconds.pop();
  }
}
