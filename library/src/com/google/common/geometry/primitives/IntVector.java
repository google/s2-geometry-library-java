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

import static java.lang.Math.max;

import com.google.common.geometry.primitives.Ints.IntBiConsumer;
import com.google.common.geometry.primitives.Ints.IntComparator;
import com.google.common.geometry.primitives.Ints.IntConsumer;
import com.google.common.geometry.primitives.Ints.IntList;
import com.google.common.geometry.primitives.Ints.IntPredicate;
import com.google.common.geometry.primitives.Ints.IntSequence;
import com.google.common.geometry.primitives.Ints.MutableIntList;
import com.google.common.geometry.primitives.Ints.OfInt;
import com.google.common.primitives.Ints;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.AbstractList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.Random;

/**
 * IntVector is a general purpose implementation of {@link MutableIntList}, and analogous to {@code
 * List<Integer>} for primitive ints, except that nulls are not permitted. It is intended to be a
 * replacement for {@code std::vector<int>} when porting C++ to Java, and high performance is a
 * priority.
 *
 * <p>This implementation stores values in an array which is doubled in size when it runs out of
 * space, so an IntVector may use 2x as much space as an array of exactly the required size.
 * Doubling size while growing is somewhat expensive, but may often be avoided or minimized by
 * constructing an IntVector of the appropriate size, or using {@link #ensureCapacity(int)} before
 * adding many elements.
 *
 * <p>Benchmarks in {@link IntVectorBenchmark} compare IntVector to raw arrays and to an {@code
 * ArrayList<Integer>}. IntVector is nearly as fast as a raw int[] array for many access patterns,
 * including random updates, but for simple linear loops, arrays may be twice as fast as an
 * IntVector. However, IntVector takes less than half the time as an ArrayList to grow by repeatedly
 * adding elements. Similarly, it takes less than half the time as an ArrayList to fill and then
 * read the contents, and one sixth of the time to make a series of random updates.
 *
 * <p>For fast sorting in natural order, Arrays.sort() can be used directly on the IntVector data
 * using {@link #asArray()}. This is always the fastest option, (over twice as fast as
 * ArrayList.sort(naturalOrder())) but be careful to only sort the part of the returned array that
 * contains valid data, i.e. Arrays.sort(vector.data, 0, vector.size()). A {@link #sort()} method is
 * provided which does exactly that, for convenience.
 *
 * <p>For sorting with a Comparator, IntVector provides a {@link #sort(IntComparator)} method
 * implemented using the fastest currently available stable sort. This is the recommended option
 * for most clients. However, the best sort performance might depend on situation details. IntVector
 * provides asList(), which provides an AbstractList that may be sorted with the highly optimized
 * List.sort(). List.sort() is also a stable sort, i.e. does not reorder elements that compare as
 * equal. The only downside is it boxes ints. The other option is to use {@link IntSorter}, which
 * avoids boxing but usually does more comparisons, and despite avoiding boxing currently benchmarks
 * as slower. Also, the current implementation of IntSorter is not a stable sort.
 *
 * <p>Besides the speed advantages over {@code ArrayList<Integer>}, in more complex algorithms, the
 * avoidance of boxing Integers and consequent memory and garbage collector pressure should provide
 * additional advantages.
 */
public class IntVector implements MutableIntList {
  /** The default initial capacity. */
  private static final int DEFAULT_INITIAL_CAPACITY = 16;

  /** Storage of elements in the vector. Only elements [0 .. numElements-1] are actually used. */
  private int[] data;

  /** The number of elements in the vector. */
  private int numElements = 0;

  /** Constructs an empty IntVector. */
  public IntVector() {
    data = new int[DEFAULT_INITIAL_CAPACITY];
  }

  /** Copy constructor. */
  public IntVector(IntVector other) {
    this.data = new int[other.data.length];
    System.arraycopy(other.data, 0, this.data, 0, other.numElements);
    this.numElements = other.numElements;
  }

  /**
   * Constructs an IntVector with an initial size and sufficient capacity for that size. Private to
   * avoid confusion of size and capacity: use either the {@link IntVector#ofSize()} or
   * {@link IntVector#ofCapacity()} factory methods.
   */
  private IntVector(int initialSize) {
    int powerOfTwo = Integer.highestOneBit(initialSize);
    powerOfTwo = (powerOfTwo == initialSize) ? powerOfTwo : powerOfTwo << 1;
    data = new int[max(DEFAULT_INITIAL_CAPACITY, powerOfTwo)];
    numElements = initialSize;
  }

  /** Factory method returning an empty IntVector. */
  public static IntVector empty() {
    return new IntVector();
  }

  /** Factory method returning an IntVector initialized with the given size. */
  public static IntVector ofSize(int size) {
    return new IntVector(size);
  }

  /** Factory method returning an empty IntVector with at least the given capacity. */
  public static IntVector ofCapacity(int capacity) {
    IntVector vec = new IntVector();
    vec.ensureCapacity(capacity);
    return vec;
  }

  // TODO(torrey): Add an IntVector.of(int[] ints) factory method which just assigns the given ints
  // to data, without copying. Add unit test coverage to ensure this works as expected with arrays
  // that are not sized as a power of two. It would be nice if it was possible to provide an
  // asArray() or release() method that truncated the internal data array to length numElements and
  // returned it without copying. With those two methods, converting between arrays and IntVectors
  // would be cheap. Unfortunately Java cannot truncate an array. Also, implement shrinkToFit() that
  // reduces the data size to exactly the number of elements required. Make other changes as needed
  // to handle non-power-of-two data array sizes.

  /** Factory method returning an IntVector initialized with the given ints. */
  public static IntVector of(int... ints) {
    return IntVector.copyOf(ints);
  }

  /**
   * Factory method returning an IntVector initialized with the same size as, and containing a copy
   * of the given values.
   */
  public static IntVector copyOf(int[] values) {
    IntVector vec = new IntVector();
    vec.resize(values.length);
    System.arraycopy(values, 0, vec.data, 0, values.length);
    return vec;
  }

  /**
   * Factory method returning an IntVector initialized with the same size as, and containing a copy
   * of the given values.
   */
  public static IntVector copyOf(double... values) {
    IntVector vec = IntVector.ofCapacity(values.length);
    Arrays.stream(values).mapToInt(d -> (int) d).forEach(vec::push);
    return vec;
  }

  /**
   * Factory method returning an IntVector initialized with the same size as, and containing a copy
   * of the given IntSequence 'values'.
   */
  public static IntVector copyOf(IntSequence values) {
    IntVector vec = new IntVector();
    vec.ensureCapacity(values.size());
    values.forEach((IntConsumer) vec::push);
    return vec;
  }

  /**
   * Factory method returning an IntVector initialized with the same size as, and containing a copy
   * of the given Iterable of Integer 'values'.
   */
  public static IntVector copyOf(Iterable<Integer> values) {
    IntVector vec = new IntVector();
    for (Integer value : values) {
      vec.push(value);
    }
    return vec;
  }

  /** Replace the current contents of this IntVector with a copy of the given IntSequence. */
  public void copy(IntSequence sequence) {
    clear();
    addAll(sequence);
  }

  @Override
  public String toString() {
    return debugString();
  }

  /** Sets the number of elements to zero. Keeps the existing data storage. */
  @Override
  public void clear() {
    numElements = 0;
  }

  /** Fills all elements of this IntVector with the given value. */
  @Override
  public void fill(int value) {
    Arrays.fill(data, value);
  }

  /** Reverses the ordering of all the elements in the IntVector. */
  public void reverse() {
    Ints.reverse(data, 0, numElements);
  }

  /**
   * Reverses the ordering of elements in the IntVector, from 'fromIndex' (inclusive) to 'toIndex'
   * (exclusive).
   */
  public void reverse(int fromIndex, int toIndex) {
    validateFromAndToIndex(fromIndex, toIndex);
    Ints.reverse(data, fromIndex, toIndex);
  }

  /**
   * Fills this IntVector with consecutive integers beginning with zero, keeping the current
   * number of elements but replacing all the values.
   */
  @Override
  public void fillConsecutive() {
    for (int i = 0; i < numElements; i++) {
      data[i] = i;
    }
  }

  /** The number of elements in the IntVector. */
  @Override
  public int size() {
    return numElements;
  }

  @Override
  public boolean isEmpty() {
    return (numElements == 0);
  }

  /**
   * The total number of elements that can be stored before reallocating the array will be needed.
   */
  public int capacity() {
    return data.length;
  }

  /**
   * Decrease the size of the data array to the smallest power of two equal or greater than the
   * current number of elements, or DEFAULT_INITIAL_CAPACITY, whichever is larger.
   */
  public void shrink() {
    int powerOfTwo = Integer.highestOneBit(numElements);
    powerOfTwo = (powerOfTwo == numElements) ? powerOfTwo : powerOfTwo << 1;
    int newDataLength = max(DEFAULT_INITIAL_CAPACITY, powerOfTwo);
    if (newDataLength == data.length) {
      return;
    }
    data = Arrays.copyOf(data, newDataLength);
  }

  /**
   * Ensure that at least 'newCapacity' elements can be stored in the backing data array, increasing
   * the capacity by allocating a new array and copying existing contents to it, if required. Does
   * not change the number of elements in the vector. Returns true if a new array was allocated,
   * false if the existing array had sufficient capacity.
   *
   * <p>The new capacity will be the smallest power of two equal or greater than 'newCapacity', or
   * the current capacity, or DEFAULT_INITIAL_CAPACITY, whichever is larger.
   */
  @CanIgnoreReturnValue
  public boolean ensureCapacity(int newCapacity) {
    int minSize = max(newCapacity, numElements);
    int powerOfTwo = Integer.highestOneBit(minSize);
    powerOfTwo = (powerOfTwo == newCapacity) ? powerOfTwo : powerOfTwo << 1;
    int newDataLength = max(DEFAULT_INITIAL_CAPACITY, powerOfTwo);
    if (newDataLength <= data.length) {
      return false;
    }

    // Copy the existing data into a new, larger array.
    data = Arrays.copyOf(data, newDataLength);
    return true;
  }

  private void validateFromAndToIndex(int fromIndex, int toIndex) {
    if (fromIndex < 0 || fromIndex >= numElements) {
      throw new IndexOutOfBoundsException(
          "fromIndex " + fromIndex + " out of bounds on vector of size " + numElements);
    }
    if (toIndex > numElements) {
      throw new IndexOutOfBoundsException(
          "toIndex " + toIndex + " out of bounds on vector of size " + numElements);
    }
    int size = toIndex - fromIndex;
    if (size < 0) {
      throw new IllegalArgumentException("toIndex less than fromIndex");
    }
  }

  // TODO(torrey): Consider adding a mutableSubList() method. Also, improve unit test coverage of
  // subList behavior.
  /**
   * Returns a view of the contents of this IntVector from {@code fromIndex} (inclusive) to {@code
   * toIndex} (exclusive).
   */
  @Override
  public IntList subList(int fromIndex, int toIndex) {
    validateFromAndToIndex(fromIndex, toIndex);

    return new IntList() {
      final int lo = fromIndex; // inclusive
      final int hi = toIndex; // exclusive

      @Override
      public int get(int index) {
        return data[index + lo];
      }

      @Override
      public OfInt intIterator() {
        return new OfInt() {
          int index = lo;

          @Override
          public boolean hasNext() {
            return index < hi;
          }

          @Override
          public int nextInt() {
            if (!hasNext()) {
              throw new IndexOutOfBoundsException();
            }
            return data[index++];
          }
        };
      }

      @Override
      public void forEach(IntConsumer action) {
        for (int i = lo; i < hi; i++) {
          action.accept(data[i]);
        }
      }

      @Override
      public void forEach(IntBiConsumer action) {
        for (int i = lo; i < hi; i++) {
          action.accept(i - lo, data[i]);
        }
      }

      @Override
      public int size() {
        return hi - lo;
      }
    };
  }

  /**
   * Provides an {@code AbstractList<Integer>} view of this IntVector: modifications to the
   * IntVector modify the AbstractList, and vice versa.
   */
  @Override
  public AbstractList<Integer> asList() {
    return new AbstractList<>() {
      @Override
      public Integer get(int index) {
        if (index > numElements) {
          // TODO(torrey): Here and below should be IndexOutOfBoundsException to match the Java.List
          // API, but J2CL doesn't provide that.
          throw new NoSuchElementException("Index out of bounds: " + index);
        }
        return data[index];
      }

      @Override
      public Integer set(int index, Integer value) {
        if (index > numElements) {
          throw new NoSuchElementException("Index out of bounds: " + index);
        }
        int p = data[index];
        data[index] = value;
        return p;
      }

      @Override
      public int size() {
        return numElements;
      }
    };
  }

  /**
   * Returns an int[] view of the data in this IntVector. Modifications to the returned array
   * modify the IntVector and vice versa. Note that the returned array is typically larger
   * than size(), so the user must use vector.size() rather than vector.asArray().length.
   */
  public int[] asArray() {
    return data;
  }

  /**
   * Returns a new int[] containing a copy of the contents of this IntVector and with length
   * equal to this IntVector's size().
   */
  @Override
  public int[] toArray() {
    return Arrays.copyOf(data, numElements);
  }

  /**
   * Increases the size of the vector to the requested 'newSize' if the current size is smaller. New
   * elements have initial value zero. Truncates the vector if the current size is greater. Does
   * nothing if the current size is equal to the requested newSize.
   */
  @Override
  public void resize(int newSize) {
    if (newSize < numElements) {
      truncate(newSize);
    }
    if (newSize > numElements) {
      enlarge(newSize);
    }
  }

  /**
   * Does nothing if "newSize" is larger than or equal to the current size(). Otherwise, removes all
   * elements from this IntVector, from the element at the index "newSize" to the end, decreasing
   * the size of the vector to "newSize". Does not shrink the underlying data storage.
   */
  @Override
  public void truncate(int newSize) {
    if (newSize >= numElements) {
      return;
    }
    numElements = newSize;
  }

  /**
   * Does nothing if "newSize" is less than or equal to the current size(). Otherwise, adds elements
   * with initial value zero to the end of the vector, increasing the size of the vector to
   * "newSize".
   */
  public void enlarge(int newSize) {
    if (newSize < numElements) {
      return;
    }
    // The backing array is being reused and must be cleared.
    ensureCapacity(newSize);
    Arrays.fill(data, numElements, newSize, 0);
    numElements = newSize;
  }

  @Override
  public void add(int value) {
    if (numElements == data.length) {
      ensureCapacity(numElements + 1);
    }
    data[numElements] = value;
    numElements++;
  }

  @Override
  public void push(int value) {
    add(value);
  }

  @Override
  @CanIgnoreReturnValue
  public int pop() {
    int head = data[numElements - 1];
    numElements--;
    return head;
  }

  @Override
  public int peek() {
    return data[numElements - 1];
  }

  /** Returns the value at the last index in the vector. */
  public int back() {
    return data[numElements - 1];
  }

  /** Increments the value at index and returns the result. */
  @Override
  @CanIgnoreReturnValue
  public int incrementAt(int index) {
    if (index < 0 || index >= numElements) {
      throw new ArrayIndexOutOfBoundsException(
          "incrementAt index " + index + " out of bounds on vector of size " + numElements);
    }
    data[index] = data[index] + 1;
    return data[index];
  }

  @Override
  @CanIgnoreReturnValue
  public int incrementAt(int index, int amount) {
    if (index < 0 || index >= numElements) {
      throw new ArrayIndexOutOfBoundsException(
          "incrementAt index " + index + " out of bounds on vector of size " + numElements);
    }
    data[index] += amount;
    return data[index];
  }

  /** Decrements the value at index and returns the result. */
  @Override
  @CanIgnoreReturnValue
  public int decrementAt(int index) {
    if (index < 0 || index >= numElements) {
      throw new ArrayIndexOutOfBoundsException(
          "decrementAt index " + index + " out of bounds on vector of size " + numElements);
    }
    data[index] = data[index] - 1;
    return data[index];
  }

  @Override
  @CanIgnoreReturnValue
  public int decrementAt(int index, int amount) {
    if (index < 0 || index >= numElements) {
      throw new ArrayIndexOutOfBoundsException(
          "decrementAt index " + index + " out of bounds on vector of size " + numElements);
    }
    data[index] -= amount;
    return data[index];
  }

  @Override
  public int get(int index) {
    if (index < 0 || index >= numElements) {
      throw new ArrayIndexOutOfBoundsException(
          "get at index " + index + " out of bounds on vector of size " + numElements);
    }
    return data[index];
  }

  @Override
  public void set(int index, int value) {
    if (index < 0 || index >= numElements) {
      throw new ArrayIndexOutOfBoundsException(
          "set at index " + index + " out of bounds on vector of size " + numElements);
    }
    data[index] = value;
  }

  /**
   * Sort the vector contents using the natural ordering of integers. For more sorting options, such
   * as using a custom comparator, see the class javadoc discussion.
   */
  @Override
  public void sort() {
    Arrays.sort(data, 0, size());
  }

  /**
   * Sorts the vector's contents using the provided IntComparator. This is a stable sort; that is,
   * will not change the relative ordering of elements that compare as equal according to the given
   * IntComparator.
   */
  public void sort(IntComparator comparator) {
    Collections.sort(asList(), comparator);
  }

  /** Returns true if this IntVector is sorted according to the natural ordering of integers. */
  public boolean isSorted() {
    if (numElements <= 1) {
      return true;
    }
    for (int i = 1; i < numElements; i++) {
      if (data[i - 1] > data[i]) {
        return false;
      }
    }
    return true;
  }


  /** Shuffles the order of the vector randomly, using the provided Random generator. */
  public void shuffle(Random random) {
    for (int i = 0; i < size(); i++) {
      swap(i, random.nextInt(size()));
    }
  }

  /**
   * Does a linear search through the vector to find and return the lowest index where the stored
   * int is equal to the given 'value'. If the value is not found, returns -1.
   */
  public int firstIndexOf(int value) {
    for (int i = 0; i < numElements; i++) {
      if (data[i] == value) {
        return i;
      }
    }
    return -1;
  }

  /**
   * Removes duplicate adjacent elements from this vector, copying values downward and reducing the
   * number of elements.
   *
   * <p>NOTE: unique() does not sort the vector before removing duplicate adjacent elements. Use
   * {@link #sort()} if required.
   */
  @Override
  public void unique() {
    if (numElements <= 1) {
      return;
    }
    int dst = 0;
    for (int src = 1; src < numElements; src++) {
      if (data[src] != data[dst]) {
        dst++;
        data[dst] = data[src];
      }
    }
    numElements = dst + 1;
  }

  /**
   * Performs a binary search in this vector and returns the lowest index containing an element
   * greater than or equal to 'key' using the natural ordering of integers. If all elements are less
   * than 'key', returns 'size()'.
   *
   * <P>NOTE: This vector must be pre-sorted with natural ordering. Use {@link #sort()} if required.
   */
  public int lowerBound(int key) {
    return lowerBound(
        0,
        numElements,
        // targetIsGreater returns true if the target 'key' is greater than the value at 'index'.
        index -> Integer.compare(key, data[index]) > 0);
  }

  /**
   * Performs a binary search in this vector and returns the lowest index containing an element
   * greater than or equal to 'key', as determined by the given IntComparator. If all elements are
   * less than 'key', returns 'size()'.
   *
   * <p>NOTE: This vector must be pre-sorted according to the given IntComparator. Use
   * {@link #sort(IntComparator)} if required.
   */
  public int lowerBound(int key, IntComparator comparator) {
    return lowerBound(
        0,
        numElements,
        // targetIsGreater returns true if the target 'key' is greater than the value at 'index'.
        index -> comparator.compare(key, data[index]) > 0);
  }

  /**
   * Performs a binary search in the subrange of vector defined by 'low' (inclusive) to 'high'
   * (exclusive) and returns the lowest index containing an element greater than or equal to 'key'
   * using the natural ordering of integers. If all elements in the subrange are less than 'key',
   * returns 'high'.
   *
   * <P>NOTE: The [low, high) subrange of this vector must be pre-sorted with natural ordering.
   */
  public int lowerBound(int key, int low, int high) {
    return lowerBound(low, high, index -> key > data[index]);
  }

  /**
   * Performs a binary search in the subrange of vector defined by 'low' (inclusive) to 'high'
   * (exclusive) and returns the lowest index containing an element greater than or equal to 'key',
   * as determined by the given IntComparator. If all elements are less than 'key', returns 'high'.
   *
   * <p>NOTE: The [low, high) subrange of this vector must be pre-sorted according to the given
   * IntComparator.
   */
  public int lowerBound(int key, int low, int high, IntComparator comparator) {
    return lowerBound(low, high, index -> comparator.compare(key, data[index]) > 0);
  }

  // TODO(torrey): The static lowerBound and upperBound methods below are copied directly from
  // S2ShapeUtil to avoid a circular dependency of primitives -> geometry -> primitives. These
  // should be somewhere else, perhaps Guava, or some other open-sourced utility library.

  /**
   * Returns the lowest index in the range {@code [low, high)} not smaller than a target. If every
   * value is smaller, returns {@code high}.
   */
  public static int lowerBound(int low, int high, IntPredicate targetIsGreater) {
    while (low < high) {
      int middle = low + (high - low) / 2;
      if (targetIsGreater.test(middle)) {
        low = middle + 1;
      } else {
        high = middle;
      }
    }
    return low;
  }

  /**
   * Performs a binary search in this vector and returns the lowest index containing an element
   * greater than 'key', using the natural ordering of integers. If no elements are greater than
   * 'key', returns size().
   *
   * <p>NOTE: This vector must be pre-sorted with natural ordering. Use {@link #sort()} if required.
   */
  public int upperBound(int key) {
    return upperBound(0, numElements, index -> key < data[index]);
  }

  /**
   * Performs a binary search in this vector and returns the lowest index containing an element
   * greater than 'key', as determined by the given IntComparator. If no element compares as greater
   * than 'key', returns size().
   *
   * <p>NOTE: This vector must be pre-sorted according to the given IntComparator. Use
   * {@link #sort(IntComparator)} if required.
   */
  public int upperBound(int key, IntComparator comparator) {
    return upperBound(0, numElements, index -> comparator.compare(key, data[index]) < 0);
  }

  /**
   * Performs a binary search in the subrange of vector defined by 'low' (inclusive) to 'high'
   * (exclusive)  and returns the lowest index containing an element
   * greater than 'key', using the natural ordering of integers. If no elements are greater than
   * 'key', returns size().
   *
   * <p>NOTE: The [low, high) subrange of this vector must be pre-sorted with natural ordering.
   */
  public int upperBound(int key, int low, int high) {
    return upperBound(low, high, index -> Integer.compare(key, data[index]) < 0);
  }

  /**
   * Performs a binary search in the subrange of vector defined by 'low' (inclusive) to 'high'
   * (exclusive) and returns the lowest index containing an element greater than 'key', as
   * determined by the given IntComparator. If no element compares as greater than 'key', returns
   * size().
   *
   * <p>NOTE: The [low, high) subrange of this vector must be pre-sorted according to the given
   * IntComparator.
   */
  public int upperBound(int key, int low, int high, IntComparator comparator) {
    return upperBound(
        low,
        high,
        // targetIsSmaller returns true if the target 'key' is less than the value at 'index'.
        index -> comparator.compare(key, data[index]) < 0);
  }

  /**
   * Returns the lowest index in the range {@code [low, high)} greater than a target. If every value
   * is less than or equal to the target, returns {@code high}.
   */
  public static int upperBound(int low, int high, IntPredicate targetIsSmaller) {
    while (low < high) {
      int middle = low + (high - low) / 2;
      if (targetIsSmaller.test(middle)) {
        high = middle;
      } else {
        low = middle + 1;
      }
    }
    return low;
  }
  /**
   * Performs a right rotation of this IntVector of 'distance' places, so that the first element is
   * moved to index 'distance'. A negative 'distance' will rotate left. An element at index
   * {@code i} is moved to index {@code (distance + i) mod size()}.
   */
  public void rotate(int distance) {
    if (numElements < 2 || distance == 0) {
      return;
    }
    Ints.rotate(data, distance, 0, numElements);
  }

  /**
   * Exchanges the content of this IntVector with the content of the 'other' IntVector. The size
   * of the vectors may differ. This is a constant-time operation, regardless of the size of the
   * vectors.
   */
  public void swap(IntVector other) {
    int[] tmpData = data;
    int tmpNumElements = numElements;

    data = other.data;
    numElements = other.numElements;

    other.data = tmpData;
    other.numElements = tmpNumElements;
  }

  @Override
  public boolean less(int leftIndex, int rightIndex) {
    return data[leftIndex] < data[rightIndex];
  }

  /**
   * Swaps the elements at index 'a' and index 'b', without bounds checks. Call
   * {@link IntList#rangeCheck(int, int)} to check bounds before a sequence of swap() operations.
   */
  @Override
  public void swap(int indexA, int indexB) {
    int tmp = data[indexA];
    data[indexA] = data[indexB];
    data[indexB] = tmp;
  }

  @Override
  public void forEach(IntConsumer action) {
    for (int i = 0; i < numElements; i++) {
      action.accept(data[i]);
    }
  }

  @Override
  public void forEach(IntBiConsumer action) {
    for (int i = 0; i < numElements; i++) {
      action.accept(i, data[i]);
    }
  }

  @Override
  public void addAll(IntSequence values) {
    values.forEach((IntConsumer) this::add);
  }

  @Override
  public void addAll(List<Integer> values) {
    for (Integer value : values) {
      add(value);
    }
  }

  /**
   * Inserts the given 'value' into this IntVector at position 'pos', which must be between 0 and
   * {@code size()} inclusive. The vector elements currently at 'pos' to the end will be shifted
   * right by 1.
   */
  public void insert(int pos, int value) {
    int newSize = numElements + 1;
    ensureCapacity(newSize);
    System.arraycopy(data, pos, data, pos + 1, numElements - pos);
    data[pos] = value;
    numElements = newSize;
  }

  /**
   * Inserts the given 'values' into this IntVector beginning at position 'start', which must be
   * between 0 and {@code size()} inclusive. The vector elements currently at 'start' to the end
   * will be shifted right by {@code values.size()}.
   */
  public void insert(int start, IntSequence values) {
    // Enlarge the underlying array if needed.
    int newSize = numElements + values.size();
    ensureCapacity(newSize);

    // Move the current contents of 'data' from 'start' onward upward by values.size().
    System.arraycopy(data, start, data, start + values.size(), numElements - start);

    // Add the 'values' beginning at 'start', overwriting current data.
    numElements = start;
    values.forEach(this::add);
    numElements = newSize;
  }

  @Override
  public OfInt intIterator() {
    return new OfInt() {
      int position = 0;

      @Override
      public boolean hasNext() {
        return position < numElements;
      }

      @Override
      public int nextInt() {
        return data[position++];
      }
    };
  }
}
