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

import com.google.common.geometry.primitives.Sorter.SortableCollection;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.AbstractList;
import java.util.Comparator;
import java.util.List;

/**
 * A pull model for accessing encoded data. Java's Collection classes are inefficient for large
 * numbers of small objects due to boxing and heap fragmentation. This API avoids those problems
 * by repeatedly decoding from a memory-efficient layout into a mutable instance of a more usable
 * object model. For example, instead of storing a collection of objects of 3 ints each, (12 bytes
 * of real data + 24 bytes of overhead) we can store 3 parallel arrays to have almost no overhead.
 * Then the parallel arrays can be decoded into a mutable object of 3 ints for use.
 *
 * The root interface, "Pullable" supports iteration, similar to Java Iterable. PullCollection and
 * PullList, also defined here, are comparable to Java's Collection and List, respectively.
 */
public interface Pullable<T> {
  /**
   * PullIterator provides iteration over the entries of a Pullable. Implementations take an
   * instance of the Mutable value type to fill in. Each call to pull() copies the next value into
   * that instance.
   */
  public interface PullIterator {
    /**
     * Returns false if there is no next value available, otherwise updates the mutable T with a
     * copy of the next element in the PullIterable and returns true.
     */
    boolean pull();
  }

  /** Returns a PullIterator over the collection that will repeatedly fill the given value. */
  PullIterator iterator(T value);

  /** Repeatedly pulls a new value into the given "result" and then runs the given action.  */
  default void forEach(T result, Runnable action) {
    for (PullIterator it = iterator(result); it.pull(); ) {
      action.run();
    }
  }

  /**
   * PullCollection extends {@link Pullable}. It is similar to a java.util.Collection, but owns and
   * allocates space for its own elements, rather than containing references to objects allocated
   * elsewhere. Element values are copied in and out of a PullCollection.
   */
  public interface PullCollection<T> extends Pullable<T> {
    /**
     * Returns an instance of T which is not a member of the collection. Implementations may keep a
     * pool of these objects for reuse, and the returned value is not necessarily in a default
     * state. Ideally, the returned value should be passed to {@link #destroyElement(T)} when it is
     * no longer needed. However, a typical implementation may simply return {@code new T();}.
     */
    T newElement();

    /** Accepts a T previously returned by newElement() for reuse or cleanup. */
    default void destroyElement(T value) {}

    /**
     * Ensure that at least 'capacity' entries can be stored in the backing storage without
     * additional later allocation. Returns true if storage was increased, false if the existing
     * storage had sufficient capacity. Does not change the number of elements in the collection.
     *
     * <p>May throw UnsupportedOperationException if the implementation does not support increasing
     * capacity.
     */
    @CanIgnoreReturnValue
    boolean ensureCapacity(int capacity);

    /** The number of elements in the collection. */
    int size();

    /** Returns true if the collection contains no elements. */
    default boolean isEmpty() {
      return size() == 0;
    }

    /** Removes all elements from the collection. */
    void clear();

    /** Adds a copy of the given value to the collection. */
    void add(T value);

    /** Adds copies of all the given values to the collection. */
    default void addAll(Pullable<T> values) {
      T val = newElement();
      for (PullIterator it = values.iterator(val); it.pull(); ) {
        add(val);
      }
      destroyElement(val);
    }
  }

  /**
   * A PullList extends PullCollection to provide random access to elements by index, like a
   * java.util.List.
   */
  public interface PullList<T> extends PullCollection<T> {
    /**
     * Copies the value of the element at the given index into the given value. Throws
     * IndexOutOfBoundsException if the index is out of range.
     */
    void get(int index, T value);

    /**
     * Sets the value at the given index to a copy of the given value. Throws
     * IndexOutOfBoundsException if the index is out of range.
     */
    void set(int index, T value);

    /** Returns an iterator over the list in order by index. */
    @Override
    default PullIterator iterator(T value) {
      return new PullIterator() {
        // The index of the next element to pull.
        int index = 0;

        @Override
        public boolean pull() {
          if (index >= size()) {
            return false;
          }
          get(index, value);
          index++;
          return true;
        }
      };
    }

    /**
     * Copies the value at the given fromIndex to the value at the given toIndex. Throws
     * IndexOutOfBoundsException if either index is out of range.
     */
    default void copy(int fromIndex, int toIndex) {
      T tmp = newElement();
      get(fromIndex, tmp);
      set(toIndex, tmp);
      destroyElement(tmp);
    }

    /**
     * Swaps the elements at index 'a' and index 'b'. Throws IndexOutOfBoundsException if either
     * index is out of range.
     */
    default void swap(int indexA, int indexB) {
      T tmp = newElement();
      get(indexA, tmp);
      copy(indexB, indexA);
      set(indexB, tmp);
      destroyElement(tmp);
    }

    /** Sorts the list using the given comparator. This is not necessarily a stable sort. */
    default void sort(Comparator<T> comparator) {
      T left = newElement();
      T right = newElement();
      Sorter.sort(new SortableCollection() {
        @Override
        public int size() {
          return PullList.this.size();
        }

        @Override
        public void truncate(int end) {
          PullList.this.truncate(end);
        }

        @Override
        public boolean less(int leftIndex, int rightIndex) {
          PullList.this.get(leftIndex, left);
          PullList.this.get(rightIndex, right);
          return comparator.compare(left, right) < 0;
        }

        @Override
        public void swap(int indexA, int indexB) {
          PullList.this.swap(indexA, indexB);
        }
      });
      destroyElement(left);
      destroyElement(right);
    }

    /**
     * Makes this list have exactly 'newSize' elements. If the current size is smaller, adds new
     * elements to the end of the list. Truncates the list if the current size is greater. Does
     * nothing if the current size is equal to the requested newSize. The implementation may reuse
     * previously allocated storage and if so, is not required to reinitialize those values to
     * default values.
     *
     * <p>May throw UnsupportedOperationException if increasing capacity would be required but the
     * implementation does not support increasing capacity.
     */
    default void resize(int newSize) {
      if (newSize < size()) {
        truncate(newSize);
      } else if (newSize > size()) {
        enlarge(newSize);
      }
    }

    /**
     * Makes this list have at most 'newSize' elements. Does nothing if "newSize" is larger than or
     * equal to the current size(). Otherwise, removes entries from the list, from index "newSize"
     * to the end. Does not necessarily shrink the underlying data storage.
     */
    void truncate(int newSize);

    /**
     * Makes this list have at least 'newSize' elements. Does nothing if "newSize" is less than or
     * equal to the current size(). Otherwise, increases the size of the list. The implementation
     * may reuse previously allocated storage and is not required to reinitialize those values to
     * default values.
     *
     * <p>May throw UnsupportedOperationException if increasing capacity would be required but the
     * implementation does not support increasing capacity.
     */
    void enlarge(int newSize);

    /** Provides a List view of the PullList. Convenient but inefficient. */
    default List<T> asList() {
      return new AbstractList<T>() {

        @Override
        public T get(int index) {
          T value = newElement();
          PullList.this.get(index, value);
          return value;
        }

        @Override
        public T set(int index, T newValue) {
          T previousValue = get(index);
          PullList.this.set(index, newValue);
          return previousValue;
        }

        @Override
        public int size() {
          return PullList.this.size();
        }
      };
    }
  }
}
