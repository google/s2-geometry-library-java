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

import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.function.Supplier;

/**
 * A PooledList reuses elements to avoid allocations, increasing efficiency for lists that are
 * repeatedly cleared and re-filled. It contains elements of the parameterized type T, which must
 * be mutable. Objects in a PooledList are created by and owned by the PooledList, and are all the
 * same type, as opposed to a {@code List<T>} which can contain a mixture of different types that
 * all extend or implement T.
 *
 * <p>IMPORTANT: Avoid holding long-lived references to elements in a PooledList, as those objects
 * are owned by the PooledList and may be reused if the PooledList is modified. Specifically, if you
 * hold a reference to an element which is then removed from the PooledList by clear(), remove(),
 * removeLast(), or resize(), a subsequent call to add() may reuse the element. Then the reference
 * you are holding is to an object with different contents and a different position.
 *
 * <p>Likewise, avoid using elements in both PooledList and Java Collections, as this may result in
 * the collection holding references to elements that are unexpectedly modified when reused.
 *
 * <p>PooledList implements Iterable, but not List, as its mutation model is different. In
 * particular, PooledList does not implement an {@code add(T)} method that accepts an element to
 * add. Instead, it provides an {@link #add()} method that increases the size of the list by one,
 * and returns the added element, which is at the end of the list. Note that the previous content of
 * the element is _not_ cleared. To put elements into the PooledList, mutate elements returned by
 * {@link #add()} to the desired values. Elements returned by add() are reused if elements are
 * available in the pool, otherwise new elements are created using the supplier provided to the
 * constructor.
 *
 * <p>PooledList does provide {@link #asList()}, which returns a readonly List view for ease of
 * interoperability, but the elements in the returned List are still owned by the PooledList and may
 * be reused. Use with caution.
 */
public final class PooledList<T> implements Iterable<T> {
  /**
   * The elements are stored in an ArrayList, which may have a size greater than the number of
   * elements in the PooledList.
   */
  private final List<T> elements = new ArrayList<>();

  /** The supplier of new elements. */
  private final Supplier<T> supplier;

  /** The number of elements in the PooledList. */
  private int size = 0;

  /** Constructs a new PooledList that will use the given supplier to construct new elements. */
  public PooledList(Supplier<T> supplier) {
    this.supplier = supplier;
  }

  /**
   * Returns a read-only view of the elements in the PooledList as a List. IMPORTANT: As described
   * in the class javadoc, avoid keeping long-lived references to the List elements, as they are
   * reused if the PooledList is modified.
   */
  public List<T> asList() {
    return new AbstractList<T>() {
      @Override
      public T get(int index) {
        return PooledList.this.get(index);
      }

      @Override
      public int size() {
        return PooledList.this.size();
      }
    };
  }

  /**
   * Returns the item at the given index. IMPORTANT: As described in the class javadoc, avoid long-
   * lived references to the returned element, as it is reused if the PooledList is modified.
   */
  public T get(int index) {
    rangeCheck(index);
    return elements.get(index);
  }

  /**
   * Returns the item at the end of the list. IMPORTANT: As described in the class javadoc, avoid
   * long-lived references to the returned element, as it is reused if the PooledList is modified.
   */
  public T back() {
    return get(size() - 1);
  }

  /**
   * Returns the item at the beginning of the list. IMPORTANT: As described in the class javadoc,
   * avoid long-lived references to the returned element, as it is reused if the PooledList is
   * modified.
   */
  public T front() {
    return get(0);
  }

  /**
   * Increases the size of the list by one and returns the element at the end of the list. The
   * returned element should then be initialized. Reuses existing elements if available, otherwise
   * uses the supplier to construct a new element.
   *
   * <p>IMPORTANT: As described in the class javadoc, avoid long-lived references to the returned
   * element, as it is reused if the PooledList is modified.
   */
  public T add() {
    if (size == elements.size()) {
      T element = supplier.get();
      elements.add(element);
      size++;
      return element;
    }
    return elements.get(size++);
  }

  /** Removes the last element from the list, returning it to the pool for reuse. */
  public void removeLast() {
    size--;
  }

  /**
   * Removes the element at the given index from this PooledList, shifting any subsequent elements
   * to the left. The removed element is returned to the pool for reuse.
   */
  public void remove(int index) {
    rangeCheck(index);
    // TODO(torrey): If there are a lot of unused entries in the pool it would be worthwhile to
    // avoid shifting all of them, by shifting only the elements in the [index, size) range.
    T element = elements.remove(index);
    // Put the removed element back in the pool.
    elements.add(element);
    size--;
  }

  /** Returns the number of elements in the PooledList. */
  public int size() {
    return size;
  }

  /** Returns true if the PooledList has size zero. */
  public boolean isEmpty() {
    return size == 0;
  }

  /** Sets the size of the PooledList to zero, but keeps elements in the pool for reuse. */
  public void clear() {
    size = 0;
  }

  /**
   * Sets the size of the PooledList. If the new size is less than the current size, unneeded
   * elements are returned to the pool. If the new size is greater, existing elements in the pool
   * are reused if possible, otherwise new elements are created using the supplier. Note that if the
   * list size is increased, all the newly added elements need to be initialized.
   */
  public void resize(int size) {
    while (size > this.size) {
      T unused = add();
    }
    if (size < this.size) {
      this.size = size;
    }
  }

  /** Sorts the elements of the PooledList using the given comparator. */
  public void sort(Comparator<T> comparator) {
    Collections.sort(elements.subList(0, size), comparator);
  }

  /**
   * Throws an IndexOutOfBoundsException if the given index is not within the bounds of the list.
   */
  public void rangeCheck(int index) {
    if (index >= size || index < 0) {
      throw new IndexOutOfBoundsException("Index: " + index + ", Size: " + size);
    }
  }

  /** Swaps the elements at the given indices by exchanging references. */
  public void swap(int leftIndex, int rightIndex) {
    rangeCheck(leftIndex);
    rangeCheck(rightIndex);
    Collections.swap(elements, leftIndex, rightIndex);
  }

  /**
   * Returns an iterator over the elements in the PooledList. IMPORTANT: As described in the class
   * javadoc, avoid long-lived references to these elements, as they are reused if the PooledList is
   * modified.
   */
  @Override
  public Iterator<T> iterator() {
    return new Iterator<T>() {
      private int index = 0;

      @Override
      public boolean hasNext() {
        return index < size;
      }

      @Override
      public T next() {
        if (index >= size) {
          throw new NoSuchElementException();
        }
        return elements.get(index++);
      }
    };
  }
}
