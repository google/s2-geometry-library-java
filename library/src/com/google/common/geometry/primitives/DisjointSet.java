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

import com.google.common.base.Preconditions;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * A disjoint set (AKA union-find set, AKA merge-find set) is a data structure that stores a
 * partition of a set into disjoint subsets. It allows us to efficiently add, merge and find a
 * representative member of a set.
 *
 * <p>In practice, disjoint sets are used for connected-component like algorithms where we want to
 * build up relationships between subsets and query the final connected pieces.
 *
 * <p>This implementation uses both path compression and union-by-size, so both the union() and
 * findRoot() operations have O(a(N)) amortized complexity, where a(n) is the inverse Ackermann
 * function. For any value that we'll ever care about, a(n) is less than 5, which is, for all
 * practical purposes, constant.
 */
@SuppressWarnings("Assertion")
public class DisjointSet<T> {
  /** Map from elements to their indices in the elements, parents, and sizes lists. */
  private final Map<T, Integer> elementIndices;
  /** The inverse of elementIndices, maps from index back to element. */
  private final ArrayList<T> elements;
  /** Parallel array with elements: the index of the parent of each element. */
  private final IntVector parents;
  /** Parallel array with elements: the number of elements that have this element as their root. */
  private final IntVector sizes;

  /** Creates a new disjoint set with default capacity. */
  public DisjointSet() {
    elementIndices = new HashMap<>();
    elements = new ArrayList<>();
    parents = new IntVector();
    sizes = new IntVector();
  }

  /** Creates a new disjoint set with an initial capacity. */
  public DisjointSet(int capacity) {
    // To ensure that we don't have to resize the map, scale according to the default load factor.
    elementIndices = new HashMap<>((int) Math.ceil(capacity / 0.75));
    elements = new ArrayList<>(capacity);
    parents = IntVector.ofCapacity(capacity);
    sizes = IntVector.ofCapacity(capacity);
  }

  /** Reserves capacity to avoid repeated resizing. */
  public void reserve(int capacity) {
    // elementIndices is skipped as Map doesn't provide ensureCapacity.
    elements.ensureCapacity(capacity);
    parents.ensureCapacity(capacity);
    sizes.ensureCapacity(capacity);
  }

  /**
   * Adds a new element to the set. The parent of the element is itself, that is, it's disjoint from
   * all the other elements of the set.
   *
   * <p>If the element is already in the set, no changes are made and false is returned. Otherwise
   * returns true. Immediately after calling add(val), findRoot(val) == val.
   */
  @CanIgnoreReturnValue
  public boolean add(T val) {
    if (elementIndices.containsKey(val)) {
      return false;
    }

    // Add the element to the list of elements, and record its index.
    int elementIndex = elements.size();
    elements.add(val);
    elementIndices.put(val, elementIndex);

    // Each element is initially its own parent.
    parents.add(elementIndex);
    sizes.add(1);

    Preconditions.checkState(parents.get(elementIndex) == elementIndex);
    return true;
  }

  /**
   * Returns the root of the given element. If the element isn't in the set, returns null.
   */
  public T findRoot(T val) {
    if (!elementIndices.containsKey(val)) {
      return null;
    }
    // Get index of the given element, follow parent chain to root, and return that root element.
    return elements.get(findRootIndex(elementIndices.get(val)));
  }

  /** Returns the index of the root for the element at the given initialIndex. */
  private int findRootIndex(int elementIndex) {
    int parentIndex = parents.get(elementIndex);
    // If this element has itself as its parent, it is the root.
    if (parentIndex == elementIndex) {
      return parentIndex;
    }
    // Otherwise, follow the parent chain.
    int rootIndex = findRootIndex(parentIndex);
    parents.set(elementIndex, rootIndex); // Path compression.
    return rootIndex;
  }

  /**
   * Performs a Union of two subsets. Makes the parent of a's subset equal to the parent of b's
   * subset. If either a or b isn't in the set, returns false and does not modify the set. Otherwise
   * returns true.
   */
  @CanIgnoreReturnValue
  public boolean union(T a, T b) {
    if (!elementIndices.containsKey(a) || !elementIndices.containsKey(b)) {
      return false;
    }
    int aRootIndex = findRootIndex(elementIndices.get(a));
    int bRootIndex = findRootIndex(elementIndices.get(b));
    if (aRootIndex == bRootIndex) {
      return true; // Already in the same set.
    }
    // Try to balance the tree: Make the smaller subtree the child of the larger.
    if (sizes.get(aRootIndex) < sizes.get(bRootIndex)) {
      parents.set(aRootIndex, bRootIndex);
      sizes.incrementAt(bRootIndex, sizes.get(aRootIndex));
    } else {
      parents.set(bRootIndex, aRootIndex);
      sizes.incrementAt(aRootIndex, sizes.get(bRootIndex));
    }
    return true;
  }

  /** Returns the total number of elements in the set. */
  public int size() {
    return elements.size();
  }

  /** Clears the set so that size() == 0. */
  public void clear() {
    elementIndices.clear();
    elements.clear();
    parents.clear();
    sizes.clear();
  }
}
