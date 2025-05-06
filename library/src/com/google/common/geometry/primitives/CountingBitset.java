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

import java.util.BitSet;

/**
 * Similar to java.util.Bitset, but has a specific size and maintains a count of the total number of
 * bits that have been set. Unlike java.util.Bitset, CountingBitset does not grow on demand, but
 * may be reset to a new size. Having a specific size allows CountingBitset to provide a
 * {@link #full()} method indicating that all bits are set. (java.util.BitSet grows underlying
 * storage on demand, and provides the appearance of a set of indefinite size, limited only by the
 * range of Integer and available memory).
 *
 * <p>The count is idempotent (i.e. setting the same bit twice will only add one to the count), so
 * it's safe to use to track visiting a set of e.g. edges in an out-of-order fashion.
 */
public class CountingBitset {
  private final BitSet bits;
  private int size = 0;
  private int numSet = 0;

  /** Creates a new CountingBitset with size zero. Must be reset to a new size before use. */
  public CountingBitset() {
    this(0);
  }

  /** Creates a new CountingBitset with the given size in bits. All bits are initially false. */
  public CountingBitset(int size) {
    bits = new BitSet(size);
    this.size = size;
  }

  /** Clears the content of this CountingBitset and resets the size to the given size. */
  public void reset(int size) {
    clear();
    this.size = size;
  }

  /** Clears all set bits in this CountingBitset, retaining the size.*/
  public void clear() {
    bits.clear();
    numSet = 0;
  }

  /** Returns the size of this CountingBitset in bits. */
  public int bits() {
    return size;
  }

  /** Returns the number of bits set in this CountingBitset. */
  public int numSet() {
    return numSet;
  }

  /** Returns true if this CountingBitset has no bits set. */
  public boolean empty() {
    return numSet == 0;
  }

  /** Returns true if this CountingBitset has all bits set. */
  public boolean full() {
    return numSet == size;
  }

  /**
   * Returns the current value of the bit at the given bitIndex.
   *
   * @throws IndexOutOfBoundsException if bitIndex is negative or greater than or equal to the size.
   */
  public boolean get(int bitIndex) {
    if (bitIndex >= size) {
      throw new IndexOutOfBoundsException("bitIndex >= size: " + bitIndex);
    }
    // The BitSet checks that bitIndex is non-negative and throws IndexOutOfBoundsException if not.
    return bits.get(bitIndex);
  }

  /**
   * Sets the value of the bit at the given bitIndex to true. If we're changing the bit's state,
   * then the count is modified, otherwise it's left as is. Thus this function is idempotent.
   *
   * @throws IndexOutOfBoundsException if bitIndex is negative or greater than or equal to the size.
   */
  public void set(int bitIndex) {
    set(bitIndex, true);
  }

  /**
   * Sets the value of the bit at the given bitIndex to the given value. If we're changing the bit's
   * state, then the count is modified, otherwise it's left as is. Thus this function is idempotent.
   *
   * @throws IndexOutOfBoundsException if bitIndex is negative or greater than or equal to the size.
   */
  public void set(int bitIndex, boolean value) {
    if (bitIndex >= size) {
      throw new IndexOutOfBoundsException("bitIndex >= size: " + bitIndex);
    }
    // The BitSet checks that bitIndex is non-negative and throws IndexOutOfBoundsException if not.
    boolean old = bits.get(bitIndex);
    if (old != value) {
      bits.set(bitIndex, value);
      numSet += value ? +1 : -1;
    }
  }
}
