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
package com.google.common.geometry;

import com.google.common.geometry.S2ShapeIndex.CellRelation;
import com.google.errorprone.annotations.CanIgnoreReturnValue;

/**
 * S2CellRangeIterator is a wrapper around an S2Iterator that caches the rangeMin() and rangeMax()
 * that each cell covers as it iterates. This lets us define range based methods such as seekTo()
 * and locate() efficiently.
 *
 * <p>Computing the rangeMax() and rangeMin() of a cell isn't expensive but it's not free either, so
 * we extend the S2Iterator interface instead of integrating this functionality there, allowing the
 * user to pay only for what they use.
 *
 * <p>An S2CellRangeIterator wraps an S2Iterator, but is also itself an S2Iterator and thus can be
 * used anywhere one is required.
 *
 * <p>TODO(user): Remove S2ShapeUtil.RangeIterator, which is an older version of this.
 */
public class S2CellRangeIterator<T extends S2Iterator.Entry> implements S2Iterator<T> {
  private final S2Iterator<T> it;
  private S2CellId rangeMin;
  private S2CellId rangeMax;

  /** Constructs a new S2CellRangeIterator positioned at the beginning. */
  public S2CellRangeIterator(S2Iterator<T> it) {
    this.it = it;
    restart();
  }

  /**
   * Builds a new S2CellRangeIterator from an S2ShapeIndex.
   *
   * <p>The index must live for the duration of the iterator.
   */
  public static S2CellRangeIterator<S2ShapeIndex.Cell> makeS2CellRangeIterator(S2ShapeIndex index) {
    return new S2CellRangeIterator<S2ShapeIndex.Cell>(index.iterator());
  }

  /** Builds a new S2CellRangeIterator from an S2Iterator */
  public static <E extends S2Iterator.Entry> S2CellRangeIterator<E> makeS2CellRangeIterator(
      S2Iterator<E> iter) {
    return new S2CellRangeIterator<>(iter);
  }

  @Override
  public S2CellRangeIterator<T> copy() {
    return new S2CellRangeIterator<>(it.copy());
  }

  /** Returns the underlying S2Iterator. */
  public S2Iterator<T> iterator() {
    return it;
  }

  /** The range min of the current cell. */
  public S2CellId rangeMin() {
    return rangeMin;
  }

  /** The range max of the current cell. */
  public S2CellId rangeMax() {
    return rangeMax;
  }

  /** For convenience and to match the C++ API, an alias for restart(). */
  public void begin() {
    restart();
  }

  @Override // S2Iterator implementation
  public void restart() {
    it.restart();
    refresh();
  }

  @Override // S2Iterator implementation
  public S2CellId id() {
    return it.id();
  }

  @Override // S2Iterator implementation
  public T entry() {
    return it.entry();
  }

  @Override // S2Iterator implementation
  @CanIgnoreReturnValue
  public boolean next() {
    boolean result = it.next();
    refresh();
    return result;
  }

  @Override // S2Iterator implementation
  @CanIgnoreReturnValue
  public boolean prev() {
    boolean result = it.prev();
    refresh();
    return result;
  }

  @Override // S2Iterator implementation
  public boolean done() {
    return it.done();
  }

  @Override // S2Iterator implementation
  public boolean atBegin() {
    return it.atBegin();
  }

  @Override // S2Iterator implementation
  public void seek(S2CellId target) {
    it.seek(target);
    refresh();
  }

  @Override // S2Iterator implementation
  public void seekForward(S2CellId target) {
    it.seekForward(target);
    refresh();
  }

  @Override // S2Iterator implementation
  public void finish() {
    it.finish();
    refresh();
  }

  @Override
  public boolean locate(S2Point target) {
    boolean result = it.locate(target);
    refresh();
    return result;
  }

  @Override // S2Iterator implementation
  public CellRelation locate(S2CellId target) {
    // Let T be the target cell id, let I = seek(T.rangeMin()) and let prev(I) be the predecessor of
    // I. If T contains any index cells, then T contains I. Similarly, if T is contained by an index
    // cell, then the containing cell is either I or prev(I). We test for containment by comparing
    // the ranges of leaf cells spanned by T, I, and prev(I).
    seek(target.rangeMin());
    if (!done()) {
      // The target is contained by the cell we landed on, so it's indexed.
      if (id().greaterOrEquals(target) && rangeMin().lessOrEquals(target)) {
        return CellRelation.INDEXED;
      }

      // The cell we landed on is contained by the target, so it's subdivided.
      if (id().lessOrEquals(target.rangeMax())) {
        return CellRelation.SUBDIVIDED;
      }
    }

    // Otherwise check the previous cell (if it exists). If it contains the target then it's
    // indexed, otherwise the target cell is disjoint.
    if (prev() && rangeMax().greaterOrEquals(target)) {
      return CellRelation.INDEXED;
    }
    return CellRelation.DISJOINT;
  }

  /** Convenience re-implementation of the above function, see it for details. */
  public CellRelation locate(S2CellRangeIterator<T> target) {
    seek(target.rangeMin());
    if (!done()) {
      // The target is contained by the cell we landed on, so it's indexed.
      if (id().greaterOrEquals(target.id()) && rangeMin().lessOrEquals(target.id())) {
        return CellRelation.INDEXED;
      }

      // The cell we landed on is contained by the target, so it's subdivided.
      if (id().lessOrEquals(target.rangeMax())) {
        return CellRelation.SUBDIVIDED;
      }
    }

    // Otherwise check the previous cell (if it exists). If it contains the target then it's
    // indexed, otherwise the target cell is disjoint.
    if (prev() && rangeMax().greaterOrEquals(target.id())) {
      return CellRelation.INDEXED;
    }
    return CellRelation.DISJOINT;
  }

  public void seekTo(S2CellRangeIterator<?> target) {
    seek(target.rangeMin());

    // If the current cell does not overlap "target", it is possible that the previous cell is the
    // one we are looking for. This can only happen when the previous cell contains "target" but has
    // a smaller S2CellId.
    if (done() || rangeMin().greaterThan(target.rangeMax())) {
      if (prev() && rangeMax().lessThan(target.id())) {
        next();
      }
    }
    refresh();
  }

  public void seekBeyond(S2CellRangeIterator<?> target) {
    seek(target.rangeMax().next());
    if (!done() && rangeMin().lessOrEquals(target.rangeMax())) {
      next();
    }
    refresh();
  }

  /**
   * Queries the relationship between two range iterators. Returns -1 if this iterator's current
   * position entirely precedes the other iterator's current position, +1 if it entirely follows,
   * and 0 if they overlap.
   */
  public int relation(S2CellRangeIterator<?> b) {
    if (rangeMax().lessThan(b.rangeMin())) {
      return -1;
    }
    if (rangeMin().greaterThan(b.rangeMax())) {
      return +1;
    }
    return 0;
  }

  private void refresh() {
    if (done()) {
      rangeMin = S2CellId.sentinel().rangeMin();
      rangeMax = S2CellId.sentinel().rangeMax();
    } else {
      rangeMin = id().rangeMin();
      rangeMax = id().rangeMax();
    }
  }
}
