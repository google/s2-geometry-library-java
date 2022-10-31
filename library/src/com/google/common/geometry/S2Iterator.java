/*
 * Copyright 2015 Google Inc.
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

import static java.lang.Math.max;

import com.google.common.base.Function;
import com.google.common.geometry.S2ShapeIndex.CellRelation;
import com.google.common.primitives.UnsignedLongs;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.List;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;

/**
 * A random access iterator that provides low-level access to entries sorted by cell ID. The
 * behavior of this iterator is more like a database cursor, where accessing properties at the
 * current position does not alter the position of the cursor. The cursor has a {@link #compareTo}
 * method to compare the value at the current position of the iterator with a given S2CellId.
 */
// TODO(user): Replace Entry<T> with a Multimap<long,T> that is space efficient, and supports
// time-efficient inserts, removes, and lookups.
@JsType
public final class S2Iterator<T extends S2Iterator.Entry> {
  /** An interface to provide the cell ID for an element in a sorted list. */
  @JsType
  public interface Entry {
    /** Returns the cell ID of this cell as a primitive. */
    long id();
  }

  /**
   * Creates an iterator given a list of entries. Package private and not public, since only S2
   * classes guarantee the necessary preconditions on {@code entries} -- that the cell IDs of each
   * entry are sorted in ascending order.
   */
  static <T extends S2Iterator.Entry> S2Iterator<T> create(List<T> entries) {
    return new S2Iterator<T>(entries);
  }

  /**
   * Same as {@link #create(List)}, but accepts {@code seekFunction}, which is used as the
   * implementation of {@link #seek(S2CellId)}.
   *
   * @param entries the list of entries which back this iterator.
   * @param seekFunction a function which takes a target {@link S2CellId} and returns an index to
   *     which this iterator will be repositioned.
   */
  static <T extends S2Iterator.Entry> S2Iterator<T> create(
      List<T> entries, Function<S2CellId, Integer> seekFunction) {
    return new S2Iterator<T>(entries, seekFunction);
  }

  /** Creates a new iterator with the same entries and position as {@code it}. */
  static <T extends S2Iterator.Entry> S2Iterator<T> copy(S2Iterator<T> it) {
    S2Iterator<T> copy = new S2Iterator<T>(it.entries, it.seekFunction);
    copy.pos = it.pos;
    return copy;
  }

  /** Returns a copy of this iterator, positioned as this iterator is. */
  public S2Iterator<T> copy() {
    S2Iterator<T> it = new S2Iterator<T>(entries, seekFunction);
    it.pos = pos;
    return it;
  }

  private final List<T> entries;
  private final Function<S2CellId, Integer> seekFunction;
  private int pos;

  /**
   * Create a new iterator based on the given list of entries. Results are undefined if the entries
   * are not in ascending sorted order.
   *
   * @param entries the list of entries which back this iterator.
   */
  S2Iterator(List<T> entries) {
    this.entries = entries;
    this.seekFunction =
        (target) -> {
          int start = 0;
          int end = entries.size() - 1;
          while (start <= end) {
            int mid = (start + end) / 2;
            long id = entries.get(mid).id();
            int result = UnsignedLongs.compare(id, target.id());
            if (result > 0) {
              end = mid - 1;
            } else if (result < 0) {
              start = mid + 1;
            } else if (start != mid) {
              end = mid;
            } else {
              return mid;
            }
          }
          return start;
        };
  }

  /**
   * Same as {@link #S2Iterator(List)}, but accepts {@code seekFunction}, which is used as the
   * implementation of {@link #seek(S2CellId)}.
   *
   * @param entries the list of entries which back this iterator.
   * @param seekFunction a function which takes a target {@link S2CellId} and returns an index to
   *     which this iterator will be repositioned.
   */
  S2Iterator(List<T> entries, Function<S2CellId, Integer> seekFunction) {
    this.entries = entries;
    this.seekFunction = seekFunction;
  }

  /** Positions the iterator so that {@link #atBegin()} is true. */
  public void restart() {
    pos = 0;
  }

  /** Returns the comparison from the current iterator cell to the given cell ID. */
  public int compareTo(S2CellId cellId) {
    return UnsignedLongs.compare(entry().id(), cellId.id());
  }

  /** Returns true if {@code o} is an {@link S2Iterator} with equal entries and position. */
  @Override
  public boolean equals(Object o) {
    return o instanceof S2Iterator<?> && equalIterators((S2Iterator<?>) o);
  }

  @Override
  public int hashCode() {
    return 31 * pos + entries.hashCode();
  }

  /** Returns true if these iterators have the same entries and position. */
  public boolean equalIterators(S2Iterator<?> it) {
    return entries == it.entries && pos == it.pos;
  }

  /** Returns the cell id for the current cell. */
  public S2CellId id() {
    return new S2CellId(entry().id());
  }

  /** Returns the current entry. */
  public T entry() {
    // assert !done();
    return entries.get(pos);
  }

  /** Returns the center of the cell (used as a reference point for shape interiors.) */
  public S2Point center() {
    return id().toPoint();
  }

  /**
   * If {@code pos} is equal to the number of cells in the index, does not move the iterator, and
   * returns false. Otherwise, advances the iterator to the next cell in the index and returns true.
   */
  @CanIgnoreReturnValue
  public boolean next() {
    if (pos < entries.size()) {
      pos++;
      return true;
    }
    return false;
  }

  /**
   * If {@code pos} is equal to 0, does not move the iterator and returns false. Otherwise,
   * positions the iterator at the previous cell in the index and returns true.
   */
  @CanIgnoreReturnValue
  public boolean prev() {
    if (pos > 0) {
      pos--;
      return true;
    }
    return false;
  }

  /** Returns true if the iterator is positioned past the last index cell. */
  public boolean done() {
    return pos == entries.size();
  }

  /** Returns true if the iterator is positioned at the first index cell. */
  public boolean atBegin() {
    return pos == 0;
  }

  /**
   * Positions the iterator at the first cell with {@code id() >= target}, or at the end of the
   * index if no such cell exists.
   */
  public void seek(S2CellId target) {
    pos = seekFunction.apply(target);
  }

  /**
   * Advances the iterator to the next cell with {@code id() >= target}. If the iterator is
   * {@link #done()} or already satisfies {@code id() >= target}, there is no effect.
   */
  public void seekForward(S2CellId target) {
    if (!done() && compareTo(target) < 0) {
      int tmpPos = pos;
      seek(target);
      pos = max(pos, tmpPos + 1);
    }
  }

  /** Positions the iterator so that {@link #done()} is true. */
  public void finish() {
    pos = entries.size();
  }

  /**
   * Positions the iterator at the index cell containing "target" and returns true, or if no such
   * cell exists in the index, the iterator is positioned arbitrarily and this method returns false.
   *
   * <p>The resulting index position is guaranteed to contain all edges that might intersect the
   * line segment between {@code targetPoint} and {@link #center()}.
   */
  @JsIgnore // No method overloading for J2CL.
  public boolean locate(S2Point targetPoint) {
    // Let I be the first cell not less than T, where T is the leaf cell containing "targetPoint".
    // Then if T is contained by an index cell, then the containing cell is either I or I'. We
    // test for containment by comparing the ranges of leaf cells spanned by T, I, and I'.
    S2CellId target = S2CellId.fromPoint(targetPoint);
    seek(target);
    if (!done() && id().rangeMin().lessOrEquals(target)) {
      return true;
    }
    if (!atBegin()) {
      prev();
      if (id().rangeMax().greaterOrEquals(target)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Positions the iterator at the index cell containing the given cell, if possible, and returns
   * the {@link CellRelation} that describes the relationship between the index and the given target
   * cell:
   *
   * <ul>
   *   <li>Returns {@link CellRelation#INDEXED} if the iterator was positioned at an index cell that
   *       is equal to or contains the given cell. I.e. the given target exists in the index as a
   *       leaf cell.
   *   <li>Returns {@link CellRelation#SUBDIVIDED} if the iterator was positioned at the first of
   *       one or more cells contained by the given target cell. I.e. the target does not exist in
   *       the index, but the first of its descendants was selected.
   *   <li>Returns {@link CellRelation#DISJOINT} if the iterator had to be positioned arbitrarily
   *       because the given target cell does not intersect any of the index's cells.
   * </ul>
   */
  public CellRelation locate(S2CellId target) {
    // Let T be the target, let I be the first cell not less than T.rangeMin(), and let I' be the
    // predecessor of I.  If T contains any index cells, then T contains I.  Similarly, if T is
    // contained by an index cell, then the containing cell is either I or I'.  We test for
    // containment by comparing the ranges of leaf cells spanned by T, I, and I'.
    seek(target.rangeMin());
    if (!done()) {
      if (id().greaterOrEquals(target) && id().rangeMin().lessOrEquals(target)) {
        return CellRelation.INDEXED;
      }
      if (id().lessOrEquals(target.rangeMax())) {
        return CellRelation.SUBDIVIDED;
      }
    }
    if (!atBegin()) {
      prev();
      if (id().rangeMax().greaterOrEquals(target)) {
        return CellRelation.INDEXED;
      }
    }
    return CellRelation.DISJOINT;
  }

  /** Set this iterator to the position given by the other iterator. */
  public void position(S2Iterator<T> it) {
    // assert entries == it.entries;
    pos = it.pos;
  }
}
