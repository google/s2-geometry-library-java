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

import com.google.common.base.Function;
import com.google.common.geometry.S2ShapeIndex.CellRelation;
import com.google.common.primitives.UnsignedLongs;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;
import it.unimi.dsi.fastutil.objects.ObjectBidirectionalIterator;
import java.util.List;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;

/**
 * A random access iterator that provides low-level access to entries sorted by cell ID. The
 * behavior of this iterator is more like a database cursor, where accessing properties at the
 * current position does not alter the position of the cursor. The cursor has a {@link #compareTo}
 * method to compare the value at the current position of the iterator with a given S2CellId.
 */
@JsType
@SuppressWarnings("Assertion")
public interface S2Iterator<T extends S2Iterator.Entry> {
  /** An interface to provide the cell ID for an element in a sorted list. */
  @JsType
  interface Entry {
    /** Returns the cell ID of this cell as a primitive. */
    long id();
  }

  /**
   * Creates an iterator given a list of entries. Package private and not public, since only S2
   * classes guarantee the necessary preconditions on {@code entries} -- that the cell IDs of each
   * entry are sorted in ascending order.
   */
  static <T extends S2Iterator.Entry> ListIterator<T> fromList(List<T> entries) {
    return new ListIterator<>(entries);
  }

  /** Returns a copy of this iterator, positioned as this iterator is. */
  S2Iterator<T> copy();

  /** An {@link S2Iterator} based on a list. */
  class ListIterator<T extends Entry> implements S2Iterator<T> {
    private final List<T> entries;
    protected int pos;

    /**
     * Create a new iterator based on the given list of entries. Results are undefined if the
     * entries are not in ascending sorted order.
     *
     * @param entries the list of entries which back this iterator.
     */
    ListIterator(List<T> entries) {
      this.entries = entries;
    }

    @Override
    public S2Iterator<T> copy() {
      ListIterator<T> it = new ListIterator<T>(entries);
      it.pos = pos;
      return it;
    }

    @Override
    public void restart() {
      pos = 0;
    }

    @Override
    public T entry() {
      assert !done();
      return entries.get(pos);
    }

    @Override
    public boolean next() {
      if (pos < entries.size()) {
        pos++;
        return true;
      }
      return false;
    }

    @Override
    public boolean prev() {
      if (pos > 0) {
        pos--;
        return true;
      }
      return false;
    }

    @Override
    public boolean done() {
      return pos == entries.size();
    }

    @Override
    public boolean atBegin() {
      return pos == 0;
    }

    @Override
    public void seek(S2CellId target) {
      seek(0, target);
    }

    @Override
    public void seekForward(S2CellId target) {
      seek(pos, target);
    }

    private void seek(int start, S2CellId target) {
      int end = entries.size() - 1;
      while (start <= end) {
        pos = (start + end) / 2;
        int result = id().compareTo(target);
        if (result > 0) {
          end = pos - 1;
        } else if (result < 0) {
          start = pos + 1;
        } else if (start != pos) {
          end = pos;
        } else {
          return;
        }
      }
      pos = start;
    }

    @Override
    public void finish() {
      pos = entries.size();
    }
  }

  /** * An {@link S2Iterator} implemented on an {@link ObjectAVLTreeSet AVL tree}. */
  class AvlSetIterator<T extends Entry> implements S2Iterator<T> {
    /** The entries of this set. */
    private final ObjectAVLTreeSet<T> entries;

    /** Returns an entry that is less than all other entries with the given cell. */
    private final Function<S2CellId, T> min;

    /** The current unbounded iterator position in this set. Positioned after 'entry'. */
    private ObjectBidirectionalIterator<T> setIter;

    /** The current entry of this S2Iterator. Null if the iterator is done. */
    private T entry;

    /**
     * Creates a new AVL tree-based S2Iterator.
     *
     * @param entries the entries to wrap
     * @param min returns an entry that is less than all other entries with the given cell
     */
    public AvlSetIterator(ObjectAVLTreeSet<T> entries, Function<S2CellId, T> min) {
      this.entries = entries;
      this.min = min;
      restart();
    }

    @Override
    public void restart() {
      setIter = entries.iterator(); // at this point, 'it' has 'next' but 'curr' is null.
      entry = setIter.hasNext() ? setIter.next() : null;
    }

    @Override
    public void finish() {
      this.setIter = entries.isEmpty() ? entries.iterator() : entries.iterator(entries.last());
      entry = null;
    }

    /**
     * Positions this iterator at the given entry and remembers the current entry. The given entry
     * must be one of the entries in the set: otherwise an assertion error will be thrown (if
     * assertions are enabled).
     */
    private void set(T entry) {
      assert entries.contains(entry);
      this.setIter = entries.iterator(entry);
      this.entry = entry;
    }

    @Override
    public S2Iterator<T> copy() {
      AvlSetIterator<T> result = new AvlSetIterator<>(entries, min);
      result.set(entry);
      assert result.isEqualTo(this);
      return result;
    }

    @Override
    public T entry() {
      return entry;
    }

    // TODO(torrey): comparing entries with equals() is non-trivial, can it be avoided or faster?
    @Override
    public boolean atBegin() {
      return entries.isEmpty() || entries.first().equals(entry);
    }

    @Override
    public boolean next() {
      if (setIter.hasNext()) {
        entry = setIter.next();
        return true;
      }
      entry = null;
      return false;
    }

    /**
     * If there is a previous entry, move to it and return true. Otherwise, do not move and return
     * false.
     */
    @Override
    public boolean prev() {
      if (!setIter.hasPrevious()) {
        return false;
      }

      // done() is a special case.
      if (entry == null) {
        entry = setIter.previous();
        return true;
      }

      // The setIter is just after 'entry'. Position setIter just before 'entry' so we can find out
      // if there is a previous entry.
      setIter.previous();

      if (setIter.hasPrevious()) {
        // There is a previous entry. Get it, which positions setIter just before it.
        entry = setIter.previous();
        // Reposition setIter just after 'entry', where we want it.
        setIter.next();
        return true;
      } else {
        // No previous entry. Reposition setIter back to the original position and return false.
        setIter.next();
        return false;
      }
    }

    @Override
    public boolean done() {
      return entry == null;
    }

    @Override
    public void seek(S2CellId target) {
      this.setIter = entries.iterator(min.apply(target));
      this.entry = setIter.hasNext() ? setIter.next() : null;
    }

    /**
     * Returns true if this AvlSetIterator is equal to the given other AvlSetIterator: they are
     * iterating over the same AVLTreeSet instance, and are positioned at the same entry.
     */
    @SuppressWarnings("ReferenceEquality")
    public boolean isEqualTo(AvlSetIterator<T> other) {
      if (entries != other.entries) {
        return false;
      }
      if (min != other.min) {
        return false;
      }
      if ((entry == null) != (other.entry == null)) {
        return false;
      }
      if (entry != null && entry != other.entry) {
        return false;
      }
      if (setIter.hasNext() != other.setIter.hasNext()) {
        return false;
      }
      if (setIter.hasPrevious() != other.setIter.hasPrevious()) {
        return false;
      }
      return true;
    }
  }

  /** Positions the iterator so that {@link #atBegin()} is true. */
  void restart();

  /** Returns the comparison from the current iterator cell to the given cell ID. */
  default int compareTo(S2CellId cellId) {
    return UnsignedLongs.compare(entry().id(), cellId.id());
  }

  /** Returns the cell id for the current cell. */
  default S2CellId id() {
    return new S2CellId(entry().id());
  }

  /** Returns the current entry. */
  T entry();

  /** Returns the center of the cell (used as a reference point for shape interiors.) */
  default S2Point center() {
    return id().toPoint();
  }

  /**
   * If {@code pos} is equal to the number of cells in the index, does not move the iterator, and
   * returns false. Otherwise, advances the iterator to the next cell in the index and returns true.
   */
  @CanIgnoreReturnValue
  boolean next();

  /**
   * If {@code pos} is equal to 0, does not move the iterator and returns false. Otherwise,
   * positions the iterator at the previous cell in the index and returns true.
   */
  @CanIgnoreReturnValue
  boolean prev();

  /** Returns true if the iterator is positioned past the last index cell. */
  boolean done();

  /** Returns true if the iterator is positioned at the first index cell. */
  boolean atBegin();

  /**
   * Positions the iterator at the first cell with {@code id() >= target}, or at the end of the
   * index if no such cell exists.
   */
  void seek(S2CellId target);

  /**
   * Advances the iterator to the next cell with {@code id() >= target}. If the iterator is {@link
   * #done()} or already satisfies {@code id() >= target}, there is no effect.
   */
  default void seekForward(S2CellId target) {
    if (!done() && compareTo(target) < 0) {
      seek(target);
    }
  }

  /** Positions the iterator so that {@link #done()} is true. */
  void finish();

  /**
   * Positions the iterator at the index cell containing "target" and returns true, or if no such
   * cell exists in the index, the iterator is positioned arbitrarily and this method returns false.
   *
   * <p>The resulting index position is guaranteed to contain all edges that might intersect the
   * line segment between {@code targetPoint} and {@link #center()}.
   */
  @JsIgnore // No method overloading for J2CL.
  default boolean locate(S2Point targetPoint) {
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
  default CellRelation locate(S2CellId target) {
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
}
