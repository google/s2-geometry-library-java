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

import static java.util.Comparator.naturalOrder;

import com.google.common.geometry.S2Iterator.AvlSetIterator;
import com.google.common.primitives.UnsignedLongs;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import it.unimi.dsi.fastutil.objects.ObjectAVLTreeSet;
import java.util.Comparator;
import java.util.Objects;
import org.jspecify.annotations.Nullable;

/**
 * S2PointIndex maintains an index of points. Each point has some associated client-supplied data,
 * such as an index or object the point was taken from, useful to map query results back to another
 * data structure. The data may be null. Currently the index is an ordered set, sorted first by the
 * leaf S2CellId corresponding to the point, then the point, then the data.
 *
 * <p>The type D of the data is not required to be Comparable, but if it is, you can construct the
 * S2PointIndex with {@code S2PointIndex.forComparableData()}, and then identical points with
 * different non-null data values will be ordered accordingly. Otherwise, Entries with identical
 * points are sorted nulls-first, and otherwise have undefined order.
 *
 * <p>The class supports adding or removing points dynamically, and provides a seekable iterator
 * interface for navigating the index.
 *
 * <p>You can use this class in conjunction with {@link S2ClosestPointQuery} to find the closest
 * index points to a given query point. For example:
 *
 * {@snippet :
 * void test(List<S2Point> points, S2Point target) {
 *   // The generic type allows auxiliary data to be attached to each point
 *   // In this case, attach the original index of the point.
 *   S2PointIndex<Integer> index = new S2PointIndex();
 *   for (int i = 0; i < points.size(); i++) {
 *     index.add(points.get(i), i);
 *   }
 *   S2ClosestPointQuery<Integer> query = new S2ClosestPointQuery<>(index);
 *   query.findClosestPoint(target);
 *   if (query.num_points() > 0) {
 *     // query.point(0) is the closest point (result 0).
 *     // query.distance(0) is the distance to the target.
 *     // query.data(0) is the auxiliary data (the array index set above).
 *     doSomething(query.point(0), query.data(0), query.distance(0));
 *   }
 * }
 * }
 *
 * <p>Alternatively, you can access the index directly using the iterator interface. For example,
 * here is how to iterate through all the points in a given S2CellId "targetId":
 *
 * {@snippet :
 * S2Iterator<S2PointIndex.Entry<Integer>> it = index.iterator();
 * it.seek(targetId.rangeMin());
 * for (; !it.done() && it.compareTo(targetId.rangeMax()) <= 0; it.next()) {
 *   doSomething(it.entry());
 * }
 * }
 *
 * <p>Points can be added or removed from the index at any time by calling add() or remove(), but
 * doing so invalidates existing iterators. New iterators must be created.
 *
 * <p>This class is not thread-safe. TODO(user): Provide a thread-safe S2PointIndex.
 * TODO(user): Make this a subtype of S2Region, so that it can also be used to efficiently
 * compute coverings of a collection of S2Points.
 */
public final class S2PointIndex<D> {
  private final ObjectAVLTreeSet<Entry<D>> entries;

  /**
   * Creates a new S2PointIndex for a Data type D that is Comparable. If D is not Comparable, use
   * the default constructor instead.
   */
  public static <D extends Comparable<D>> S2PointIndex<D> forComparableData() {
    return new S2PointIndex<D>(Entry.stableOrder());
  }

  /** Creates a new S2PointIndex. */
  public S2PointIndex() {
    entries = new ObjectAVLTreeSet<>(Entry.order());
  }

  /**
   * Creates a new S2PointIndex that compares Entries with the given Comparator. The only intended
   * use is by {@link #forComparableData()}.
   */
  private S2PointIndex(Comparator<Entry<D>> comparator) {
    entries = new ObjectAVLTreeSet<>(comparator);
  }

  /** Returns the number of points in the index. */
  public int numPoints() {
    return entries.size();
  }

  /** Returns true if the index is empty. */
  public boolean isEmpty() {
    return entries.isEmpty();
  }

  /** Returns a new iterator over the cells of this index. */
  public S2Iterator<Entry<D>> iterator() {
    return new AvlSetIterator<>(entries, id -> new Entry<>(id, null, null));
  }

  /**
   * @deprecated No longer needed and has no effect.
   */
  @Deprecated
  public void applyUpdates() {}

  /** As {@link #add(Entry)}, but more convenient. */
  public void add(S2Point point, D data) {
    add(createEntry(point, data));
  }

  /** Adds a new entry to the index. Invalidates all iterators; clients must create new ones. */
  public void add(Entry<D> entry) {
    entries.add(entry);
  }

  /** As {@link #remove(Entry)}, but more convenient. */
  @CanIgnoreReturnValue
  public boolean remove(S2Point point, D data) {
    return remove(createEntry(point, data));
  }

  /**
   * Removes the given entry from the index, and returns whether the given entry was present and
   * removed. Both the "point" and "data" fields must match the point to be removed. Invalidates all
   * iterators; clients must create new ones.
   */
  @CanIgnoreReturnValue
  public boolean remove(Entry<D> entry) {
    return entries.remove(entry);
  }

  /**
   * Resets the index to its original empty state. Invalidates all iterators; clients must create
   * new ones.
   */
  public void reset() {
    entries.clear();
  }

  /**
   * Convenience method to create an index entry from the given point and data value. The data may
   * be null, but the point must not be.
   */
  public static <D> Entry<D> createEntry(S2Point point, D data) {
    return new Entry<>(S2CellId.fromPoint(point), point, data);
  }

  /**
   * An S2Iterator-compatible pair of S2Point with associated client data of a given type.
   *
   * <p>Equality and hashing are based on the point and data value. The natural order of this type
   * is only by the leaf cell that contains the point, which is <strong>not</strong> consistent with
   * equals. The provided Comparators {@link #order()} and {@link #stableOrder()} also consider the
   * point and data values and are consistent with equals.
   */
  public static class Entry<D> implements S2Iterator.Entry, Comparable<Entry<D>> {

    /**
     * A comparator for {@code Entry<D>} when D is Comparable. Considers the id, point and data
     * values in that order. Sorting with this Comparator is stable. See also {@link #order()},
     * which does not require D to be Comparable.
     */
    public static <D extends Comparable<D>> Comparator<Entry<D>> stableOrder() {
      return (a, b) -> compare(a, b, naturalOrder());
    }

    /**
     * A comparator for {@code Entry<D>} when D is not Comparable. Considers the id, point and data
     * values in that order. If data values are considered, they are only checked for equality, so
     * if multiple Entries with the same point and id have different data values, sorting with this
     * comparator will be unstable. See also {@link #stableOrder()}, which requires D to be
     * Comparable.
     */
    public static <D> Comparator<Entry<D>> order() {
      return (a, b) -> compare(a, b, (aData, bData) -> aData.equals(bData) ? 0 : -1);
    }

    private final long id;
    private final S2Point point;
    private final @Nullable D data;

    public Entry(S2CellId cellId, S2Point point, D data) {
      this.id = cellId.id();
      this.point = point;
      this.data = data;
    }

    @Override
    public long id() {
      return id;
    }

    public S2Point point() {
      return point;
    }

    public D data() {
      return data;
    }

    @Override
    public boolean equals(Object other) {
      if (other instanceof Entry) {
        Entry<?> e = (Entry<?>) other;
        return point.equalsPoint(e.point) && Objects.equals(data, e.data);
      } else {
        return false;
      }
    }

    @Override
    public int hashCode() {
      return point.hashCode() * 31 + (data == null ? 0 : data.hashCode());
    }

    @Override
    public int compareTo(Entry<D> other) {
      return UnsignedLongs.compare(id, other.id);
    }

    @Override
    public String toString() {
      return point.toDegreesString() + " : " + data;
    }

    /**
     * Helper for {@link #order()} and {@link #stableOrder()}. Compares id, point, and data values
     * in that order, finally using the given Comparator to break ties if a.data and b.data are both
     * non-null.
     */
    @SuppressWarnings("ReferenceEquality")
    private static <D> int compare(Entry<D> a, Entry<D> b, Comparator<D> comparator) {
      if (a == b) {
        return 0;
      }

      int cmp = UnsignedLongs.compare(a.id, b.id);
      if (cmp != 0) {
        return cmp;
      }

      if (a.point == null && b.point == null) {
        cmp = 0;
      } else if (a.point == null) {
        return -1;
      } else if (b.point == null) {
        return 1;
      } else {
        cmp = a.point.compareTo(b.point);
      }
      if (cmp != 0) {
        return cmp;
      }

      if (a.data == null && b.data == null) {
        return 0;
      } else if (a.data == null) {
        return -1;
      } else if (b.data == null) {
        return 1;
      }
      return comparator.compare(a.data, b.data);
    }
  }
}
