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

import com.google.common.base.Objects;
import com.google.common.collect.Lists;
import com.google.common.primitives.UnsignedLongs;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.Collections;
import java.util.List;

/**
 * S2PointIndex maintains an index of points sorted by leaf S2CellId. Each point has some associated
 * client-supplied data, such as an index or object the point was taken from, useful to map query
 * results back to another data structure.
 *
 * <p>The class supports adding or removing points dynamically, and provides a seekable iterator
 * interface for navigating the index.
 *
 * <p>You can use this class in conjunction with {@link S2ClosestPointQuery} to find the closest
 * index points to a given query point. For example:
 *
 * <pre>{@code
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
 * }</pre>
 *
 * <p>Alternatively, you can access the index directly using the iterator interface. For example,
 * here is how to iterate through all the points in a given S2CellId "targetId":
 *
 * <pre>{@code
 * S2Iterator<S2PointIndex.Entry<Integer>> it = index.iterator();
 * it.seek(targetId.rangeMin());
 * for (; !it.done() && it.compareTo(targetId.rangeMax()) <= 0; it.next()) {
 *   doSomething(it.entry());
 * }
 * }</pre>
 *
 * <p>Points can be added or removed from the index at any time by calling add() or remove(), but
 * doing so invalidates existing iterators. New iterators must be created.
 *
 * <p>This class is not thread-safe.
 * TODO(user): Make this a subtype of S2Region, so that it can also be used to efficiently
 * compute coverings of a collection of S2Points.
 */
public final class S2PointIndex<D> {
  private final List<Entry<D>> entries = Lists.newArrayList();
  private boolean sorted = true;

  /** Returns the number of points in the index. */
  public int numPoints() {
    return entries.size();
  }

  /**
   * Returns a new iterator over the cells of this index, after sorting entries by cell ID if any
   * modifications have been made since the last iterator was created.
   */
  public S2Iterator<Entry<D>> iterator() {
    if (!sorted) {
      Collections.sort(entries);
      sorted = true;
    }
    return S2Iterator.create(entries);
  }

  /** As {@link #add(Entry)}, but more convenient. */
  public void add(S2Point point, D data) {
    add(createEntry(point, data));
  }

  /** Adds a new entry to the index. Invalidates all iterators; clients must create new ones. */
  public void add(Entry<D> entry) {
    sorted = false;
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
    sorted = true;
    entries.clear();
  }

  /** Convenience method to create an index entry from the given point and data value. */
  public static <D> Entry<D> createEntry(S2Point point, D data) {
    return new Entry<>(S2CellId.fromPoint(point), point, data);
  }

  /**
   * An S2Iterator-compatible pair of S2Point with associated client data of a given type.
   *
   * <p>Equality and hashing are based on the point and data value. The natural order of this type
   * is by the leaf cell that contains the point, which is <strong>not</strong> consistent with
   * equals.
   */
  public static class Entry<D> implements S2Iterator.Entry, Comparable<Entry<D>> {
    private final long id;
    private final S2Point point;
    private final D data;

    private Entry(S2CellId cellId, S2Point point, D data) {
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
        return point.equalsPoint(e.point) && Objects.equal(data, e.data);
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
      return new S2LatLng(point) + ": " + data;
    }
  }
}
