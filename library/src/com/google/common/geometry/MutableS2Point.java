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

import com.google.common.geometry.primitives.Pullable.PullList;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.AbstractList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.List;

/**
 * MutableS2Point is an interface to an XYZ coordinate. They can be set from an S2Point or from
 * coordinates. These don't typically exist as individual objects, (although {@link
 * MutableS2PointImpl} provides a simple implementation), but are instead used to access the
 * contents of a {@link MutableS2PointList} which stores a list of points in a memory efficient
 * array.
 */
@SuppressWarnings("Assertion")
public interface MutableS2Point extends Comparable<MutableS2Point> {
  /** Returns the x coordinate of this point. */
  double getX();

  /** Returns the y coordinate of this point. */
  double getY();

  /** Returns the z coordinate of this point. */
  double getZ();

  /** Sets the coordinates of this point to equal the given S2Point. */
  void set(S2Point other);

  /** Sets the coordinates of this point to the given X, Y, and Z values. */
  void set(double x, double y, double z);

  /** Returns true if this MutableS2Point is currently equal to the given S2Point. */
  default boolean isEqualTo(S2Point other) {
    return getX() == other.getX() && getY() == other.getY() && getZ() == other.getZ();
  }

  /** Returns true if this MutableS2Point is currently equal to the given MutableS2Point. */
  default boolean isEqualTo(MutableS2Point other) {
    return getX() == other.getX() && getY() == other.getY() && getZ() == other.getZ();
  }

  /**
   * Returns the comparison of this point to the given point. Comparison is first by X, breaking
   * ties in X with Y, and then breaking ties in X and Y with Z.
   */
  @Override
  default int compareTo(MutableS2Point other) {
    if (getX() < other.getX()) {
      return -1;
    }
    if (other.getX() < getX()) {
      return 1;
    }
    if (getY() < other.getY()) {
      return -1;
    }
    if (other.getY() < getY()) {
      return 1;
    }
    if (getZ() < other.getZ()) {
      return -1;
    }
    if (other.getZ() < getZ()) {
      return 1;
    }
    return 0;
  }

  /** A trivial implementation of MutableS2Point. */
  public static class MutableS2PointImpl implements MutableS2Point {
    private double x;
    private double y;
    private double z;

    @Override
    public void set(double x, double y, double z) {
      this.x = x;
      this.y = y;
      this.z = z;
    }

    @Override
    public void set(S2Point other) {
      x = other.getX();
      y = other.getY();
      z = other.getZ();
    }

    @Override
    public double getZ() {
      return z;
    }

    @Override
    public double getY() {
      return y;
    }

    @Override
    public double getX() {
      return x;
    }

    @Override
    public String toString() {
      return "(" + x + ", " + y + ", " + z + ")";
    }
  }

  /** A visitor of points as X,Y,Z values. */
  public interface PointVisitor {
    /** Receives a point as its X,Y,Z values. Returns true to keep visiting, or false to stop. */
    boolean accept(double x, double y, double z);
  }

  /** A visitor of indexed points, as an index and X,Y,Z values. */
  public interface PointOffsetVisitor {
    /**
     * Receives an indexed point as the index and X,Y,Z values at that index. Returns true to keep
     * visiting, or false to stop.
     */
    boolean accept(int index, double x, double y, double z);
  }

  /**
   * MutableS2PointList is a list of MutableS2Points, stored as a single array of doubles. It does
   * not initialize its contents, and reuses the underlying array unless it must be grown to
   * increase capacity. The owner is responsible for initializing all points before use.
   */
  public static class MutableS2PointList implements PullList<MutableS2Point> {
    /** The default initial capacity, as the number of points. */
    static final int DEFAULT_CAPACITY = 16;

    /**
     * The points, stored as consecutive X,Y,Z triples. The point with index 'i' has X,Y,Z values at
     * indices 3*i, 3*i+1, and 3*i+2. The array will often be larger than required for the current
     * size.
     */
    private double[] coordinates;

    /** The current number of points stored in the list. */
    private int size;

    /** The current capacity of the list. Always equal to coordinates.length / 3. */
    private int capacity;

    /** Creates a new, empty MutableS2PointList with a default initial capacity. */
    public MutableS2PointList() {
      this(new double[DEFAULT_CAPACITY * 3], 0);
    }

    /**
     * Constructs a new MutableS2PointList with a given initial size and capacity. All the points
     * will be S2.ORIGIN, i.e. all zeros.
     */
    public MutableS2PointList(int size) {
      this(new double[size * 3], size);
    }

    /**
     * Copy constructor. The new list will have the same size as the 'other' list, but capacity
     * exactly equal to the size.
     */
    public MutableS2PointList(MutableS2PointList other) {
      this(Arrays.copyOf(other.coordinates, other.size() * 3), other.size());
    }

    /**
     * Constructs a MutableS2PointList taking ownership of the given coordinates array and with the
     * given size. The given coordinates array must not be modified by the caller afterwards, so
     * this is private.
     */
    private MutableS2PointList(double[] coordinates, int size) {
      assert coordinates.length % 3 == 0;
      assert coordinates.length >= size * 3;
      this.coordinates = coordinates;
      this.size = size;
      this.capacity = coordinates.length / 3;
    }

    /** Creates a new, empty MutableS2PointList with a given capacity. */
    public static MutableS2PointList ofCapacity(int capacity) {
      MutableS2PointList list = new MutableS2PointList();
      list.ensureCapacity(capacity);
      return list;
    }

    /** Creates a new MutableS2PointList of capacity and size both equal to 2. */
    public static MutableS2PointList pair() {
      return new MutableS2PointList(2);
    }

    @CanIgnoreReturnValue
    @Override // Defined by PullCollection
    public boolean ensureCapacity(int requiredCapacity) {
      if (requiredCapacity <= capacity) {
        return false;
      }
      capacity = Integer.highestOneBit(requiredCapacity);
      capacity = (capacity == requiredCapacity) ? capacity : capacity << 1;
      coordinates = Arrays.copyOf(coordinates, capacity * 3);
      return true;
    }

    @Override // Defined by PullList
    public void enlarge(int newSize) {
      ensureCapacity(newSize);
      if (newSize > size) {
        size = newSize;
      }
    }

    @Override // Defined by PullList. Does not reduce capacity.
    public void truncate(int newSize) {
      if (newSize < size) {
        size = newSize;
      }
    }

    @Override // Defined by PullCollection
    public int size() {
      return size;
    }

    /**
     * Returns the current capacity of the list, i.e. how many points it can hold before
     * reallocating the internal array will be required.
     */
    public int capacity() {
      return capacity;
    }

    @Override // Defined by PullCollection
    public MutableS2Point newElement() {
      // TODO(torrey): Implement destroyElement and keep a pool of these for reuse.
      return new MutableS2PointImpl();
    }

    @Override // Defined by PullCollection
    public void clear() {
      size = 0;
    }

    /**
     * Sorts this MutableS2PointList using the standard lexicographical ordering of points, as
     * defined by {@link MutableS2Point#compareTo(MutableS2Point)}.
     */
    public void sort() {
      sort(Comparator.naturalOrder());
    }

    /** Throws IndexOutOfBoundsException if the index is out of range. */
    private void rangeCheck(int index) {
      if (index >= size || index < 0) {
        throw new IndexOutOfBoundsException("index " + index + " out of bounds for size " + size);
      }
    }

    /** Returns a new S2Point with the same coordinates as the MutableS2Point at the given index. */
    public S2Point getImmutable(int i) {
      rangeCheck(i);
      return new S2Point(coordinates[i * 3], coordinates[i * 3 + 1], coordinates[i * 3 + 2]);
    }

    /**
     * Returns true if the MutableS2Point at the given index in this list is equal to the given
     * S2Point.
     */
    public boolean isEqualTo(int index, S2Point other) {
      rangeCheck(index);
      int i = index * 3;
      return (coordinates[i + 0] == other.getX()
          && coordinates[i + 1] == other.getY()
          && coordinates[i + 2] == other.getZ());
    }

    /**
     * Returns true if the MutableS2Point at the given index in this list is equal to the given
     * MutableS2Point.
     */
    public boolean isEqualTo(int index, MutableS2Point other) {
      rangeCheck(index);
      int i = index * 3;
      return (coordinates[i + 0] == other.getX()
          && coordinates[i + 1] == other.getY()
          && coordinates[i + 2] == other.getZ());
    }

    /** Returns true if the MutableS2Points in this list at 'aIndex' and 'bIndex' are equal. */
    public boolean isEqualTo(int aIndex, int bIndex) {
      int ia = aIndex * 3;
      int ib = bIndex * 3;
      return (coordinates[ia + 0] == coordinates[ib + 0]
          && coordinates[ia + 1] == coordinates[ib + 1]
          && coordinates[ia + 2] == coordinates[ib + 2]);
    }

    @Override // Defined by PullList
    public void get(int index, MutableS2Point value) {
      rangeCheck(index);
      int i = index * 3;
      value.set(coordinates[i + 0], coordinates[i + 1], coordinates[i + 2]);
    }

    /**
     * Returns a view of the MutableS2Point at the given index in this list. The contents may be
     * changed via the returned MutableS2Point or by other means.
     *
     * <p>NOTE: The given index is range checked when get() is first called, but if this
     * MutableS2PointList is cleared or truncated to make size <= index, the returned MutableS2Point
     * then refers to an entry that is beyond (size - 1). That is NOT currently checked, and
     * operations on that MutableS2Point will not obviously fail, because the underlying array is
     * never reduced in size.
     */
    public MutableS2Point get(int index) {
      rangeCheck(index);
      return new MutableS2Point() {
        private final int i = index * 3;

        @Override
        public double getX() {
          return coordinates[i];
        }

        @Override
        public double getY() {
          return coordinates[i + 1];
        }

        @Override
        public double getZ() {
          return coordinates[i + 2];
        }

        @Override
        public void set(S2Point other) {
          coordinates[i + 0] = other.getX();
          coordinates[i + 1] = other.getY();
          coordinates[i + 2] = other.getZ();
        }

        @Override
        public void set(double x, double y, double z) {
          coordinates[i + 0] = x;
          coordinates[i + 1] = y;
          coordinates[i + 2] = z;
        }

        @Override
        public String toString() {
          return "("
              + coordinates[i + 0]
              + ", "
              + coordinates[i + 1]
              + ", "
              + coordinates[i + 2]
              + ")";
        }
      };
    }

    @Override // Defined by PullList
    public void set(int index, MutableS2Point point) {
      rangeCheck(index);
      int i = index * 3;
      coordinates[i + 0] = point.getX();
      coordinates[i + 1] = point.getY();
      coordinates[i + 2] = point.getZ();
    }

    /**
     * Sets the coordinates of the MutableS2Point at the given index in this list to equal the
     * coordinates of the given S2Point.
     */
    public void set(int index, S2Point point) {
      rangeCheck(index);
      int i = index * 3;
      coordinates[i + 0] = point.getX();
      coordinates[i + 1] = point.getY();
      coordinates[i + 2] = point.getZ();
    }

    /** Copies the point at indexA in this list to indexB, without bounds checks. */
    @Override // PullList's default implementation, for possibly greater efficiency
    public void copy(int indexA, int indexB) {
      int ia = indexA * 3;
      int ib = indexB * 3;
      coordinates[ib + 0] = coordinates[ia + 0];
      coordinates[ib + 1] = coordinates[ia + 1];
      coordinates[ib + 2] = coordinates[ia + 2];
    }

    @Override // Defined by PullCollection
    public void add(MutableS2Point value) {
      add(value.getX(), value.getY(), value.getZ());
    }

    /** Adds a copy of the given S2Point to the end of the list. */
    public void add(S2Point value) {
      add(value.x, value.y, value.z);
    }

    /** Adds a point with the given values to the list, increasing the size by one. */
    public void add(double x, double y, double z) {
      ensureCapacity(size + 1);
      int i = size * 3;
      coordinates[i + 0] = x;
      coordinates[i + 1] = y;
      coordinates[i + 2] = z;
      size++;
    }

    /**
     * Provides an {@code List<S2Point>} view of this MutableS2PointList. Modifications to the
     * underlying MutableS2PointList modify the returned List, and vice versa.
     */
    public List<S2Point> asPointList() {
      return new AbstractList<>() {
        @Override
        public S2Point get(int index) {
          return getImmutable(index);
        }

        @Override
        public S2Point set(int index, S2Point value) {
          S2Point returnValue = get(index);
          MutableS2PointList.this.set(index, value);
          return returnValue;
        }

        @Override
        public int size() {
          return MutableS2PointList.this.size();
        }
      };
    }

    /** Sends each point as X,Y,Z values to the given 'action'. */
    public boolean forEach(PointVisitor action) {
      for (int index = 0; index < size; index++) {
        int i = index * 3;
        if (!action.accept(coordinates[i + 0], coordinates[i + 1], coordinates[i + 2])) {
          return false;
        }
      }
      return true;
    }

    /** Sends each point as its index and X,Y,Z values to the given 'action'. */
    public boolean forEach(PointOffsetVisitor action) {
      for (int index = 0; index < size; index++) {
        int i = index * 3;
        if (!action.accept(index, coordinates[i + 0], coordinates[i + 1], coordinates[i + 2])) {
          return false;
        }
      }
      return true;
    }

    /** Swaps the points at indexA and indexB in this list. Does not bounds check vs. size. */
    @Override // PullList's default implementation, for possibly greater efficiency
    public void swap(int indexA, int indexB) {
      int ia = indexA * 3;
      int ib = indexB * 3;
      double tx = coordinates[ia + 0];
      double ty = coordinates[ia + 1];
      double tz = coordinates[ia + 2];
      coordinates[ia + 0] = coordinates[ib + 0];
      coordinates[ia + 1] = coordinates[ib + 1];
      coordinates[ia + 2] = coordinates[ib + 2];
      coordinates[ib + 0] = tx;
      coordinates[ib + 1] = ty;
      coordinates[ib + 2] = tz;
    }

    /** Compares the points at the given indices in this list. */
    public int compare(int leftIndex, int rightIndex) {
      int i = leftIndex * 3;
      int j = rightIndex * 3;
      int cmp = Double.compare(coordinates[i + 0], coordinates[j]);
      if (cmp != 0) {
        return cmp;
      }
      cmp = Double.compare(coordinates[i + 1], coordinates[j + 1]);
      if (cmp != 0) {
        return cmp;
      }
      return Double.compare(coordinates[i + 2], coordinates[j + 2]);
    }
  }
}
