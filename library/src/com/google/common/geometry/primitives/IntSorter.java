
/*
 * Copyright 2022 Google Inc.
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

import static com.google.common.geometry.primitives.Ints.INTEGER_COMPARATOR;
import static java.lang.Math.min;

import com.google.common.geometry.primitives.Ints.IntComparator;
import com.google.common.geometry.primitives.Ints.MutableIntList;

/**
 * An implementation of QuickSort for collections of primitive ints.
 *
 * <p>Collections implementing {@link Ints.MutableIntList} may be sorted in an order defined by a
 * provided {@link Ints.IntComparator}, or for convenience, natural integer ordering. This
 * implementation is NOT a stable sort, i.e. elements comparing as equal according to the given
 * IntComparator may have their relative ordering changed by the sort.
 */
public final class IntSorter {
  private IntSorter() {}

  /** Sorts the contents of 'data' using the natural integer ordering. */
  public static void sort(MutableIntList data) {
    sort(INTEGER_COMPARATOR, data, 0, data.size() - 1);
  }

  /** Sorts the contents of data from index "left" to index "right" inclusive in natural order. */
  public static void sort(MutableIntList data, int left, int right) {
    sort(INTEGER_COMPARATOR, data, left, right);
  }

  /** Sorts the contents of 'data' using the given 'comparator'. */
  public static void sort(IntComparator comparator, MutableIntList data) {
    sort(comparator, data, 0, data.size() - 1);
  }

  /**
   * Sorts the contents of data from index "left" to index "right" inclusive using the provided
   * IntComparator.
   */
  public static void sort(IntComparator comparator, MutableIntList data, int left, int right) {
    while ((right - left) >= 8) {
      int pivot = pickPivot(comparator, data, left, right);
      // Partition data[left,right] into five regions:
      //
      //     +------+-----------+------------------+-----------+------+
      //     |  ==  |     <     |         ?        |     >     |  ==  |
      //     +------+-----------+------------------+-----------+------+
      //      |      |           |                |           |      |
      //      v      v           v                v           v      v
      //      left   pa          pb               pc          pd    right
      //
      int pa = left;
      int pb = left;
      int pc = right;
      int pd = right;
      while (true) {
        // Shrink middle region: "pb,pc"
        while ((pb <= pc) && !less(comparator, data, pivot, pb)) {
          if (!less(comparator, data, pb, pivot)) {
            // data[pb] == pivot
            data.swap(pa, pb);
            pivot = pa++;
          }
          pb++;
        }
        while ((pb <= pc) && !less(comparator, data, pc, pivot)) {
          if (!less(comparator, data, pivot, pc)) {
            // data[pc] == pivot
            data.swap(pc, pd);
            pivot = pd--;
          }
          pc--;
        }
        if (pb > pc) {
          break;
        }
        if (pb == pivot) {
          pivot = pc;
        } else if (pc == pivot) {
          pivot = pb;
        }
        data.swap(pb, pc);
        pb++;
        pc--;
      }
      // Shuffle "==" data from both ends to middle of array
      int s;
      s = min(pa - left, pb - pa);
      swapRange(data, left, pb - s, s);
      s = min(pd - pc, right - pd);
      swapRange(data, pb, right + 1 - s, s);
      // Get subranges left to sort (avoid sorting "==" ranges)
      int l1 = left;
      int r1 = left + (pb - pa) - 1;
      int l2 = right + 1 - (pd - pc);
      int r2 = right;
      // Make "left,right" be larger part and "l1,r1" be smaller part
      if ((r1 - l1) < (r2 - l2)) {
        left = l2;
        right = r2;
      } else {
        left = l1;
        right = r1;
        l1 = l2;
        r1 = r2;
      }
      // Recurse on smaller part, iterate on larger
      if (l1 < r1) {
        sort(comparator, data, l1, r1);
      }
    }
    insertionSort(comparator, data, left, right);
  }

  /** Returns true if data[leftIndex] < data[rightIndex] according to the comparator. */
  private static boolean less(
      IntComparator comparator, MutableIntList data, int leftIndex, int rightIndex) {
    return comparator.compare(data.get(leftIndex), data.get(rightIndex)) < 0;
  }

  private static void insertionSort(
      IntComparator comparator, MutableIntList data, int left, int right) {
    for (int i = left; i <= right; i++) {
      for (int j = i; j > left && less(comparator, data, j, j - 1); j--) {
        data.swap(j, j - 1);
      }
    }
  }

  private static int pickPivot(IntComparator comparator, MutableIntList data, int left, int right) {
    if ((right - left + 1) > 100) {
      // Pick median of medians
      int d = (right - left) / 8;
      int a = median(comparator, data, left + 0 * d, left + 1 * d, left + 2 * d);
      int b = median(comparator, data, left + 3 * d, left + 4 * d, left + 5 * d);
      int c = median(comparator, data, left + 6 * d, left + 7 * d, left + 8 * d);
      return median(comparator, data, a, b, c);
    } else {
      return median(comparator, data, left, (left + right) / 2, right);
    }
  }

  private static int median(IntComparator comparator, MutableIntList data, int a, int b, int c) {
    return (less(comparator, data, a, b)
        ? (less(comparator, data, b, c) ? b : (less(comparator, data, a, c) ? c : a))
        : (less(comparator, data, c, b) ? b : (less(comparator, data, c, a) ? c : a)));
  }

  private static void swapRange(MutableIntList data, int p1, int p2, int n) {
    for (int i = 0; i < n; i++) {
      data.swap(p1 + i, p2 + i);
    }
  }
}
