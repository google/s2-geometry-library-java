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

import static java.lang.Math.min;

/**
 * An implementation of QuickSort for abstract collections implementing SortableCollection, which
 * requires only 'less' and 'swap'. Sorter is not instantiated, but provides a static sort() method.
 */
public final class Sorter {
  private Sorter() {}

  /** Any implementation of SortableCollection can be sorted by Sorter. */
  public interface SortableCollection {
    /** True if the element at leftIndex should be ordered before the element at rightIndex. */
    boolean less(int leftIndex, int rightIndex);

    /** Exchanges the elements at the given indices. */
    void swap(int leftIndex, int rightIndex);

    /** Returns the number of elements in this sortable. */
    int size();

    /** Sorts this collection. */
    default void sort() {
      if (size() > 1) {
        Sorter.sort(this, 0, size() - 1);
      }
    }

    /** Eliminates all elements at index 'start' and later. */
    void truncate(int start);

    /** Removes adjacent duplicates from this collection. */
    default void unique() {
      int size = size();
      if (size <= 1) {
        return;
      }
      int dst = 0;
      for (int i = 1; i < size; i++) {
        if (less(dst, i)) {
          swap(++dst, i);
        }
      }
      truncate(dst + 1);
    }
  }

  /** Sorts the contents of "data". */
  public static void sort(SortableCollection data) {
    sort(data, 0, data.size() - 1);
  }

  /** Sorts the contents of data from index "left" to index "right" inclusive. */
  public static void sort(SortableCollection data, int left, int right) {
    while ((right - left) >= 8) {
      int pivot = pickPivot(data, left, right);
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
        while ((pb <= pc) && !data.less(pivot, pb)) {
          if (!data.less(pb, pivot)) {
            // data[pb] == pivot
            data.swap(pa, pb);
            pivot = pa++;
          }
          pb++;
        }
        while ((pb <= pc) && !data.less(pc, pivot)) {
          if (!data.less(pivot, pc)) {
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
        sort(data, l1, r1);
      }
    }
    insertionSort(data, left, right);
  }

  private static void insertionSort(SortableCollection data, int left, int right) {
    for (int i = left; i <= right; i++) {
      for (int j = i; j > left && data.less(j, j - 1); j--) {
        data.swap(j, j - 1);
      }
    }
  }

  private static int pickPivot(SortableCollection data, int left, int right) {
    if ((right - left + 1) > 100) {
      // Pick median of medians
      int d = (right - left) / 8;
      int a = median(data, left + 0 * d, left + 1 * d, left + 2 * d);
      int b = median(data, left + 3 * d, left + 4 * d, left + 5 * d);
      int c = median(data, left + 6 * d, left + 7 * d, left + 8 * d);
      return median(data, a, b, c);
    } else {
      return median(data, left, (left + right) / 2, right);
    }
  }

  private static int median(SortableCollection data, int a, int b, int c) {
    return (data.less(a, b)
        ? (data.less(b, c) ? b : (data.less(a, c) ? c : a))
        : (data.less(c, b) ? b : (data.less(c, a) ? c : a)));
  }

  private static void swapRange(SortableCollection data, int p1, int p2, int n) {
    for (int i = 0; i < n; i++) {
      data.swap(p1 + i, p2 + i);
    }
  }
}
