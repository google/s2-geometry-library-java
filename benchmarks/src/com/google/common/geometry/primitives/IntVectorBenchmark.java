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

import static java.util.Comparator.naturalOrder;
import static java.util.concurrent.TimeUnit.MICROSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.geometry.benchmarks.S2BenchmarkBaseState;
import com.google.common.geometry.primitives.Ints.IntComparator;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.openjdk.jmh.annotations.Benchmark;
import org.openjdk.jmh.annotations.BenchmarkMode;
import org.openjdk.jmh.annotations.Level;
import org.openjdk.jmh.annotations.Measurement;
import org.openjdk.jmh.annotations.Mode;
import org.openjdk.jmh.annotations.OutputTimeUnit;
import org.openjdk.jmh.annotations.Param;
import org.openjdk.jmh.annotations.Scope;
import org.openjdk.jmh.annotations.Setup;
import org.openjdk.jmh.annotations.State;
import org.openjdk.jmh.annotations.Warmup;

/** Benchmarks for IntVector. */
public class IntVectorBenchmark {
  private IntVectorBenchmark() {}

  private static final IntComparator LEXICOGRAPHIC_COMPARATOR =
      (a, b) -> Integer.toString(a).compareTo(Integer.toString(b));

  /**
   * Benchmark state for comparing performance of IntVector to an ArrayList as well as a raw array.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
  public static class ArrayVsVectorVsList extends S2BenchmarkBaseState {
    static final int ITERS = 1000000;

    @Param({"256", "4096", "16384"})
    int numElements;

    // An array of random integers, each less than or equal to numElements.
    protected static final int NUM_LOCATIONS = 64;
    protected int[] locations = new int[NUM_LOCATIONS];

    // A fixed-size array of random small (< 256) integers.
    protected static final int NUM_RANDOMS = 64;
    protected int[] randoms = new int[NUM_RANDOMS];

    // The three test objects to benchmark.
    protected int[] array;
    protected IntVector vector;
    protected ArrayList<Integer> list;

    /** Reset the containers for each iteration. */
    @Setup(Level.Iteration)
    @Override
    public void setup() {
      super.setup();
      for (int i = 0; i < NUM_LOCATIONS; i++) {
        locations[i] = data.nextInt(numElements);
      }
      for (int i = 0; i < NUM_RANDOMS; i++) {
        randoms[i] = data.nextInt(256);
      }

      // Benchmarks begin with a vector, array, and list, each of size numElements, all zeros.
      vector = IntVector.ofSize(numElements);
      array = new int[numElements];
      list = new ArrayList<>(numElements);
      list.addAll(Collections.nCopies(numElements, 0));
    }

    /** Sets numElements elements sequentially and then sums them up sequentially. */
    @Benchmark
    public long linearFillAndReadArray() {
      for (int i = 0; i < numElements; i++) {
        array[i] = randoms[i % NUM_RANDOMS];
      }
      long sum = 0;
      for (int i = 0; i < numElements; i++) {
        sum += array[i];
      }
      return sum;
    }

    /** Sets numElements elements sequentially and then sums them up sequentially. */
    @Benchmark
    public long linearFillAndReadVector() {
      vector.clear();
      for (int i = 0; i < numElements; i++) {
        vector.push(randoms[i % NUM_RANDOMS]);
      }
      long sum = 0;
      for (int i = 0; i < numElements; i++) {
        sum += vector.get(i);
      }
      return sum;
    }

    /** Sets numElements elements sequentially and then sums them up sequentially. */
    @Benchmark
    public long linearFillAndReadList() {
      list.clear();
      for (int i = 0; i < numElements; i++) {
        list.add(randoms[i % NUM_RANDOMS]);
      }
      long sum = 0;
      for (int i = 0; i < numElements; i++) {
        sum += list.get(i);
      }
      return sum;
    }

    /** Construct an empty IntVector, grow it by numElements. */
    @Benchmark
    public int constructAndGrowVector() {
      IntVector empty = IntVector.ofSize(0);
      for (int i = 0; i < numElements; i++) {
        empty.add(randoms[i % NUM_RANDOMS]);
      }
      return empty.get(locations[numElements % NUM_LOCATIONS]);
    }

    /** Construct an empty List, grow it by numElements. */
    @Benchmark
    public int constructAndGrowList() {
      ArrayList<Integer> empty = new ArrayList<>();
      for (int i = 0; i < numElements; i++) {
        empty.add(randoms[i % NUM_RANDOMS]);
      }
      return empty.get(locations[numElements % NUM_LOCATIONS]);
    }

    /** Reads and updates ITERS elements at random locations. */
    @Benchmark
    public int randomUpdateArray() {
      int loc = 0;
      for (int i = 0; i < ITERS; i++) {
        loc = locations[i % NUM_LOCATIONS];
        array[loc] = array[loc] + randoms[i % NUM_RANDOMS];
      }
      return array[loc];
    }

    /** Reads and updates ITERS elements at random locations. */
    @Benchmark
    public int randomUpdateVector() {
      int loc = 0;
      for (int i = 0; i < ITERS; i++) {
        loc = locations[i % NUM_LOCATIONS];
        vector.set(loc, vector.get(loc) + randoms[i % NUM_RANDOMS]);
      }
      return vector.get(loc);
    }

    /** Reads and updates ITERS elements at random locations. */
    @Benchmark
    public int randomUpdateVectorAsList() {
      List<Integer> asList = vector.asList();
      int loc = 0;
      for (int i = 0; i < ITERS; i++) {
        loc = locations[i % NUM_LOCATIONS];
        asList.set(loc, asList.get(loc) + randoms[i % NUM_RANDOMS]);
      }
      return vector.get(loc);
    }

    /** Reads and updates ITERS elements at random locations. */
    @Benchmark
    public int randomUpdateList() {
      int loc = 0;
      for (int i = 0; i < ITERS; i++) {
        loc = locations[i % NUM_LOCATIONS];
        list.set(loc, list.get(loc) + randoms[i % NUM_RANDOMS]);
      }
      return list.get(loc);
    }
  }

  /**
   * Benchmark state for comparing sorting performance of IntSorter on an IntVector vs. standard
   * Java Arrays.sort and Collections.sort.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(MICROSECONDS)
  @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
  public static class SortArrayVsVectorVsList extends S2BenchmarkBaseState {
    @Param({"256", "4096", "16384"})
    int numElements;

    // The test objects to benchmark.
    protected int[] array;
    protected IntVector vector;
    protected ArrayList<Integer> list;

    // Wrap the given array in an AbstractList for Collections.sort.
    protected static AbstractList<Integer> getWrapper(final int[] intArray) {
      return new AbstractList<>() {
        @Override
        public Integer get(int index) {
          return intArray[index];
        }

        @Override
        public Integer set(int index, Integer value) {
          int p = intArray[index];
          intArray[index] = value;
          return p;
        }

        @Override
        public int size() {
          return intArray.length;
        }
      };
    }

    /** Randomize the contents for each trial. */
    @Setup(Level.Trial)
    @Override
    public void setup() {
      super.setup();

      // Each benchmark trial begins with a vector, array, and list containing the same random ints.
      vector = IntVector.ofSize(numElements);
      array = new int[numElements];
      list = new ArrayList<>(numElements);

      for (int i = 0; i < numElements; i++) {
        int r = data.nextInt();
        vector.set(i, r);
        array[i] = r;
        list.add(r);
      }
    }

    // ARRAY

    // Sorts a raw array using {@link Arrays#sort(int[])}.
    @Benchmark
    public int naturalSortRawArrayWithArraysSort() {
      Arrays.sort(array);
      return array[0] + array[numElements - 1];
    }
    // Sorts a raw array wrapped by an AbstractList using {@link List#sort()}.
    @Benchmark
    public int naturalSortWrappedArrayWithListSort() {
      AbstractList<Integer> wrapper = getWrapper(array);
      wrapper.sort(naturalOrder());
      return wrapper.get(0) + wrapper.get(numElements - 1);
    }
    // Sorts a raw array wrapped by an AbstractList using lexicographic ordering.
    @Benchmark
    public int lexicographicSortWrappedArrayWithListSort() {
      AbstractList<Integer> wrapper = getWrapper(array);
      wrapper.sort(LEXICOGRAPHIC_COMPARATOR);
      return wrapper.get(0) + wrapper.get(numElements - 1);
    }

    // INT VECTOR

    /**
     * IntVector.sort() uses {@link Arrays#sort()}. About 100x faster than IntSorter.
     */
    @Benchmark
    public int naturalSortVectorWithArraysSort() {
      vector.sort();
      return vector.get(0) + vector.get(numElements - 1);
    }

    /** Sorts an IntVector with {@link IntSorter#sort(MutableIntList)}. */
    @Benchmark
    public int naturalSortVectorWithIntSorter() {
      IntSorter.sort(vector);
      return vector.get(0) + vector.get(numElements - 1);
    }

    /** Sorts an IntVector as an AbstractList using {@link List#sort()}. */
    @Benchmark
    public int naturalSortWrappedVectorWithListSort() {
      vector.asList().sort(naturalOrder());
      return vector.get(0) + vector.get(numElements - 1);
    }

    /**
     * Sorts an IntVector with {@link IntSorter#sort(IntComparator, IntList)} and a
     * LEXICOGRAPHIC_COMPARATOR.
     */
    @Benchmark
    public int lexicographicSortVectorWithIntSorter() {
      IntSorter.sort(LEXICOGRAPHIC_COMPARATOR, vector);
      return vector.get(0) + vector.get(numElements - 1);
    }

    /**
     * Sorts an IntVector as an AbstractList using {@link List#sort()} and a
     * LEXICOGRAPHIC_COMPARATOR.
     */
    @Benchmark
    public int lexicographicSortWrappedVectorWithListSort() {
      vector.asList().sort(LEXICOGRAPHIC_COMPARATOR);
      return vector.get(0) + vector.get(numElements - 1);
    }

    // ARRAY LIST

    /** Sorts an ArrayList with {@link ArrayList#sort()}. */
    @Benchmark
    public int naturalSortList() {
      list.sort(naturalOrder());
      return list.get(0) + list.get(numElements - 1);
    }

    /** Sorts an ArrayList lexicograpically with {@link ArrayList#sort()}. */
    @Benchmark
    public int lexicographicSortList() {
      list.sort(LEXICOGRAPHIC_COMPARATOR);
      return list.get(0) + list.get(numElements - 1);
    }
  }
}
