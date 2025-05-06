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

import com.google.common.geometry.primitives.Sorter.SortableCollection;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.util.AbstractList;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;
import java.util.stream.IntStream;

/**
 * Interfaces for working with collections of primitive ints.
 *
 * <p>These interfaces are only to be used internally by S2, and may be removed and replaced by
 * Guava or FastUtil or something else in the future.
 */
public final class Ints {
  private Ints() {}

  /** How many elements should {@link IntSequence#debugString()} include? */
  static final int DEBUG_STRING_MAX_ELEMENTS = 16;

  /** An Android-compatible replacement for java.util.function.IntConsumer */
  public interface IntConsumer {
    void accept(int element);
  }

  /** A consumer of pairs of integers. */
  public interface IntBiConsumer {
    void accept(int key, int value);
  }

  /** An Android-compatible replacement for java.util.PrimitiveIterator.OfInt */
  public static interface OfInt {
    int nextInt();
    boolean hasNext();
    default void forEachRemaining(IntConsumer action) {
      while (hasNext()) {
        action.accept(nextInt());
      }
    }
  }

  /** An Android-compatible replacement for java.util.function.IntPredicate */
  public static interface IntPredicate {
    boolean test(int value);
  }

  /**
   * IntIterators provides useful static implementations of the {@link OfInt} iterator interface.
   */
  public static final class IntIterators {
    private IntIterators() {}

    /** An 'OfInt' iterator that contains no elements. */
    public static final OfInt EMPTY =
        new OfInt() {
          @Override
          public boolean hasNext() {
            return false;
          }

          @Override
          public int nextInt() {
            throw new NoSuchElementException("nextInt() called on IntIterator.EMPTY");
          }

          @Override
          public void forEachRemaining(IntConsumer action) {}
        };

    /** Creates a new 'OfInt' iterator over the single given int value. */
    public static OfInt of(int value) {
      return new OfInt() {
        boolean hasNext = true;

        @Override
        public boolean hasNext() {
          return hasNext;
        }

        @Override
        public int nextInt() {
          if (hasNext) {
            hasNext = false;
            return value;
          }
          throw new NoSuchElementException();
        }

        @Override public void forEachRemaining(IntConsumer action) {
          if (hasNext) {
            hasNext = false;
            action.accept(value);
          }
        }
      };
    }

    /** Creates a new 'OfInt' iterator over the given array of ints. */
    public static OfInt of(int... ints) {
      return new OfInt() {
        int i = 0;

        @Override
        public boolean hasNext() {
          return i < ints.length;
        }

        @Override
        public int nextInt() {
          return ints[i++]; // throws NoSuchElementException past the end of the array
        }

        @Override
        public void forEachRemaining(IntConsumer action) {
          while (i < ints.length) {
            action.accept(i++);
          }
        }
      };
    }

    /** Creates a new 'OfInt' iterator over the given boxed {@code Iterable<Integer>}. */
    public static OfInt of(Iterable<Integer> ints) {
      return new OfInt() {
        Iterator<Integer> iter = ints.iterator();

        @Override
        public boolean hasNext() {
          return iter.hasNext();
        }

        @Override
        public int nextInt() {
          return iter.next();
        }

        @Override
        public void forEachRemaining(IntConsumer action) {
          while (iter.hasNext()) {
            action.accept(iter.next());
          }
        }
      };
    }
  }

  /**
   * Returns an IntConsumer that places ints into the given 'intArray' in order. The given array
   * must be large enough to hold all the elements.
   */
  public static IntConsumer arraySink(final int[] intArray) {
    return new IntConsumer() {
      int index = 0;

      @Override
      public void accept(int element) {
        intArray[index++] = element;
      }
    };
  }

  /**
   * IntSequence provides an ordered collection of primitive ints with known size. It is similar to
   * an {@code Iterable<Integer>} and could easily extend that interface, but specifically chooses
   * not to, so as to ensure that users do not accidentally use boxing iteration. Elements of a
   * Sequence may be iterated with {@link #intIterator()} or visited with {@link
   * #forEach(IntConsumer)}.
   *
   * <p>The Sequence interface does not provide operations that change the size or contents of the
   * Sequence, but implementations may be mutable: {@link MutableIntList} is one example. See also
   * {@link ImmutableIntSequence}. Sequence does not provide random access. {@link IntList} extends
   * Sequence to provide that capability.
   *
   * <p>As an interface, Sequence cannot provide implementations of equals(), hashcode(), and
   * toString() that override the implementation on Object. Immutable implementations should
   * implement equals() and hashcode() using MIXING_HASHER, as {@link ImmutableIntSequence} does.
   * Mutable implementations should probably not implement hashCode() or equals. All implementations
   * may implement toString() using debugString().
   */
  public interface IntSequence {
    /** Returns an 'OfInt' iterator over the ints in this IntSequence. */
    OfInt intIterator();

    /** Calls action.accept() for every int in this IntSequence in order. */
    void forEach(IntConsumer action);

    /**
     * Calls action.accept() for every index and value at that index in this IntSequence in order.
     */
    void forEach(IntBiConsumer action);

    /** Returns the number of elements in this IntSequence. */
    int size();

    /** Returns true if this IntSequence contains no elements. */
    default boolean isEmpty() {
      return size() == 0;
    }

    /** Returns true if this IntSequence contains the given value. This is a linear search. */
    default boolean contains(int value) {
      for (OfInt iter = intIterator(); iter.hasNext(); ) {
        if (iter.nextInt() == value) {
          return true;
        }
      }
      return false;
    }

    /** Returns an IntStream of this IntSequence. */
    default IntStream stream() {
      return IntStream.generate(intIterator()::nextInt).limit(size());
    }

    /** Returns a new int[] containing a copy of the contents of this IntSequence. */
    default int[] toArray() {
      final int[] result = new int[size()];
      forEach(arraySink(result));
      return result;
    }

    /** Creates a new IntSequence wrapping the given {@code List<Integer>}. */
    public static IntSequence of(List<Integer> ints) {
      return new IntSequence() {
        @Override
        public OfInt intIterator() {
          return IntIterators.of(ints);
        }

        @Override
        public void forEach(IntConsumer action) {
          for (int i : ints) {
            action.accept(i);
          }
        }

        @Override
        public void forEach(IntBiConsumer action) {
          for (int i = 0; i < ints.size(); i++) {
            action.accept(i, ints.get(i));
          }
        }

        @Override
        public int size() {
          return ints.size();
        }

        @Override
        public String toString() {
          return debugString();
        }
      };
    }

    /**
     * Returns a String representation of up to the first DEBUG_STRING_MAX_ELEMENTS elements of this
     * IntSequence. Mainly intended for use in unit tests.
     */
    default String debugString() {
      StringBuilder sb = new StringBuilder();
      sb.append("IntSequence of ").append(size()).append(" elements");

      if (isEmpty()) {
        sb.append(".");
        return sb.toString();
      }

      OfInt iter = intIterator();
      sb.append(": [").append(iter.nextInt());
      int count = 1;
      while (iter.hasNext() && count++ < DEBUG_STRING_MAX_ELEMENTS) {
        sb.append(", ").append(iter.nextInt());
      }
      if (size() > DEBUG_STRING_MAX_ELEMENTS) {
        sb.append("...");
      }
      sb.append("]");
      return sb.toString();
    }
  }

  /**
   * An abstract immutable IntSequence class, and some convenience implementations thereof.
   *
   * <p>ImmutableIntSequence extends IntSequence to provide equals() and hashCode(). IntSequences
   * are considered equal if they have the same elements in the same order. Similarly, hash codes of
   * IntSequences should only depend on the order and value of the values, not on the type or
   * implementation details of the IntSequence. The MixingHasher implementation of IntSequenceHasher
   * guarantees this.
   */
  public abstract static class ImmutableIntSequence implements IntSequence {
    @Override
    public boolean equals(Object obj) {
      if (this == obj) {
        return true;
      }
      if (obj instanceof IntSequence) {
        return MIXING_HASHER.equals(this, (IntSequence) obj);
      }
      return false;
    }

    @Override
    public int hashCode() {
      return MIXING_HASHER.hashCode(this);
    }

    @Override
    public String toString() {
      return debugString();
    }

    /** A constant ImmutableIntSequence that contains no elements. */
    public static final ImmutableIntSequence EMPTY =
        new ImmutableIntSequence() {
          @Override
          public OfInt intIterator() {
            return IntIterators.EMPTY;
          }

          @Override
          public void forEach(IntConsumer action) {}

          @Override
          public void forEach(IntBiConsumer action) {}

          @Override
          public int size() {
            return 0;
          }
        };

    /** Creates a new ImmutableIntSequence containing the single given value. */
    public static ImmutableIntSequence of(int value) {
      return new ImmutableIntSequence() {
        @Override
        public OfInt intIterator() {
          return IntIterators.of(value);
        }

        @Override
        public void forEach(IntConsumer action) {
          action.accept(value);
        }

        @Override
        public void forEach(IntBiConsumer action) {
          action.accept(0, value);
        }

        @Override
        public int size() {
          return 1;
        }
      };
    }

    /** Creates a new ImmutableIntSequence wrapper for the given ints. */
    public static ImmutableIntSequence of(int... ints) {
      return new ImmutableIntSequence() {
        @Override
        public OfInt intIterator() {
          return IntIterators.of(ints);
        }

        @Override
        public void forEach(IntConsumer action) {
          for (int i : ints) {
            action.accept(i);
          }
        }

        @Override
        public void forEach(IntBiConsumer action) {
          for (int i = 0; i < ints.length; i++) {
            action.accept(i, ints[i]);
          }
        }

        @Override
        public int size() {
          return ints.length;
        }
      };
    }

    /**
     * Creates a new ImmutableIntSequence wrapping the given {@code List<Integer>}.
     *
     * <p>Note that although the ImmutableIntSequence interface does not provide mutating methods,
     * the underlying List might be mutated separately; in general this is likely to cause
     * non-deterministic behavior and should be avoided.
     */
    public static ImmutableIntSequence of(List<Integer> ints) {
      return new ImmutableIntSequence() {
        @Override
        public OfInt intIterator() {
          return IntIterators.of(ints);
        }

        @Override
        public void forEach(IntConsumer action) {
          for (int i : ints) {
            action.accept(i);
          }
        }

        @Override
        public void forEach(IntBiConsumer action) {
          for (int i = 0; i < ints.size(); i++) {
            action.accept(i, ints.get(i));
          }
        }

        @Override
        public int size() {
          return ints.size();
        }
      };
    }

    /**
     * Creates a new ImmutableIntSequence wrapping the given IntSequence.
     *
     * <p>Note that although the ImmutableIntSequence interface does not provide mutating methods,
     * the underlying List might be mutated separately; in general this is likely to cause
     * non-deterministic behavior and should be avoided.
     */
    public static ImmutableIntSequence viewOf(IntSequence ints) {
      return new ImmutableIntSequence() {
        @Override
        public OfInt intIterator() {
          return ints.intIterator();
        }

        @Override
        public void forEach(IntConsumer action) {
          ints.forEach(action);
        }

        @Override
        public void forEach(IntBiConsumer action) {
          ints.forEach(action);
        }

        @Override
        public int size() {
          return ints.size();
        }
      };
    }
  }

  /**
   * Defines equals() and hashCode() for IntSequences.
   *
   * <p>The IntSequence interface does not provide mutation methods, but implementations may be
   * mutable. Use equals() and hashCode() on mutable IntSequences with caution. As with any mutable
   * collection, the results of a call to equals() or hashcode() will become invalid if the
   * IntSequence is modified.
   */
  public interface IntSequenceHasher {
    /** Returns a hashcode for the given IntSequence of integers. */
    int hashCode(IntSequence iterable);

    /** Returns true if and only if the two IntSequences are equal. */
    boolean equals(IntSequence a, IntSequence b);
  }

  /**
   * An IntSequenceHasher which defines equality based on the order and value of the integers in the
   * sequence. Hash code computation is using "fast mixing" based on util/hash/mix.h.
   */
  public static class MixingHasher implements IntSequenceHasher {
    // Internal state is a long, but truncated to an int for output.
    private static final long MUL = 0xdc3eb94af8ab4c93L;
    private long hash;

    public MixingHasher() {}

    private void reset() {
      hash = 84;
    }

    private void mix(long val) {
      hash *= MUL;
      hash = ((hash << 19) | (hash >>> (64 - 19))) + val; // unsigned right shift
    }

    @Override
    public int hashCode(IntSequence intIter) {
      reset();
      intIter.forEach((IntConsumer) this::mix);
      return (int) hash;
    }

    /**
     * IntSequences are considered equal if they produce the same values in the same order when
     * iterated over.
     */
    @Override
    public boolean equals(IntSequence a, IntSequence b) {
      if (a.size() != b.size()) {
        return false;
      }
      OfInt aIter = a.intIterator();
      OfInt bIter = b.intIterator();
      while (aIter.hasNext()) {
        if (!bIter.hasNext()) {
          return false;
        }
        if (aIter.nextInt() != bIter.nextInt()) {
          return false;
        }
      }
      return !bIter.hasNext();
    }
  }

  public static final MixingHasher MIXING_HASHER = new MixingHasher();

  /** Analogous to a {@code Comparator<Integer>}, but for primitive ints. */
  public interface IntComparator extends Comparator<Integer> {
    /**
     * Compares its two arguments for order. Returns a negative integer, zero, or a positive integer
     * when "a" is less than, equal to, or greater than "b", respectively.
     */
    int compare(int a, int b);

    @Override
    default int compare(Integer a, Integer b) {
      return compare(a.intValue(), b.intValue());
    }
  }

  /** Analogous to a {@code Comparator<Pair<Integer, Integer>>}, but without boxing. */
  public interface IntPairComparator {
    /**
     * Compares two pairs for order. Returns a negative integer, zero, or a positive integer when
     * the pair "(a0, a1)" is less than, equal to, or greater than "(b0, b1)", respectively.
     */
    int comparePair(int a0, int a1, int b0, int b1);
  }

  /** The normal comparator for sorting integers in order. */
  public static final IntComparator INTEGER_COMPARATOR = new IntComparator() {
    @Override
    public int compare(int a, int b) {
      return Integer.compare(a, b);
    }
  };

  /**
   * IntList provides the non-mutating methods of {@link MutableIntList}. It extends IntSequence,
   * and adds random access to values by index.
   */
  public interface IntList extends IntSequence {
    /**
     * Returns the element at the specified position in this list.
     *
     * @throws IndexOutOfBoundsException if the index is negative, or greater than or equal to
     *     size().
     */
    int get(int index);

    /**
     * Checks that the given 'fromIndex' and 'toIndex' are valid indexes for this IntList, both
     * considered inclusively, and that toIndex >= fromIndex.
     *
     * @throws IndexOutOfBoundsException if either index is negative, or greater than or equal to
     *     size().
     */
    default void rangeCheck(int fromIndex, int toIndex) {
      if (fromIndex > toIndex) {
        throw new IllegalArgumentException(
            "fromIndex " + fromIndex + " greater than toIndex " + toIndex);
      }
      if (fromIndex < 0) {
        throw new ArrayIndexOutOfBoundsException("fromIndex " + fromIndex + " out of bounds.");
      }
      if (toIndex >= size()) {
        throw new ArrayIndexOutOfBoundsException(
            "toIndex " + toIndex + " out of bounds on IntList of size " + size());
      }
    }

    /**
     * Returns true if this IntList currently contains the same ints in the same order as the given
     * IntList 'other'. Does not override Object.equals() because most implementations of IntList
     * are mutable.
     */
    default boolean isEqualTo(IntList other) {
      if (size() != other.size()) {
        return false;
      }
      for (int i = 0; i < size(); i++) {
        if (get(i) != other.get(i)) {
          return false;
        }
      }
      return true;
    }
  }

  // TODO(user): Consider removing MutableIntList as an interface that is separate from the
  // IntVector implementation, if IntVector is and will be the only implementation.
  /**
   * MutableIntList is analogous to {@code List<Integer>} for primitive ints, except that nulls are
   * not permitted.
   *
   * <p>MutableIntList extends IntList, which extends IntSequence. It inherits iteration methods and
   * size(). It adds methods inspired by java.util.List and the C++ standard library's {@code
   * vector<int>}. Mutable IntList capabilities include:
   *
   * <ul>
   *   <li>Methods to grow or truncate to a specified size, or clear the contents.
   *   <li>Methods for append values from various sources to the end of an MutableIntList.
   *   <li>Methods to use an MutableIntList as a stack: The top of the stack is at the end of the
   *       list, so pop() returns the value at (size() - 1).
   *   <li>Methods to initialize or fill an MutableIntList with a constant value or consecutive
   *       values.
   *   <li>An MutableIntList may be sorted using IntSorter and an IntComparator.
   * </ul>
   *
   * MutableIntList is a general purpose container, but was particularly intended to be a
   * replacement for {@code std::vector<int>} when porting C++ to Java. High performance was a
   * priority. See {@link IntVector}, which implements IntList, for performance details and
   * benchmark results.
   */
  public interface MutableIntList extends IntList, SortableCollection {
    /** Appends the specified element to the end of this list, increasing the size by one. */
    void add(int e);

    /**
     * Appends all of the elements in the given IntSequence of ints to the end of this list,
     * increasing the size by the number of elements added.
     */
    void addAll(IntSequence ints);

    /**
     * Appends all of the Integers in the given List to the end of this list, increasing the size by
     * the number of elements added.
     */
    void addAll(List<Integer> values);

    /** Removes all of the elements from this MutableIntList, making the size zero. */
    void clear();

    /** Sets all elements of the MutableIntList to the given value. The size is not changed. */
    void fill(int value);

    /**
     * Initialize the list with consecutive integers beginning with zero. The size is not changed.
     */
    void fillConsecutive();

    /**
     * Removes all elements from the list from the given 'fromIndex' (inclusive) to the end, if the
     * current size is larger than fromIndex. Does nothing if the current size is smaller than
     * fromIndex. The size of the list will be the smaller of the current size or 'fromIndex' after
     * this operation.
     */
    @Override
    void truncate(int fromIndex);

    /**
     * Increases the size of the list to the given 'size', if the list is not already at least that
     * size. Truncates the list to the given 'size', if the list is larger. Does nothing if the
     * given size is equal to the current size.
     */
    void resize(int size);

    /**
     * Changes the element at the specified position in this list to the specified value.
     *
     * @throws IndexOutOfBoundsException if the index is negative, or greater than or equal to
     *     size().
     */
    void set(int index, int value);

    /** Swaps the value of the data elements at the two indices. */
    @Override
    void swap(int index1, int index2);

    /**
     * Returns a view of the contents of this IntList from {@code fromIndex} (inclusive) to {@code
     * toIndex} (exclusive).
     */
    IntSequence subList(int fromIndex, int toIndex);

    /** Increments the value at 'index' and returns the result. */
    @CanIgnoreReturnValue
    int incrementAt(int index);

    /** Increments the value at 'index' by 'amount' and returns the result. */
    @CanIgnoreReturnValue
    int incrementAt(int index, int amount);

    /** Decrements the value at 'index' and returns the result. */
    @CanIgnoreReturnValue
    int decrementAt(int index);

    /** Decrements the value at 'index' by 'amount' and returns the result. */
    @CanIgnoreReturnValue
    int decrementAt(int index, int amount);

    // Stack methods:

    /** Adds the given value to the top of the stack, growing it by one. */
    void push(int value);
    /** Removes and returns the value on the top of the stack, shrinking it by one. */
    int pop();
    /** Returns the value on the top of the stack. Doesn't modify the stack. */
    int peek();

    /**
     * Provides an {@code AbstractList<Integer>} view of this MutableIntList: modifications to the
     * MutableIntList modify the AbstractList, and vice versa.
     */
    default AbstractList<Integer> asList() {
      final MutableIntList ints = this;
      return new AbstractList<>() {
        @Override
        public Integer get(int index) {
          return ints.get(index);
        }

        @Override
        public Integer set(int index, Integer value) {
          int p = ints.get(index);
          ints.set(index, value);
          return p;
        }

        @Override
        public int size() {
          return ints.size();
        }
      };
    }
  }
}
