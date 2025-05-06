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

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.geometry.primitives.Ints.ImmutableIntSequence;
import com.google.common.geometry.primitives.Ints.IntBiConsumer;
import com.google.common.geometry.primitives.Ints.IntConsumer;
import com.google.common.geometry.primitives.Ints.IntSequence;
import com.google.common.geometry.primitives.Ints.IntSequenceHasher;
import com.google.common.geometry.primitives.Ints.MixingHasher;
import com.google.common.geometry.primitives.Ints.OfInt;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArraySet;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntSet;
import java.util.List;

/**
 * SequenceLexicon is a collection of immutable sequences of primitive ints (i.e., tuples of varying
 * length). It automatically eliminates duplicate sequences, and assigns "sequence ids" to the
 * resulting unique sequences. Sequence ids are sequentially increasing integers. Sequences may be
 * retrieved from the lexicon by id. See also {@link IdSetLexicon}.
 *
 * <p>SequenceLexicon stores sequences using arrays of primitive ints, which are more compact and
 * faster than an {@code ArrayList<Integer>}.
 */
public class SequenceLexicon {
  /**
   * All the elements of all the sequences are stored in this two-dimensional array. The first array
   * index is a block number, the second index is an offset into the block. Blocks are allocated as
   * required. When the number of blocks required exceeds the capacity of the data[] array, it is
   * reallocated.
   */
  private int[][] data;
  /** The total number of sequence elements stored in the data array. */
  private int numElements;
  /** The number of blocks allocated in the data array. */
  private int numBlocks;

  /**
   * Just like the elements of the sequences, the index of the first element of each sequence is
   * stored in a two dimensional array of blocks which grows as required.
   */
  private int[][] begins;
  /** The total number of sequences stored in the lexicon. */
  private int numSequences;
  /** The number of blocks allocated in the begins array. */
  private int numBeginsBlocks;

  /**
   * The default initial capacity for the number of sequences stored in the lexicon. Exceeding this
   * causes a reallocation in IntMultimap which is somewhat expensive for time, but the IntMultimap
   * requires a minimum of 6 x DEFAULT_NUM_SEQUENCES bytes for storage.
   */
  private static final int DEFAULT_NUM_SEQUENCES = 128; // TODO(torrey): Determine what a good
  // default size is for this IntMultimap. Consider adding another constructor for SequenceLexicon
  // that allows the expected number of sequences and other size parameters to be specified.

  /**
   * The initial size of the first dimension of the data array: the number of blocks that can be
   * allocated.
   */
  private static final int DEFAULT_NUM_BLOCKS = 32;

  /** How many sequence entries per block? Must be a power of two. */
  private static final int BLOCK_SIZE = 2048;
  /** Mask an index to get offset within a block. */
  private static final int BLOCK_MASK = BLOCK_SIZE - 1;
  /** Shifting an index right by this many bits provides the block number. */
  private static final int BLOCK_SHIFT = 31 - Integer.numberOfLeadingZeros(BLOCK_SIZE);

  /** Map from sequence hashcodes to ids of sequences in this Lexicon with that hash, if any. */
  private final Int2ObjectMap<IntSet> hashCodeToSequenceIds;

  /**
   * The IntSequenceHasher for this SequenceLexicon. Different implementations are allowed for
   * testing.
   */
  private final IntSequenceHasher hasher;

  /** For testing. If true, an exception will be thrown if an IntSequenceHasher collision occurs. */
  private final boolean throwOnCollision;

  /** Constructs a new SequenceLexicon. */
  public SequenceLexicon() {
    hashCodeToSequenceIds = new Int2ObjectOpenHashMap<>();
    hasher = new MixingHasher();
    throwOnCollision = false;
    clear();
  }

  /**
   * Copy constructor. Does a deep copy of the data, and takes a reference to the IntSequenceHasher.
   * (This is safe because {@link Ints.IntSequenceHasher#hashCode()} results should not depend on
   * persistent state.)
   */
  public SequenceLexicon(SequenceLexicon other) {
    this.data = new int[other.data.length][];
    for (int i = 0; i < other.numBlocks; i++) {
      data[i] = new int[other.data[i].length];
      System.arraycopy(other.data[i], 0, data[i], 0, other.data[i].length);
    }
    this.numElements = other.numElements;
    this.numBlocks = other.numBlocks;

    this.begins = new int[other.begins.length][];
    for (int i = 0; i < other.numBeginsBlocks; i++) {
      this.begins[i] = new int[other.begins[i].length];
      System.arraycopy(other.begins[i], 0, begins[i], 0, other.begins[i].length);
    }
    this.numSequences = other.numSequences;
    this.numBeginsBlocks = other.numBeginsBlocks;

    this.hashCodeToSequenceIds = new Int2ObjectOpenHashMap<>();
    for (Int2ObjectMap.Entry<IntSet> entry : other.hashCodeToSequenceIds.int2ObjectEntrySet()) {
      IntArraySet copy = new IntArraySet();
      copy.addAll(entry.getValue());
      hashCodeToSequenceIds.put(entry.getIntKey(), copy);
    }

    this.hasher = other.hasher;
    this.throwOnCollision = other.throwOnCollision;
  }

  /**
   * For testing only. Constructs a new SequenceLexicon that uses a provided IntSequenceHasher, and
   * optionally throws an exception if two different sequences hash to the same value.
   */
  @VisibleForTesting
  SequenceLexicon(IntSequenceHasher testHasher, boolean throwOnCollision) {
    hashCodeToSequenceIds = new Int2ObjectOpenHashMap<>(DEFAULT_NUM_SEQUENCES);
    hasher = testHasher;
    this.throwOnCollision = throwOnCollision;
    clear();
  }

  /** Clears all data from the lexicon. */
  public void clear() {
    data = new int[DEFAULT_NUM_BLOCKS][];
    numElements = 0;
    numBlocks = 0;

    begins = new int[DEFAULT_NUM_BLOCKS][];
    // Allocate the first block of the "begins" array to track where the first sequence will begin.
    begins[0] = new int[BLOCK_SIZE];
    begins[0][0] = 0; // Sequence #0 will start at element 0.
    numSequences = 0;
    numBeginsBlocks = 1;

    hashCodeToSequenceIds.clear();
  }

  /**
   * Returns true if this SequenceLexicon currently contains the same sequences with the same ids as
   * the given SequenceLexicon 'other'. As SequenceLexicon assigns sequence ids consecutively when
   * sequences are added, this means the sequences were added in the same order. Does not override
   * Object.equals() because SequenceLexicon is mutable.
   */
  public boolean isEqualTo(SequenceLexicon other) {
    if ((numElements != other.numElements) || (numSequences != other.numSequences)) {
      return false;
    }

    for (int i = 0; i < numSequences; i++) {
      if (!hasher.equals(sequence(i), other.sequence(i))) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns the sequence id of the given IntSequence 'newSequence'. If it is not already present in
   * the lexicon, it is added and assigned a new id. Ids are assigned sequentially starting from
   * zero.
   */
  @CanIgnoreReturnValue
  public int add(IntSequence newSequence) {
    // TODO(torrey): Add unit tests for adding a zero-length sequence.
    // Determine if this sequence is already in the lexicon by checking existing sequences with
    // the same hash. Normally, there is at most one.
    int hashCode = hasher.hashCode(newSequence);
    IntSet existingSequences = hashCodeToSequenceIds.get(hashCode);
    if (existingSequences != null) {
      IntIterator existingSequencesIter = existingSequences.iterator();
      while (existingSequencesIter.hasNext()) {
        int existingId = existingSequencesIter.nextInt();
        ImmutableIntSequence existingSequence = sequence(existingId);
        // Check that existingSequence and newSequence are actually the same. They usually are.
        if (hasher.equals(existingSequence, newSequence)) {
          return existingId;
        }
        // Different Sequences with the same hash! This is probably a test, but otherwise try the
        // next sequence.
        if (throwOnCollision) {
          throw new IllegalStateException(
              "IntSequenceHasher collision: adding new sequence:\n"
                  + newSequence.debugString()
                  + "Collided with existing sequence:\n"
                  + existingSequence.debugString());
        }
      }
    } else {
      existingSequences = new IntArraySet();
      hashCodeToSequenceIds.put(hashCode, existingSequences);
    }

    // The newSequence isn't already in the lexicon. Assign it the next id.
    int newId = numSequences;
    existingSequences.add(newId);

    // Store the elements of the sequence in the order provided, updating numElements as we go.
    // TODO(torrey): This could likely be faster if addElement was inlined here and unrolled to not
    // compute blockNum, check space remaining, etc. for every element.
    newSequence.forEach((IntConsumer) this::addElement);
    // Record where the *next* sequence's elements will begin, and update numSequences.
    addNextBegin(numElements);
    return newId;
  }

  /** Like add(IntSequence) but takes a {@code List<Integer>}. */
  @CanIgnoreReturnValue
  public int add(List<Integer> values) {
    return add(ImmutableIntSequence.of(values));
  }

  /**
   * Adds the given element to the end of the data array, adding another block if needed, and
   * reallocating to get space for more blocks if needed.
   */
  private void addElement(int element) {
    int blockNum = numElements >>> BLOCK_SHIFT;
    if (numBlocks == blockNum) {
      if (data.length == blockNum) {
        reallocateData(data.length + DEFAULT_NUM_BLOCKS);
      }
      data[numBlocks++] = new int[BLOCK_SIZE];
    }
    data[blockNum][numElements & BLOCK_MASK] = element;
    numElements++;
  }

  /**
   * Adds the given value "begin" to the 'begins' array, adding another block if needed, and
   * reallocating to get more space for blocks if needed.
   */
  private void addNextBegin(int index) {
    // We are always storing the beginning index of the next sequence, at (numSequences + 1), so
    // unlike addElement() above, numSequences is incremented at the beginning rather than the end.
    numSequences++;
    int blockNum = numSequences >>> BLOCK_SHIFT;
    if (blockNum == numBeginsBlocks) {
      if (begins.length == blockNum) {
        reallocateBegins(data.length + 32);
      }
      begins[numBeginsBlocks++] = new int[BLOCK_SIZE];
    }
    begins[blockNum][numSequences & BLOCK_MASK] = index;
  }

  /** Returns the number of value sequences in the lexicon. */
  public int size() {
    return numSequences;
  }

  /** Package private for testing: the block size. */
  @VisibleForTesting
  int blockSize() {
    return BLOCK_SIZE;
  }

  /** Package private for testing: the total number of blocks allocated for begins. */
  @VisibleForTesting
  int numBeginsBlocks() {
    return numBeginsBlocks;
  }

  /** Package private for testing: the length of the begins array first dimension. */
  @VisibleForTesting
  int beginsLength() {
    return begins.length;
  }

  /** Package private for testing: the total number of sequence elements stored. */
  @VisibleForTesting
  int numElements() {
    return numElements;
  }

  /** Package private for testing: the total number of blocks allocated for data. */
  @VisibleForTesting
  int numBlocks() {
    return numBlocks;
  }

  /** Package private for testing: the length of the data array first dimension. */
  @VisibleForTesting
  int dataLength() {
    return data.length;
  }

  /** Returns the index of the first element of the given sequence. */
  private int getBegin(int sequenceId) {
    return begins[sequenceId >> BLOCK_SHIFT][sequenceId & BLOCK_MASK];
  }

  /** Returns the sequence element at the given index. */
  private int getElement(int index) {
    return data[index >> BLOCK_SHIFT][index & BLOCK_MASK];
  }

  /** Increase the size of the first dimension of the data array to hold 'blocks' blocks. */
  private void reallocateData(int blocks) {
    int[][] newData = new int[blocks][];
    // Copy references to the existing blocks.
    System.arraycopy(data, 0, newData, 0, numBlocks);
    data = newData;
  }

  /** Increase the size of the first dimension of the begins array to hold 'blocks' blocks. */
  private void reallocateBegins(int blocks) {
    int[][] newBegins = new int[blocks][];
    // Copy references to the existing blocks.
    System.arraycopy(begins, 0, newBegins, 0, numBeginsBlocks);
    begins = newBegins;
  }

  /**
   * Returns a ImmutableIntSequence accessor for the value sequence with the given id, or throws
   * IllegalArgumentException if there is no sequence in the lexicon with that id.
   */
  public ImmutableIntSequence sequence(int sequenceId) {
    Preconditions.checkArgument(sequenceId < numSequences && sequenceId >= 0);

    return new ImmutableIntSequence() {
      // An empty sequence has begin == end.
      private final int begin = getBegin(sequenceId);
      private final int end = getBegin(sequenceId + 1);

      @Override
      public int size() {
        return end - begin;
      }

      @Override
      public void forEach(IntConsumer action) {
        for (int i = begin; i < end; i++) {
          action.accept(getElement(i));
        }
      }

      @Override
      public void forEach(IntBiConsumer action) {
        for (int i = begin; i < end; i++) {
          action.accept(i - begin, getElement(i));
        }
      }

      @Override
      public OfInt intIterator() {
        return new OfInt() {
          int position = begin;

          @Override
          public boolean hasNext() {
            return position < end;
          }

          @Override
          public int nextInt() {
            Preconditions.checkState(position < end);
            return getElement(position++);
          }

          @Override
          public void forEachRemaining(IntConsumer action) {
            while (position < end) {
              action.accept(getElement(position++));
            }
          }
        };
      }
    };
  }
}
