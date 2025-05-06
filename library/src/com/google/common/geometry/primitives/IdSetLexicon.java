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

import com.google.common.base.Preconditions;
import com.google.common.geometry.primitives.Ints.ImmutableIntSequence;
import com.google.common.geometry.primitives.Ints.IntBiConsumer;
import com.google.common.geometry.primitives.Ints.IntConsumer;
import com.google.common.geometry.primitives.Ints.IntSequence;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.Collection;
import java.util.NoSuchElementException;

/**
 * IdSetLexicon is a class for compactly representing immutable sets of non-negative integers such
 * as array indices ("id sets"). It is especially suitable when either (1) there are many duplicate
 * sets, or (2) there are many singleton or empty sets. See also {@link SequenceLexicon}.
 *
 * <p>Each distinct id set is mapped to a 32-bit integer. Empty and singleton sets take up no
 * additional space whatsoever; the set itself is represented by the unique id assigned to the set
 * and the lexicon doesn't actually store anything for such sets.
 *
 * <p>Duplicate sets are automatically eliminated. Note also that id sets are referred to using
 * 32-bit integers rather than object references.
 *
 * IdSetLexicon is similar to SequenceLexicon, except:
 *
 * <ol>
 *   <li>Empty and singleton sets are represented implicitly; they use no space.
 *   <li>Sets are represented rather than sequences; values are reordered to be in sorted order and
 *       duplicate values are removed.
 *   <li>The values must be 32-bit non-negative integers (only).
 * </ol>
 */
public class IdSetLexicon {

  /** An immutable set of non-negative ints, in sorted order. */
  public abstract static class IdSet extends ImmutableIntSequence {
    /** Returns the first (smallest) element in this set. */
    public abstract int first();
  }

  /**
   * An implementation of IdSet that is not actually stored in an IdSetSetLexicon. This is mainly
   * intended for use in testing. For example, a unit test can build an expected TestIdSet and check
   * equality vs. an actual IdSet returned from an IdSetLexicon.
   *
   * <p>TestIdSet automatically sorts the elements provided to its constructors and ensures they are
   * unique, as required by the IdSet contract. The equals() implementation is inherited from
   * {@link ImmutableIntSequence}.
   */
  public static class TestIdSet extends IdSet {
    private final IntVector elements;

    /** Constructs an empty TestIdSet. */
    public TestIdSet() {
      elements = IntVector.empty();
    }

    /** Constructs a TestIdSet from a copy of the provided Collection of Integers. */
    public TestIdSet(Collection<Integer> ids) {
      elements = IntVector.copyOf(ids);
      elements.sort();
      elements.unique();
    }

    /** Constructs a TestIdSet from a copy of the provided array of Integers. */
    public TestIdSet(int... ids) {
      elements = IntVector.copyOf(ids);
      elements.sort();
      elements.unique();
    }

    @Override
    public OfInt intIterator() {
      return elements.intIterator();
    }

    @Override
    public void forEach(IntConsumer action) {
      elements.forEach(action);
    }

    @Override
    public void forEach(IntBiConsumer action) {
      elements.forEach(action);
    }

    @Override
    public int size() {
      return elements.size();
    }

    @Override
    public int first() {
      return elements.get(0);
    }
  }

  // Non-singleton, non-empty IdSets are stored in this SequenceLexicon.
  private final SequenceLexicon idSets;

  /** The unique fixed id for the empty set. The empty set takes up no space. */
  public static final int EMPTY_SET_ID = Integer.MIN_VALUE;

  /** Constructs an empty IdSetLexicon. */
  public IdSetLexicon() {
    idSets = new SequenceLexicon();
  }

  /** Copy constructor. */
  public IdSetLexicon(IdSetLexicon other) {
    idSets = new SequenceLexicon(other.idSets);
  }

  /** Clears all data from this IdSetLexicon. */
  public void clear() {
    idSets.clear();
  }

  /**
   * Returns true if this IdSetLexicon currently contains the same IdSets with the same ids as the
   * given IdSetLexicon 'other'. As IdSetLexicon assigns IdSet ids consecutively when IdSets are
   * added, this means that the IdSets were added in the same order. Does not override
   * Object.equals() because IdSetLexicon is mutable.
   */
  public boolean isEqualTo(IdSetLexicon other) {
    return idSets.equals(other.idSets);
  }

  /**
   * Convenience method that returns the unique id for a IdSetLexicon containing the given single
   * value.
   *
   * @throws IllegalArgumentException if "id" is negative.
   */
  public int addSingleton(int id) {
    Preconditions.checkArgument(id >= 0);
    return id;
  }

  /**
   * Add the given integers to the lexicon as a new set, if it is not already present. Returns the
   * unique id representing the set. Duplicates in the given "ids" are removed.
   *
   * @throws IllegalArgumentException if any value in "ids" is negative.
   */
  public int add(IntVector ids) {
    if (ids.isEmpty()) {
      // Empty sets have a special id chosen not to conflict with other ids.
      return EMPTY_SET_ID;
    }

    if (ids.size() == 1) {
      // Singleton sets are represented by their element.
      Preconditions.checkArgument(ids.get(0) >= 0);
      return ids.get(0);
    }

    // Sort and deduplicate.
    ids.sort();
    // No negative elements are allowed.
    Preconditions.checkArgument(ids.get(0) >= 0);
    ids.unique();

    // Might have become a singleton after deduplication.
    if (ids.size() == 1) {
      return ids.get(0);
    }

    // Non-singleton sets are represented by the bitwise complement of the id returned by
    // SequenceLexicon. The set is stored as a sequence in sorted order.
    return ~idSets.add(ids);
  }

  /**
   * Add the given 'ids' to the lexicon as a new, deduplicated set, if it is not already present.
   * Returns the unique id representing the set.
   *
   * @throws IllegalArgumentException if any value in "ids" is negative.
   */
  public int add(IntSequence ids) {
    return add(IntVector.copyOf(ids));
  }

  /**
   * Add the given Integers to the lexicon as a new, deduplicated set, if it is not already present.
   * Returns the unique id representing the set.
   *
   * @throws IllegalArgumentException if any value in "ids" is negative.
   */
  public int add(Iterable<Integer> ids) {
    return add(IntVector.copyOf(ids));
  }

  /**
   * Returns the IdSet corresponding to an id previously returned by add(), or a singleton or empty
   * IdSet.
   *
   * <p>Note that for empty or singleton sets, there is no requirement or check that the requested
   * set ever *was* previously added, as they are not actually stored at all: such IdSets are
   * generated on demand.
   *
   * @throws IllegalArgumentException if the given setId is for a multi-element set that was not
   *     previously added to the lexicon.
   */
  public IdSet idSet(int setId) {
    if (setId >= 0) {
      // A singleton set: generate on demand.
      return new IdSet() {
        int singletonValue = setId;

        @Override
        public OfInt intIterator() {
          return new OfInt() {
            int position = 0;

            @Override
            public int nextInt() {
              if (position == 0) {
                position++;
                return singletonValue;
              }
              throw new NoSuchElementException("nextInt() called on a singleton IdSet");
            }

            @Override
            public boolean hasNext() {
              return position == 0;
            }

            @Override
            public void forEachRemaining(IntConsumer consumer) {
              if (position == 0) {
                consumer.accept(singletonValue);
                position++;
              }
            }
          };
        }

        @Override
        public int first() {
          return singletonValue;
        }

        @Override
        public int size() {
          return 1;
        }

        @Override
        public boolean isEmpty() {
          return false;
        }

        @Override
        public void forEach(IntConsumer consumer) {
          consumer.accept(singletonValue);
        }

        @Override
        public void forEach(IntBiConsumer consumer) {
          consumer.accept(0, singletonValue);
        }
      };
    } else if (setId == EMPTY_SET_ID) {
      return new IdSet() {
        @Override
        public OfInt intIterator() {
          return new OfInt() {
            @Override
            public int nextInt() {
              throw new NoSuchElementException("nextInt() called on an empty IdSet");
            }

            @Override
            public boolean hasNext() {
              return false;
            }

            @Override
            public void forEachRemaining(IntConsumer consumer) {}
          };
        }

        @Override
        public int first() {
          throw new NoSuchElementException("first() called on an empty IdSet");
        }

        @Override
        public int size() {
          return 0;
        }

        @Override
        public boolean isEmpty() {
          return true;
        }

        @Override
        public void forEach(IntConsumer consumer) {}

        @Override
        public void forEach(IntBiConsumer consumer) {}
      };
    } else {
      return new IdSet() {
        // If the underlying sequence doesn't exist, throws an IllegalArgumentException.
        private final IntSequence sequence = idSets.sequence(~setId);

        @Override
        public OfInt intIterator() {
          return sequence.intIterator();
        }

        @Override
        public int first() {
          // TODO(torrey): If Sequence implemented get(), it would be more efficient to use get(0);
          return sequence.intIterator().nextInt();
        }

        @Override
        public int size() {
          return sequence.size();
        }

        @Override
        public boolean isEmpty() {
          return false;
        }

        @Override
        public void forEach(IntConsumer consumer) {
          sequence.forEach(consumer);
        }

        @Override
        public void forEach(IntBiConsumer consumer) {
          sequence.forEach(consumer);
        }
      };
    }
  }
}
