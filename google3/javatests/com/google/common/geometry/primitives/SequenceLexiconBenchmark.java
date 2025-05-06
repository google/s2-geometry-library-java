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

import static java.util.concurrent.TimeUnit.NANOSECONDS;
import static java.util.concurrent.TimeUnit.SECONDS;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.benchmarks.S2BenchmarkBaseState;
import com.google.common.geometry.primitives.Ints.IntSequence;
import com.google.common.geometry.primitives.Ints.OfInt;
import java.util.ArrayList;
import java.util.Random;
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

/** Benchmarks for SequenceLexicon. */
public class SequenceLexiconBenchmark {
  private SequenceLexiconBenchmark() {}

  /**
   * Benchmark state for measuring performance of adding sequences to a SequenceLexicon.
   */
  @State(Scope.Thread)
  @BenchmarkMode(Mode.AverageTime)
  @OutputTimeUnit(NANOSECONDS)
  @Warmup(iterations = 3, time = 10, timeUnit = SECONDS)
  @Measurement(iterations = 5, time = 10, timeUnit = SECONDS)
  public static class RandomSequences extends S2BenchmarkBaseState {
    private static final int SEED = 123455;
    private final Random rand = new Random(SEED);

    // How many sequences to test with?
    private static final int NUM_SEQUENCES = 1 << 20; // Over 1 million
    private static final int MASK_SEQUENCES = NUM_SEQUENCES - 1;

    // There could be duplicate sequences, so use half as many ids.
    private static final int NUM_IDS = 1 << 19;
    private static final int MASK_IDS = NUM_IDS - 1;

    // How long are the sequences?
    @Param({"1", "2", "10"})
    int sequenceLength;

    /** Generates a list of random Integers of the given "length". */
    public ImmutableList<Integer> generateRandomSequence(int length) {
      ImmutableList.Builder<Integer> b = new ImmutableList.Builder<>();
      for (int i = 0; i < length; i++) {
        b.add(rand.nextInt());
      }
      return b.build();
    }

    // numSequences pregenerated random sequences.
    protected ArrayList<ImmutableList<Integer>> randomSequences = new ArrayList<>();

    // NUM_IDS sequence ids.
    protected ArrayList<Integer> sequenceIds = new ArrayList<>();
    // Which sequenceId to use next.
    protected int index;

    // A lexicon with no sequences at the beginning of an iteration.
    protected SequenceLexicon emptyLex;
    // A lexicon with "numSequences" sequences at the beginning of an iteration.
    protected SequenceLexicon fullLex;

    // At the beginning of each trial, generate NUM_SEQUENCES random sequences.
    @Setup(Level.Trial)
    public void setupSequences() {
      // Print warnings if boxing methods are called, ie. next() instead of nextInt().
      System.setProperty("org.openjdk.java.util.stream.tripwire", "true");

      for (int i = 0; i < NUM_SEQUENCES; i++) {
        randomSequences.add(generateRandomSequence(sequenceLength));
      }
    }

    // At the beginning of each iteration, create new SequenceLexicons, one empty and one containing
    // NUM_SEQUENCES sequences.
    @Setup(Level.Iteration)
    public void setupLexicons() {
      index = 0;

      emptyLex = new SequenceLexicon();
      fullLex = new SequenceLexicon();

      // Pre-add all the random sequences to fullLex, and save the ids for doing lookups.
      for (int i = 0; i < NUM_SEQUENCES; i++) {
        int id = fullLex.add(randomSequences.get(i));
        sequenceIds.add(id);
      }

      // Ensure there are at least the required number of unique sequences.
      Preconditions.checkState(fullLex.size() > NUM_IDS);
    }

    /**
     * Measures the time required to add random sequences to a lexicon that begins the iteration
     * empty.
     */
    @Benchmark
    public int addSequenceStartingFromEmpty() {
      return emptyLex.add(randomSequences.get(index++ & MASK_SEQUENCES));
    }

    /**
     * Measures the time required to add random sequences to a lexicon that contains NUM_SEQUENCES
     * sequences at the beginning of the iteration.
     */
    @Benchmark
    public int addSequenceStartingFromFull() {
      return fullLex.add(randomSequences.get(index++ & MASK_SEQUENCES));
    }

    /** Measures the time required to get an existing sequence from a lexicon. */
    @Benchmark
    public IntSequence getSequence() {
      int sequenceId = sequenceIds.get(index++ & MASK_IDS);
      return fullLex.sequence(sequenceId);
    }

    /**
     * Measures the time required to get a sequence from a lexicon and iterate over its values with
     * an iterator OfInt using hasNext() and nextInt().
     */
    @Benchmark
    public long getAndIterateSequenceWithNext() {
      int sequenceId = sequenceIds.get(index++ & MASK_IDS);
      OfInt sequenceIter = fullLex.sequence(sequenceId).intIterator();
      long total = 0;
      while (sequenceIter.hasNext()) {
        total += sequenceIter.nextInt();
      }
      return total;
    }
  }
}
