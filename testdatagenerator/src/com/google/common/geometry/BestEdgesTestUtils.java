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
package com.google.common.geometry;

import static java.util.Comparator.comparingDouble;

import com.google.common.collect.Ordering;
import com.google.common.geometry.S2BestEdgesQueryBase.Result;
import java.util.Comparator;
import java.util.List;

/**
 * Utility methods for testing S2ClosestEdgeQuery, S2FurthestEdgeQuery, and their common base class
 * S2BestEdgesQueryBase.
 *
 * <p>WARNING: These APIs are not stable and not intended for use outside of S2 benchmarks and unit
 * tests.
 */
public final class BestEdgesTestUtils {

  private BestEdgesTestUtils() { }

  /** Returns the given S2BestEdgesQueryBase.{@code Result<S1ChordAngle>} as a String. */
  static String resultToString(Result<S1ChordAngle> result) {
    StringBuilder output = new StringBuilder("distance ");
    output
        .append(result.distance())
        .append(" to ")
        .append(result.shape().getClass().getName())
        .append("#")
        .append(System.identityHashCode(result.shape()))
        .append(", edge#")
        .append(result.edgeId());
    return output.toString();
  }

  /** Returns the contents of a list of S2BestEdgesQueryBase.Results as a multi-line String. */
  static String resultsToString(List<Result<S1ChordAngle>> results) {
    StringBuilder output = new StringBuilder("Result List: \n");
    int i = 0;
    for (Result<S1ChordAngle> result : results) {
      output
          .append("  Result #")
          .append(i++)
          .append(" is: ")
          .append(resultToString(result))
          .append("\n");
    }
    return output.toString();
  }

  /** Handles differences between closestEdges and furthestEdges for testing. */
  public enum MinOrMax {
    /**
     * For finding minimum distances, better is smaller and results should be ordered from smaller
     * to larger distances.
     */
    MIN {
      @Override
      public S1ChordAngle bestDistance() {
        return S1ChordAngle.ZERO;
      }

      @Override
      public S1ChordAngle worstDistance() {
        return S1ChordAngle.INFINITY;
      }

      @Override
      public Comparator<Result<S1ChordAngle>> resultOrdering() {
        return comparingDouble((Result<S1ChordAngle> result) -> result.distance().getLength2());
      }

      @Override
      public S1ChordAngle padTowardBest(S1ChordAngle value, S1ChordAngle padding) {
        return S1ChordAngle.sub(value, padding);
      }

      @Override
      public boolean distanceWithinLimit(S1ChordAngle distance, S1ChordAngle limit) {
        return distance.lessThan(limit);
      }
    },

    /**
     * For finding maximum distances, better is larger, and results should be ordered from larger to
     * smaller distances.
     */
    MAX {
      @Override
      public S1ChordAngle bestDistance() {
        return S1ChordAngle.STRAIGHT;
      }

      @Override
      public S1ChordAngle worstDistance() {
        return S1ChordAngle.NEGATIVE;
      }

      @Override
      public Comparator<Result<S1ChordAngle>> resultOrdering() {
        return comparingDouble((Result<S1ChordAngle> result) -> result.distance().getLength2())
            .reversed();
      }

      @Override
      public S1ChordAngle padTowardBest(S1ChordAngle value, S1ChordAngle padding) {
        return S1ChordAngle.add(value, padding);
      }

      @Override
      public boolean distanceWithinLimit(S1ChordAngle distance, S1ChordAngle limit) {
        return distance.greaterThan(limit);
      }
    };

    /** Best possible valid distance. */
    public abstract S1ChordAngle bestDistance();

    /** Worst distance, not a valid distance. */
    public abstract S1ChordAngle worstDistance();

    /** Expected ordering of Results by distance. */
    public abstract Comparator<Result<S1ChordAngle>> resultOrdering();

    /** Pad the given value by the given padding, in the direction of better distances. */
    public abstract S1ChordAngle padTowardBest(S1ChordAngle value, S1ChordAngle padding);

    /** Returns true if the given distance doesn't exceed (isn't worse than) the given limit. */
    public abstract boolean distanceWithinLimit(S1ChordAngle distance, S1ChordAngle limit);
  };

  /**
   * Compare two sets of best distance results, where "bruteForce" is computed via brute force
   * (i.e., considering every possible edge) and "optimized" is computed using a spatial data
   * structure. Here "maxResults" is a bound on the maximum number of items, "maxDistance" is a
   * limit on the distance to any item, and "maxError" is the maximum error allowed when selecting
   * which items are closest. Returns true if and only if the two sets of Results are equivalent.
   */
  @SuppressWarnings("ShortCircuitBoolean")
  public static boolean resultsAreEquivalent(
      MinOrMax minOrMax,
      List<Result<S1ChordAngle>> bruteForce,
      List<Result<S1ChordAngle>> optimized,
      int maxResults,
      S1ChordAngle maxDistance,
      S1ChordAngle maxError) {
    // This is a conservative bound on the error in computing the distance from the target geometry
    // to an S2Cell. Such errors can cause candidates to be pruned from the result set even though
    // they may be slightly closer.
    final S1ChordAngle kMaxPruningError = S1ChordAngle.fromRadians(1e-15);

    // Are there brute-force results that are not present in the optimized results?
    boolean optimizedSubsetOfBruteForce =
        checkResultSet(
            minOrMax,
            optimized,
            bruteForce,
            maxResults,
            maxDistance,
            maxError,
            kMaxPruningError,
            "Optimized (X) vs. Brute Force (Y): ");

    // Are there optimized results that are not present in the brute-force results?
    boolean bruteForceSubsetOfOptimized =
        checkResultSet(
            minOrMax,
            bruteForce,
            optimized,
            maxResults,
            maxDistance,
            maxError,
            S1ChordAngle.ZERO,
            "Brute Force (X) vs. Optimized (Y): ");

    // Deliberate '&', not short-circuit as we want to check both cases.
    return (optimizedSubsetOfBruteForce & bruteForceSubsetOfOptimized);
  }

  /**
   * Check that result set "x" contains all the expected results from "y", with distance better than
   * the given maxDistance, and does not include any duplicate results.
   */
  public static boolean checkResultSet(
      MinOrMax minOrMax,
      List<Result<S1ChordAngle>> x,
      List<Result<S1ChordAngle>> y,
      int maxResults,
      S1ChordAngle distanceLimit,
      S1ChordAngle maxError,
      S1ChordAngle maxPruningError,
      String label) {

    // Both sets of Results should be sorted by distance from best to worst.
    if (!Ordering.from(minOrMax.resultOrdering()).isOrdered(x)) {
      System.err.println("X Results are not ordered by distance: " + resultsToString(x));
      return false;
    }

    if (!Ordering.from(minOrMax.resultOrdering()).isOrdered(y)) {
      System.err.println("Y Results are not ordered by distance: " + resultsToString(y));
      return false;
    }

    StringBuilder out = new StringBuilder();
    out.append(minOrMax)
        .append(" distance edges.\n")
        .append("Checking for cases where X is missing something in Y for ")
        .append(label)
        .append(". There are ")
        .append(x.size())
        .append(" results in X and ")
        .append(y.size())
        .append(" results in Y.\n")
        .append("All X results: \n")
        .append(resultsToString(x))
        .append("\n")
        .append("All Y results: \n")
        .append(resultsToString(y))
        .append("\n");

    // Result set X should contain all the items from Y whose distance is better than adjustedLimit.
    S1ChordAngle adjustedLimit = minOrMax.bestDistance();

    if (x.size() < maxResults) {
      // Result set X was not limited by "maxResults", so it should contain all the items up to
      // "distanceLimit", except that a few items very near the limit may be missed because the
      // distance measurements used for pruning S2Cells are not conservative.
      if (distanceLimit.equals(minOrMax.worstDistance())) {
        adjustedLimit = distanceLimit;
      } else {
        adjustedLimit = minOrMax.padTowardBest(distanceLimit, maxPruningError);
      }
      out.append("Note x.size() < maxResults (")
          .append(maxResults)
          .append(") so x was not limited by maxResults, and should contain ")
          .append(minOrMax)
          .append(" distance results up to ")
          .append(adjustedLimit)
          .append("\n");
    } else if (!x.isEmpty()) {
      // Result set X contains only the best "maxResults" items, to within a tolerance of
      // "maxError + maxPruningError".
      S1ChordAngle xWorstDistance = x.get(x.size() - 1).distance();
      adjustedLimit =
          minOrMax.padTowardBest(minOrMax.padTowardBest(xWorstDistance, maxError), maxPruningError);
      out.append("Note x.size() has only maxResults (")
          .append(maxResults)
          .append(") results, which should be the ")
          .append(minOrMax)
          .append(" distance ones within adjustedLimit ")
          .append(adjustedLimit)
          .append("\n");
    }

    // Check that X has a result that's matching or an acceptable replacement for each result in Y.
    boolean result = true;
    for (Result<S1ChordAngle> yp : y) {
      // Count the number of entries in x that have a shape and edgeId equal to yp. Only considers
      // entries that are within the adjustedLimit computed above, so if (for example) an entry in
      // x is for a different shape or edge, but has the ~ same distance, that's acceptable as well.
      // Note that this test also catches duplicate values.
      long count =
          x.stream()
              .filter(xp -> (xp.edgeId() == yp.edgeId()) && (xp.shape() == yp.shape()))
              .count();
      if (minOrMax.distanceWithinLimit(yp.distance(), adjustedLimit) && count != 1) {
        result = false;
        out.append((count > 1 ? "Duplicate" : "Missing"))
            .append(" Y result: ")
            .append(resultToString(yp))
            .append(", X count = ")
            .append(count)
            .append(", Distance is within adjustedLimit ")
            .append(adjustedLimit)
            .append("\n");
      }
    }
    if (!result) {
      System.err.println("Test failed:\n" + out + "\n");
    }
    return result;
  }

}
