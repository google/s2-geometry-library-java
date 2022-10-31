/*
 * Copyright 2021 Google Inc.
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

import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;

/**
 * DistanceCollector is an interface for working with abstract distances, tracking the current
 * "best" distance seen over a sequence of update calls. This can be better than simply computing
 * the distance and retaining the best, because the current best can often be used to determine that
 * the new distance cannot be better more cheaply than actually determining the precise distance.
 *
 * <p>The meaning of "best" distance is left up to the implementation. The provided update methods
 * support finding a best distance between S2 points, edges, and cells. Current implementations
 * support finding minimum and maximum distances.
 */
@JsType
public interface DistanceCollector<T extends S1Distance<T>>
    extends Comparable<DistanceCollector<T>> {
  /**
   * Returns the current best distance. If distance() is called on a DistanceCollector for which
   * update() and set() have not been called since the collector was constructed or reset, the
   * distance will be an invalid value, worse than any valid distance. For example, the distance
   * will be {@link S1ChordAngle#NEGATIVE} for a collector of maximum S1ChordAngles.
   */
  public T distance();

  /** Resets this collector to the default, worst value. */
  void reset();

  /** Sets this collector to the given distance value. */
  void set(T value);

  /** This default implementation of Comparable relies on S1Distance being comparable. */
  @Override
  public default int compareTo(DistanceCollector<T> other) {
    return distance().compareTo(other.distance());
  }

  /**
   * Update this collector to the better of the current distance and the given 'other' distance.
   * Returns true if this distance was updated, false otherwise.
   */
  boolean update(T other);

  /**
   * Update this collector to the better of the current distance, and the distance between the two
   * given points. Returns true if this distance was updated, false otherwise.
   */
  @JsMethod(name = "updatePointToPoint")
  boolean update(S2Point p1, S2Point p2);

  /**
   * Update this collector to the better of the current distance, and the distance between the
   * given point and edge. Returns true if this distance was updated, false otherwise.
   */
  @JsMethod(name = "updatePointToEdge")
  boolean update(S2Point p, S2Point v0, S2Point v1);

  /**
   * Update this collector to the better of the current distance, and the distance between the two
   * given edges. Returns true if this distance was updated, false otherwise.
   */
  @JsMethod(name = "updateEdgeToEdge")
  boolean update(S2Point v0, S2Point v1, S2Point w0, S2Point w1);

  /**
   * Update this collector to the better of the current distance, and the distance between the
   * given point and cell. Returns true if this distance was updated, false otherwise.
   */
  @JsMethod(name = "updatePointToCell")
  boolean update(S2Point p, S2Cell c);

  /**
   * Update this collector to the better of the current distance, and the distance between the
   * given edge and cell. Returns true if this distance was updated, false otherwise.
   */
  @JsMethod(name = "updateEdgeToCell")
  boolean update(S2Point v0, S2Point v1, S2Cell c);

  /**
   * Update this collector to the better of the current distance, and the distance between the two
   * given cells. Returns true if this distance was updated, false otherwise.
   */
  @JsMethod(name = "updateCellToCell")
  boolean update(S2Cell c1, S2Cell c2);
}
