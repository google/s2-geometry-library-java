/*
 * Copyright 2014 Google Inc.
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

import com.google.common.annotations.GwtCompatible;

/**
 * S2Shape is an abstract base class that defines a shape. Typically it wraps some other geometric
 * object in order to provide access to its edges without duplicating the edge data. Shapes are
 * immutable once they have been indexed; to modify a shape it must be removed and then reinserted.
 * A shape can be removed from one S2ShapeIndex and then inserted into another, but it can belong to
 * only one index at a time.
 */
@GwtCompatible
public interface S2Shape {
  /**
   * Returns the number of edges in this shape.
   */
  public int numEdges();

  /**
   * Returns the edge for the given index in {@code result}. Must not return zero-length edges.
   *
   * @param index which edge to set into {@code result}, from 0 to {@link #numEdges()} - 1
   */
  public void getEdge(int index, MutableEdge result);

  /**
   * Returns true if this shape has an interior, i.e. the shape consists of one or more closed
   * non-intersecting loops.
   */
  public abstract boolean hasInterior();

  /**
   * Returns true if this shape contains {@link S2#origin()}. Should return false for shapes that do
   * not have an interior.
   */
  public abstract boolean containsOrigin();

  /**
   * A simple receiver for the endpoints of an edge.
   *
   * <>The {@link S2Edge} class is not suitable for retrieving large numbers of edges, as it often
   * triggers allocations. This class is intended to allow fast retrieval of the endpoints in a
   * single call.
   */
  public static final class MutableEdge {
    private S2Point a;
    private S2Point b;

    /**
     * Returns the leading point of the last edge retrieved via
     * {@link S2Shape#getEdge(int, MutableEdge)}, or null if no edge has been retrieved.
     */
    public S2Point a() {
      return a;
    }

    /**
     * Returns the trailing point of the last edge retrieved via
     * {@link S2Shape#getEdge(int, MutableEdge)}, or null if no edge has been retrieved.
     */
    public S2Point b() {
      return b;
    }

    /**
     * Called by implementations of {@link S2Shape#getEdge(int, MutableEdge)} to update the
     * endpoints of this mutable edge to the given values.
     */
    public void set(S2Point a, S2Point b) {
      this.a = a;
      this.b = b;
    }
  }
}
