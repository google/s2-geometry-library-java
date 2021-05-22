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
import com.google.common.base.Preconditions;
import java.util.AbstractList;
import java.util.List;

/**
 * S2Shape is an abstract base class that defines a shape. Typically it wraps some other geometric
 * object in order to provide access to its edges without duplicating the edge data.
 */
@GwtCompatible
public interface S2Shape {
  /** Returns the number of edges in this shape. */
  int numEdges();

  /**
   * Returns the edge for the given index in {@code result}. Must not return zero-length edges.
   *
   * @param index which edge to set into {@code result}, from 0 to {@link #numEdges()} - 1
   */
  void getEdge(int index, MutableEdge result);

  /**
   * Returns true if this shape has an interior, i.e. the shape consists of one or more closed
   * non-intersecting loops.
   */
  boolean hasInterior();

  /**
   * Returns true if this shape contains {@link S2#origin()}. Should return false for shapes that do
   * not have an interior.
   */
  boolean containsOrigin();

  /**
   * A simple receiver for the endpoints of an edge.
   *
   * <p><>The {@link S2Edge} class is not suitable for retrieving large numbers of edges, as it
   * often triggers allocations. This class is intended to allow fast retrieval of the endpoints in
   * a single call.
   */
  final class MutableEdge {
    /**
     * Endpoints of this edge last set by passing this instance to {@link S2Shape#getEdge(int,
     * MutableEdge)}.
     */
    S2Point a;

    S2Point b;

    /**
     * Returns the leading point of the last edge retrieved via {@link S2Shape#getEdge(int,
     * MutableEdge)}, or null if no edge has been retrieved.
     */
    public S2Point getStart() {
      return a;
    }

    /**
     * Returns the trailing point of the last edge retrieved via {@link S2Shape#getEdge(int,
     * MutableEdge)}, or null if no edge has been retrieved.
     */
    public S2Point getEnd() {
      return b;
    }

    /** Returns true iff 'point' is either endpoint of this edge. */
    public boolean isEndpoint(S2Point point) {
      return a.equalsPoint(point) || b.equalsPoint(point);
    }

    /**
     * Called by implementations of {@link S2Shape#getEdge(int, MutableEdge)} to update the
     * endpoints of this mutable edge to the given values.
     */
    public void set(S2Point start, S2Point end) {
      this.a = start;
      this.b = end;
    }
  }

  /**
   * Returns the number of contiguous edge chains in the shape. For example, a shape whose edges are
   * [AB, BC, CD, AE, EF] may consist of two chains [A, B, C, D] and [A, E, F]. Every chain is
   * assigned a chain id numbered sequentially starting from zero.
   *
   * <p>An empty shape has no chains. A full shape (which contains the entire globe) has one chain
   * with no edges. Other shapes should have at least one chain, and the sum of all valid {@link
   * #getChainLength(int) chain lengths} should equal {@link #numEdges()} (that is, edges may only
   * be used by a single chain).
   *
   * <p>Note that it is always acceptable to implement this method by returning {@link #numEdges()}
   * (i.e. every chain consists of a single edge), but this may reduce the efficiency of some
   * algorithms.
   */
  int numChains();

  /**
   * Returns the first edge id corresponding to the edge chain for the given chain id. The edge
   * chains must form contiguous, non-overlapping ranges that cover the entire range of edge ids.
   *
   * @param chainId which edge chain to return its start, from 0 to {@link #numChains()} - 1
   */
  int getChainStart(int chainId);

  /**
   * Returns the number of edge ids corresponding to the edge chain for the given chain id. The edge
   * chains must form contiguous, non-overlapping ranges that cover the entire range of edge ids.
   *
   * @param chainId which edge chain to return its length, from 0 to {@link #numChains()} - 1
   */
  int getChainLength(int chainId);

  /**
   * Returns the edge for the given chain id and offset in {@code result}. Must not return
   * zero-length edges.
   *
   * @param chainId which chain contains the edge to return, from 0 to {@link #numChains()} - 1
   * @param offset position from chain start for the edge to return, from 0 to {@link
   *     #getChainLength(int)} - 1
   */
  void getChainEdge(int chainId, int offset, MutableEdge result);
  
  /**
   * Returns the start point of the edge that would be returned by {@link S2Shape#getChainEdge},
   * or the endpoint of the last edge if {@code edgeOffset==getChainLength(chainId)}.
   */
  S2Point getChainVertex(int chainId, int edgeOffset);

  /**
   * Returns a view of the vertices in the given chain. Note {@link S2Shape#dimension 2D} shapes
   * omit the last vertex, as it's a duplicate of the first.
   */
  default List<S2Point> chain(int chain) {
    return new AbstractList<S2Point>() {
      int length = getChainLength(chain) + (dimension() & 1);
      @Override public int size() {
        return length;
      }
      @Override public S2Point get(int index) {
        return getChainVertex(chain, index);
      }
    };
  }

  /** Returns a view of the {@link #chain chains} in this shape. */
  default List<List<S2Point>> chains() {
    return new AbstractList<List<S2Point>>() {
      @Override public int size() {
        return numChains();
      }
      @Override public List<S2Point> get(int index) {
        return chain(index);
      }
    };
  }

  /**
   * Returns the dimension of the geometry represented by this shape.
   *
   * <ul>
   *   <li>0 - Point geometry. Each point is represented as a degenerate edge.
   *   <li>1 - Polyline geometry. Polyline edges may be degenerate. A shape may represent any number
   *       of polylines. Polylines edges may intersect.
   *   <li>2 - Polygon geometry. Edges should be oriented such that the polygon interior is always
   *       on the left. In theory the edges may be returned in any order, but typically the edges
   *       are organized as a collection of edge chains where each chain represents one polygon
   *       loop. Polygons may have degeneracies, e.g., degenerate edges or sibling pairs consisting
   *       of an edge and its corresponding reversed edge. A polygon loop may also be full
   *       (containing all points on the sphere); by convention this is represented as a chain with
   *       no edges.
   * </ul>
   *
   * <p>Note that this method allows degenerate geometry of different dimensions to be
   * distinguished, e.g., it allows a point to be distinguished from a polyline or polygon that has
   * been simplified to a single point.
   */
  int dimension();

  /** Returns a point referenced to, i.e. indicating containment by, this shape. */
  default ReferencePoint getReferencePoint() {
    Preconditions.checkState(dimension() == 2);
    return ReferencePoint.create(S2.origin(), containsOrigin());
  }

  /** A point with a known containment relationship. */
  abstract class ReferencePoint extends S2Point {
    private static final ReferencePoint ORIGIN_INSIDE = create(S2.origin(), true);
    private static final ReferencePoint ORIGIN_OUTSIDE = create(S2.origin(), false);

    private ReferencePoint(S2Point p) {
      super(p.x, p.y, p.z);
    }

    /** Returns true if this point is contained by the reference shape. */
    public abstract boolean contained();

    /**
     * Returns a referenced point at an arbitrary position, suitable for shapes that contain all
     * points or no points.
     */
    public static ReferencePoint create(boolean contained) {
      return contained ? ORIGIN_INSIDE : ORIGIN_OUTSIDE;
    }

    /** Creates a referenced point at position 'p', with known containment 'contained'. */
    public static ReferencePoint create(S2Point p, boolean contained) {
      if (contained) {
        return new ReferencePoint(p) {
          @Override
          public boolean contained() {
            return true;
          }
        };
      } else {
        return new ReferencePoint(p) {
          @Override
          public boolean contained() {
            return false;
          }
        };
      }
    }

    @Override
    public boolean equals(Object o) {
      return o instanceof ReferencePoint
          && super.equals(o)
          && contained() == ((ReferencePoint) o).contained();
    }
  }
}
