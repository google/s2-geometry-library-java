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

import com.google.common.base.Preconditions;
import com.google.errorprone.annotations.InlineMe;
import java.util.AbstractList;
import java.util.List;
import java.util.Objects;
import jsinterop.annotations.JsConstructor;
import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;

/**
 * S2Shape is an abstract base class that defines a shape as a collection of edges, which are
 * organized in chains.
 *
 * <p>Typically it wraps some other geometric object in order to provide access to its edges without
 * duplicating the edge data.
 */
@JsType
@SuppressWarnings("Assertion")
public interface S2Shape {
  /**
   * Returns the number of edges in this shape. 0-dimensional shapes contain points, as one
   * degenerate edge per chain, with equal start and end vertices. 1-dimensional shape chains are
   * lines, with one vertex per edge, and an additional vertex at the end of the chain.
   * 2-dimensional shape chains are loops, with one vertex per edge.
   */
  int numEdges();

  /**
   * Returns the edge for the given edgeId in {@code result}. Points are represented as degenerate
   * edges, with equal endpoints, but not all degenerate edges are points.
   *
   * @param edgeId which edge to set into {@code result}, from 0 to {@link #numEdges()} - 1
   */
  void getEdge(int edgeId, MutableEdge result);

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
   * <p>The {@link S2Edge} class is not suitable for retrieving large numbers of edges, as it often
   * triggers allocations. This class is intended to allow fast retrieval of the endpoints in a
   * single call.
   */
  @JsType
  final class MutableEdge {
    /** The start point of this edge. */
    public S2Point a = null;

    /** The end point of this edge. */
    public S2Point b = null;

    /** Creates a new MutableEdge with the given endpoints. */
    public static MutableEdge of(S2Point a, S2Point b) {
      MutableEdge e = new MutableEdge();
      e.a = a;
      e.b = b;
      return e;
    }

    /** Returns the current start point of this edge, or null if it has not been set. */
    public S2Point getStart() {
      return a;
    }

    /** Returns the current end point of this edge, or null if it has not been set. */
    public S2Point getEnd() {
      return b;
    }

    /** Returns true if this MutableEdge is actually a point, i.e. the endpoints are equal. */
    public boolean isDegenerate() {
      return a.equalsPoint(b);
    }

    /**
     * Returns true iff 'point' is either endpoint of this edge.
     *
     * @deprecated use hasEndpoint.
     */
    @Deprecated
    @InlineMe(replacement = "this.hasEndpoint(point)")
    public boolean isEndpoint(S2Point point) {
      return hasEndpoint(point);
    }

    /** Returns true iff this edge has the given 'point' as either endpoint. */
    public boolean hasEndpoint(S2Point point) {
      return a.equalsPoint(point) || b.equalsPoint(point);
    }

    /**
     * Returns true if this MutableEdge currently has the same endpoints as the 'other' MutableEdge.
     * Does not override Object.equals() because MutableEdge is mutable.
     */
    public boolean isEqualTo(MutableEdge other) {
      return a.equalsPoint(other.a) && b.equalsPoint(other.b);
    }

    /** Returns true if this MutableEdge currently has the reversed endpoints as 'other'. */
    public boolean isSiblingOf(MutableEdge other) {
      return a.equalsPoint(other.b) && b.equalsPoint(other.a);
    }

    /** Updates the endpoints of this mutable edge to the given values. */
    public void set(S2Point start, S2Point end) {
      this.a = start;
      this.b = end;
    }

    /** Exchanges the endpoints of this edge. */
    public void reverse() {
      S2Point t = a;
      a = b;
      b = t;
    }

    public String toDegreesString() {
      return (a == null ? "null" : a.toDegreesString())
          + "-"
          + (b == null ? "null" : b.toDegreesString());
    }

    @Override
    public String toString() {
      return toDegreesString();
    }
  }

  /**
   * The position of an edge within a given S2Shape's edge chains, specified as a (chainId, offset)
   * pair. Chains are numbered sequentially starting from zero, and offsets are measured from the
   * start of each chain. ChainPosition is mutable, as it is intended to be reused.
   */
  @JsType
  final class ChainPosition {
    int chainId;
    int offset;

    /** Sets this ChainPosition's chainId and offset. */
    public void set(int chainId, int offset) {
      this.chainId = chainId;
      this.offset = offset;
    }

    /** Gets this ChainPosition's chainId. */
    public int getChainId() {
      return chainId;
    }

    /** Gets this ChainPosition's offset (how far from the start of the chain it is, from 0). */
    public int getOffset() {
      return offset;
    }

    /**
     * True if this ChainPosition has the same chainId and offset as the 'other' ChainPosition. Does
     * not override Object.equals() because ChainPosition is mutable.
     */
    public boolean isEqualTo(ChainPosition other) {
      return chainId == other.chainId && offset == other.offset;
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
   * Finds the chain containing the given edge, and returns the position of that edge as a (chainId,
   * offset) pair in {@code result}.
   */
  void getChainPosition(int edgeId, ChainPosition result);

  /**
   * Returns the start point of the edge that would be returned by {@link S2Shape#getChainEdge}, or
   * the endpoint of the last edge if {@code edgeOffset==getChainLength(chainId)}.
   */
  S2Point getChainVertex(int chainId, int edgeOffset);

  /**
   * Returns the edge id of the next edge in a chain. Wraps around at the start/end of any closed
   * chains.
   *
   * <p>This is intended for one-off lookups, as it has to look up the chain for the edge every
   * time. If you want many lookups or to iterate the edges of a chain, then it's better to do that
   * directly.
   *
   * <p>Return -1 when the end of an open chain is reached. Polygon and closed polyline chains wrap
   * around to the beginning and thus never return -1, while points always do.
   */
  default int nextEdgeWrap(int edgeId) {
    assert edgeId >= 0;
    assert edgeId < numEdges();

    ChainPosition chainPos = new ChainPosition();
    getChainPosition(edgeId, chainPos);
    int chainLength = getChainLength(chainPos.chainId);

    // Polygon chains wrap around, point and polylines don't.
    int offset = chainPos.offset;
    switch (dimension()) {
      case 2:
        // Polygon chains always wrap around.
        offset = (offset + 1) % chainLength;
        break;
      case 1:
        // If we're at the end of a polyline, wrap around if it's closed.
        offset++;
        if (offset == chainLength) {
          S2Point curr = getChainVertex(chainPos.chainId, chainLength);
          S2Point next = getChainVertex(chainPos.chainId, 0);
          if (curr.equalsPoint(next)) {
            offset = 0;
          } else {
            return -1;
          }
        }
        break;
      default:
        // Points are one per chain.
        return -1;
    }

    return getChainStart(chainPos.chainId) + offset;
  }

  /**
   * Returns the edge id of previous edge in a chain. Wraps around at the start/end of any closed
   * chains/
   *
   * <p>This is intended for one-off lookups, as it has to look up the chain for the edge every
   * time. If you want many lookups or to iterate the edges of a chain, then it's better to do that
   * directly.
   *
   * <p>Return -1 when the start of an open chain is reached. Polygon and closed polyline chains
   * wrap around to the end and thus never return -1, while points always do.
   */
  default int prevEdgeWrap(int edgeId) {
    assert edgeId >= 0;
    assert edgeId < numEdges();

    ChainPosition chainPos = new ChainPosition();
    getChainPosition(edgeId, chainPos);
    int chainLength = getChainLength(chainPos.chainId);

    int offset = chainPos.offset;
    switch (dimension()) {
      case 2:
        // Polygons always wrap around.
        --offset;
        if (offset < 0) {
          offset += chainLength;
        }
        break;
      case 1:
        --offset;
        if (offset < 0) {
          S2Point curr = getChainVertex(chainPos.chainId, 0);
          S2Point prev = getChainVertex(chainPos.chainId, chainLength);
          if (prev.equalsPoint(curr)) {
            offset += chainLength;
          } else {
            return -1;
          }
        }
        break;
      default:
        // Points are one per chain.
        return -1;
    }

    return getChainStart(chainPos.chainId) + offset;
  }

  /**
   * Returns a view of the vertices in the given chain. Note {@link S2Shape#dimension 2D} shapes
   * omit the last vertex, as it's a duplicate of the first.
   */
  default List<S2Point> chain(int chain) {
    return new AbstractList<S2Point>() {
      int length = getChainLength(chain) + (dimension() & 1);

      @Override
      public int size() {
        return length;
      }

      @Override
      public S2Point get(int index) {
        return getChainVertex(chain, index);
      }
    };
  }

  /** Returns a view of the {@link #chain chains} in this shape. */
  default List<List<S2Point>> chains() {
    return new AbstractList<List<S2Point>>() {
      @Override
      public int size() {
        return numChains();
      }

      @Override
      public List<S2Point> get(int index) {
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

  /**
   * Returns true if the shape contains no points. (Note that the full polygon is represented as a
   * chain with zero edges.)
   */
  default boolean isEmpty() {
    return numEdges() == 0 && (dimension() < 2 || numChains() == 0);
  }

  /** Returns true if the shape contains all points on the sphere and has no edges. */
  default boolean isFull() {
    return numEdges() == 0 && dimension() == 2 && numChains() > 0;
  }

  /** Returns a point referenced to, i.e. indicating containment by, this shape. */
  default ReferencePoint getReferencePoint() {
    Preconditions.checkState(dimension() == 2);
    return ReferencePoint.create(S2.origin(), containsOrigin());
  }

  /** A point with a known containment relationship. */
  @JsType
  abstract class ReferencePoint {
    private static final ReferencePoint ORIGIN_INSIDE = create(S2.origin(), true);
    private static final ReferencePoint ORIGIN_OUTSIDE = create(S2.origin(), false);
    private final S2Point point;

    @JsConstructor
    private ReferencePoint(S2Point p) {
      this.point = p;
    }

    /** Returns the point referenced to. */
    public S2Point point() {
      return point;
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
    @JsMethod(name = "createfromPoint")
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

    /** Returns true if this ReferencePoint's point is equal to the given point. */
    public boolean equalsPoint(S2Point p) {
      return point.equalsPoint(p);
    }

    @Override
    public boolean equals(Object o) {
      if (!(o instanceof ReferencePoint)) {
        return false;
      }
      ReferencePoint that = (ReferencePoint) o;
      return this.point.equalsPoint(that.point) && this.contained() == that.contained();
    }

    @Override
    public int hashCode() {
      return Objects.hash(point, contained());
    }
  }
}
