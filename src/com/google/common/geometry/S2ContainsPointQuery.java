/*
 * Copyright 2017 Google Inc.
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
import com.google.common.base.Predicate;
import com.google.common.collect.AbstractIterator;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import java.util.Iterator;

/**
 * A query for whether one or more shapes in an {@link S2ShapeIndex} contain a given S2Point.
 *
 * <p>The S2ShapeIndex may contain any number of points, polylines, and/or polygons (possibly
 * overlapping). Shape boundaries are modeled with a constructor parameter, {@link S2VertexModel},
 * which defaults to {@link S2VertexModel#SEMI_OPEN}. This may be customized to control whether or
 * not shapes are considered to contain their vertices.
 *
 * <p>This class is not thread-safe. To use it in parallel, each thread should construct its own
 * instance (this is not expensive). However, note that if you need to do a large number of point
 * containment tests, it is more efficient to re-use the S2ContainsPointQuery object rather than
 * constructing a new one each time.
 */
@GwtCompatible
public class S2ContainsPointQuery {
  /** The options for building an S2ContainsPointQuery. */
  public static final class Options {
    public static final Options OPEN = new Options(S2VertexModel.OPEN);
    public static final Options SEMI_OPEN = new Options(S2VertexModel.SEMI_OPEN);
    public static final Options CLOSED = new Options(S2VertexModel.CLOSED);
    private final S2VertexModel vertexModel;

    private Options(S2VertexModel vertexModel) {
      this.vertexModel = vertexModel;
    }
    /** Returns the vertex model in this options. */
    public S2VertexModel vertexModel() {
      return vertexModel;
    }
  }

  /** A rule for whether shapes are considered to contain their vertices. */
  public enum S2VertexModel {
    /**
     * In the OPEN model, no shapes contain their vertices (not even points). Therefore
     * contains(S2Point) returns true if and only if the point is in the interior of some polygon.
     */
    OPEN,

    /**
     * In the SEMI_OPEN model, polygon point containment is defined such that if several polygons
     * tile the region around a vertex, then exactly one of those polygons contains that vertex.
     * Points and polylines still do not contain any vertices.
     */
    SEMI_OPEN,

    /** In the CLOSED model, all shapes contain their vertices (including points and polylines). */
    CLOSED;

    /**
     * Returns true if the clipped portion of a shape 'clipped' from a cell with center 'cellCenter'
     * contains the point 'p' according to this vertex model.
     */
    public boolean shapeContains(S2Point cellCenter, S2ClippedShape clipped, S2Point p) {
      boolean inside = clipped.containsCenter();
      int numEdges = clipped.numEdges();
      if (numEdges > 0) {
        // Points and polylines can be ignored unless the vertex model is CLOSED.
        S2Shape shape = clipped.shape();
        if (!shape.hasInterior() && this != S2VertexModel.CLOSED) {
          return false;
        }
        // Test containment by drawing a line segment from the cell center to the
        // given point and counting edge crossings.
        EdgeCrosser crosser = new EdgeCrosser(cellCenter, p);
        MutableEdge edge = new MutableEdge();
        boolean crossing;
        for (int i = 0; i < numEdges; ++i) {
          shape.getEdge(clipped.edge(i), edge);
          switch (crosser.robustCrossing(edge.a, edge.b)) {
            case -1:
              // Disjoint, so advance to next edge.
              continue;
            case 1:
              // Definitely crossing.
              crossing = true;
              break;
            default:
              // Shared vertex, test if we have a vertex crossing.
              // For the OPEN and CLOSED models, check whether "p" is a vertex.
              if (this != S2VertexModel.SEMI_OPEN && edge.isEndpoint(p)) {
                return this == S2VertexModel.CLOSED;
              }
              crossing = S2EdgeUtil.vertexCrossing(cellCenter, p, edge.a, edge.b);
              break;
          }
          inside ^= crossing;
        }
      }
      return inside;
    }
  }

  private final Options options;
  private final S2Iterator<S2ShapeIndex.Cell> it;

  /** Constructs a semi-open contains-point query from the given iterator. */
  public S2ContainsPointQuery(S2ShapeIndex index) {
    this(index, Options.SEMI_OPEN);
  }

  /** Constructs a contains-point query from the given iterator, with the specified options. */
  public S2ContainsPointQuery(S2ShapeIndex index, Options options) {
    this.it = index.iterator();
    this.options = options;
  }

  /** Returns the options used to build this query. */
  public Options options() {
    return options;
  }

  /**
   * Returns true if any shape in the given iterator contains {@code p} under the specified {@link
   * S2VertexModel}.
   */
  public boolean contains(S2Point p) {
    if (!it.locate(p)) {
      return false;
    }
    S2ShapeIndex.Cell cell = it.entry();
    S2Point center = it.center();
    int numClipped = cell.numShapes();
    for (int s = 0; s < numClipped; ++s) {
      if (options.vertexModel().shapeContains(center, cell.clipped(s), p)) {
        return true;
      }
    }
    return false;
  }

  /**
   * Returns true if the given shape contains {@code p} under the specified {@link S2VertexModel}.
   */
  public boolean shapeContains(S2Shape shape, S2Point p) {
    if (!it.locate(p)) {
      return false;
    }
    S2ClippedShape clipped = it.entry().findClipped(shape);
    if (clipped == null) {
      return false;
    }
    return options.vertexModel().shapeContains(it.center(), clipped, p);
  }

  /**
   * A visitor that receives each shape that contains a query point, returning true to continue
   * receiving shapes or false to terminate early.
   */
  interface ShapeVisitor extends Predicate<S2Shape> {}

  /**
   * Visits each shape that contains {@code p} under the specified {@link S2VertexModel} exactly
   * once, and returns true, or terminates early and returns false if any invocation of {@link
   * ShapeVisitor#apply(S2Shape)} returns false.
   */
  boolean visitContainingShapes(S2Point p, ShapeVisitor visitor) {
    // This function returns "false" only if the algorithm terminates early because the "visitor"
    // function returned false.
    if (!it.locate(p)) {
      return true;
    }
    S2ShapeIndex.Cell cell = it.entry();
    S2Point center = it.center();
    int numClipped = cell.numShapes();
    for (int s = 0; s < numClipped; ++s) {
      S2ClippedShape clipped = cell.clipped(s);
      if (options.vertexModel().shapeContains(center, clipped, p)
          && !visitor.apply(clipped.shape())) {
        return false;
      }
    }
    return true;
  }

  /** A convenience function that returns all the shapes that contain {@code p}. */
  public Iterable<S2Shape> getContainingShapes(final S2Point p) {
    if (!it.locate(p)) {
      return ImmutableList.of();
    } else {
      // Must copy the iterator immediately, since the Iterable may not be used until after this.it
      // has been repositioned.
      final S2ShapeIndex.Cell cell = it.entry();
      final S2Point center = it.center();
      return new Iterable<S2Shape>() {
        @Override
        public Iterator<S2Shape> iterator() {
          return new AbstractIterator<S2Shape>() {
            int i = 0;

            @Override
            protected S2Shape computeNext() {
              while (i < cell.numShapes()) {
                S2ClippedShape clipped = cell.clipped(i++);
                if (options.vertexModel().shapeContains(center, clipped, p)) {
                  return clipped.shape();
                }
              }
              return endOfData();
            }
          };
        }
      };
    }
  }

  /** A visitor that receives each edge that has some query point p as an endpoint. */
  interface EdgeVisitor {
    /**
     * Returns true if the next edge should be received, or false to terminate early.
     *
     * @param shape The shape of this edge
     * @param edgeId The edge ID in 'shape' that produced this edge
     * @param a The startpoint of the edge
     * @param b The endpoint of the edge
     */
    boolean test(S2Shape shape, int edgeId, S2Point a, S2Point b);
  }

  /**
   * Visits each edge in the index that is incident to {@code p} exactly once, and returns true, or
   * terminates early and returns false if {@code visitor} returns false. An "incident edge" is one
   * where {@code p} is one of the edge endpoints. The visitor requires the edge endpoints, and so
   * this method requires a temporary mutable edge to store edges in.
   */
  boolean visitIncidentEdges(S2Point p, EdgeVisitor visitor, MutableEdge tmp) {
    // This function returns "false" only if the algorithm terminates early because the "visitor"
    // function returned false.
    if (!it.locate(p)) {
      return true;
    }
    S2ShapeIndex.Cell cell = it.entry();
    int numClipped = cell.numShapes();
    for (int s = 0; s < numClipped; s++) {
      S2ClippedShape clipped = cell.clipped(s);
      int numEdges = clipped.numEdges();
      if (numEdges == 0) {
        continue;
      }
      S2Shape shape = clipped.shape();
      for (int i = 0; i < numEdges; i++) {
        int edgeId = clipped.edge(i);
        shape.getEdge(edgeId, tmp);
        if (tmp.isEndpoint(p) && !visitor.test(shape, edgeId, tmp.a, tmp.b)) {
          return false;
        }
      }
    }
    return true;
  }
}
