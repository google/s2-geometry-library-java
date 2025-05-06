/*
 * Copyright 2019 Google Inc.
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

import static com.google.common.collect.Iterables.getOnlyElement;
import static com.google.common.collect.Iterables.size;

import com.google.common.base.Predicates;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.geometry.S2ShapeAspect.ChainAspect;
import java.io.IOException;
import java.util.List;

/**
 * S2LaxPolylineShape represents a polyline or a collection of polylines (a multipolyline). When
 * representing a single polyline, it is similar to {@link S2Polyline}, except that consecutive
 * duplicate vertices are allowed, and the representation is slightly more compact since this class
 * does not implement {@link S2Region}.
 *
 * <p>Polylines with fewer than two vertices do not define any edges, and attempting to create an
 * S2LaxPolylineShape with a single vertex discards the vertex, resulting in an empty polyline.
 * Attempting to decode an S2LaxPolylineShape with a single vertex throws an IOException. To create
 * an S2LaxPolylineShape representing a degenerate edge, use two identical vertices.
 *
 * <p>WARNING: Multipolylines are only partially implemented. In particular, they may not be
 * encoded.
 */
public interface S2LaxPolylineShape extends S2ShapeAspect.EdgeAspect.Open {
  // TODO(user): Complete the documentation and implementation of multipolylines.

  /** A polyline with no edges. */
  static final S2LaxPolylineShape EMPTY = new SimpleArray(ImmutableList.of());

  /**
   * An instance of {@link S2Coder} which encodes/decodes {@link S2LaxPolylineShape}s in the {@code
   * FAST} format.
   */
  public S2Coder<S2LaxPolylineShape> FAST_CODER =
      S2PointVectorCoder.FAST.delegating(S2LaxPolylineShape::vertices, SimpleList::new);

  /**
   * An instance of {@link S2Coder} which encodes/decodes {@link S2LaxPolylineShape}s in the {@code
   * COMPACT} format.
   */
  public S2Coder<S2LaxPolylineShape> COMPACT_CODER =
      S2PointVectorCoder.COMPACT.delegating(S2LaxPolylineShape::vertices, SimpleList::new);

  /**
   * Creates a lax polyline from the {@link S2Polyline} by copying its data. Single-vertex
   * S2Polylines, which represent a degenerate edge, are converted to the S2LaxPolylineShape
   * representation of a degenerate edge, using two identical vertices.
   */
  public static S2LaxPolylineShape create(S2Polyline line) {
    if (line.numVertices() == 1) {
      S2Point[] pts = new S2Point[] { line.vertices().get(0), line.vertices().get(0) };
      return create(pts);
    }
    return create(line.vertices());
  }

  /**
   * Creates a new lax polyline from the given vertices. Input consisting of zero or one vertex
   * produces an empty line.
   */
  public static S2LaxPolylineShape create(Iterable<S2Point> vertices) {
    vertices = filterLine(vertices);
    return Iterables.isEmpty(vertices) ? EMPTY : new SimpleArray(vertices);
  }

  /**
   * Creates a new lax polyline from the given vertices. Takes ownership of the array of vertices,
   * which must not be modified after this call. Package private due to the ownership transfer.
   * Input consisting of zero or one vertex produces an empty line.
   */
  static S2LaxPolylineShape create(S2Point[] vertices) {
    return (vertices.length < 2) ? EMPTY : new SimpleArray(vertices);
  }

  /** As {@link #create}, but with coordinates packed into a double[]. */
  public static S2LaxPolylineShape createPacked(Iterable<S2Point> vertices) {
    vertices = filterLine(vertices);
    return Iterables.isEmpty(vertices) ? EMPTY : new SimplePacked(vertices);
  }

  /** As {@link #create}, but with vertices at the center of cell IDs, packed into a long[]. */
  public static S2LaxPolylineShape createSnapped(Iterable<S2CellId> vertices) {
    vertices = filterLine(vertices);
    return Iterables.isEmpty(vertices) ? EMPTY : new SimpleSnapped(vertices);
  }

  /** Creates a new lax multipolyline with the given lines. */
  public static S2LaxPolylineShape createMulti(Iterable<? extends Iterable<S2Point>> lines) {
    lines = filterLines(lines);
    if (Iterables.isEmpty(lines)) {
      return EMPTY;
    } else if (size(lines) == 1) {
      return new SimpleArray(getOnlyElement(lines));
    } else {
      return new MultiArray(lines);
    }
  }

  /** As {@link #create}, but with coordinates packed into a double[]. */
  public static S2LaxPolylineShape createMultiPacked(Iterable<? extends Iterable<S2Point>> lines) {
    lines = filterLines(lines);
    if (Iterables.isEmpty(lines)) {
      return EMPTY;
    } else if (size(lines) == 1) {
      return new SimplePacked(getOnlyElement(lines));
    } else {
      return new MultiPacked(lines);
    }
  }

  /** As {@link #create}, but with vertices at the center of cell IDs, packed into a long[]. */
  public static S2LaxPolylineShape createMultiSnapped(
      Iterable<? extends Iterable<S2CellId>> lines) {
    lines = filterLines(lines);
    if (Iterables.isEmpty(lines)) {
      return EMPTY;
    } else if (size(lines) == 1) {
      return new SimpleSnapped(getOnlyElement(lines));
    } else {
      return new MultiSnapped(lines);
    }
  }

  /** Returns 'input' or an empty iterable if 'input' has only one vertex. */
  static <T> Iterable<T> filterLine(Iterable<T> input) {
    return size(input) < 2 ? ImmutableList.of() : input;
  }

  static <T> Iterable<? extends Iterable<T>> filterLines(Iterable<? extends Iterable<T>> input) {
    return Iterables.filter(
        Iterables.transform(input, S2LaxPolylineShape::filterLine),
        Predicates.not(Iterables::isEmpty));
  }

  /** Canonicalize exactly empty polylines to EMPTY. */
  default Object readResolve() {
    return numVertices() == 0 ? EMPTY : this;
  }

  @Override
  default int dimension() {
    return 1;
  }

  @Override
  default boolean hasInterior() {
    return false;
  }

  @Override
  default boolean containsOrigin() {
    return false;
  }

  @Override
  default int numEdges() {
    return numVertices() == 0 ? 0 : numVertices() - numChains();
  }

  /** Returns true unless there is at least one edge in this line. */
  @Override
  default boolean isEmpty() {
    return numEdges() == 0;
  }

  /** Returns false in all cases since a polyline may never cover the entire sphere. */
  @Override
  default boolean isFull() {
    return false;
  }

  /** A polyline storing references to previously allocated S2Point instances. */
  static class SimpleArray extends ChainAspect.Simple.Array implements S2LaxPolylineShape {
    /** This constructor copies the provided Iterable into an array. */
    private SimpleArray(Iterable<S2Point> vertices) {
      super(vertices);
    }

    /** This constructor takes ownership of the provided array. */
    private SimpleArray(S2Point[] vertices) {
      super(vertices);
    }
  }

  /**
   * A polyline abstract class that requires implementation to define {@code numVertices} and {@code
   * vertex}.
   */
  abstract static class Simple extends ChainAspect.Simple implements S2LaxPolylineShape {}

  /** A polyline storing {@link S2Point}s in a {@link List} of {@link S2Point}. */
  static class SimpleList extends ChainAspect.Simple implements S2LaxPolylineShape {
    private final List<S2Point> vertices;

    /**
     * Constructs a SimpleList wrapping the given List of vertices. As described above, it is not
     * possible to construct a single-vertex S2LaxPolylineShape, so if such a list is encountered
     * here is is likely due to corrupt input to the S2LaxPolylineShape decoder.
     */
    private SimpleList(List<S2Point> vertices) throws IOException {
      if (vertices.size() == 1) {
        throw new IOException(
            "A valid S2LaxPolylineShape may not have a single vertex.");
      }
      this.vertices = vertices;
    }

    @Override
    public int numVertices() {
      return vertices.size();
    }

    @Override
    public S2Point vertex(int vertexId) {
      return vertices.get(vertexId);
    }
  }

  /** A polyline storing xyz coordinates in a single packed 'double' array. */
  static class SimplePacked extends ChainAspect.Simple.Packed implements S2LaxPolylineShape {
    private SimplePacked(Iterable<S2Point> vertices) {
      super(vertices);
    }
  }

  /** A polyline storing cell IDs in a single 'long' array. */
  static class SimpleSnapped extends ChainAspect.Simple.Snapped implements S2LaxPolylineShape {
    private SimpleSnapped(Iterable<S2CellId> vertices) {
      super(vertices);
    }
  }

  /** A multi polyline storing references to previously allocated S2Point instances. */
  static class MultiArray extends ChainAspect.Multi.Array implements S2LaxPolylineShape {
    private MultiArray(Iterable<? extends Iterable<S2Point>> chains) {
      super(chains);
    }
  }

  /** A multi polyline storing xyz coordinates in a single packed 'double' array. */
  static class MultiPacked extends ChainAspect.Multi.Packed implements S2LaxPolylineShape {
    MultiPacked(Iterable<? extends Iterable<S2Point>> chains) {
      super(chains);
    }
  }

  /** A multi polyline storing cell IDs in a single 'long' array. */
  static class MultiSnapped extends ChainAspect.Multi.Snapped implements S2LaxPolylineShape {
    MultiSnapped(Iterable<? extends Iterable<S2CellId>> chains) {
      super(chains);
    }
  }
}
