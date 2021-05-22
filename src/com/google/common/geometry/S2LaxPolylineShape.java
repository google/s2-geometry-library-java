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

import com.google.common.annotations.GwtIncompatible;
import com.google.common.base.Predicates;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.S2ShapeAspect.ChainAspect;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

/**
 * S2LaxPolylineShape represents a polyline. It is similar to {@link S2Polyline} except that
 * duplicate vertices are allowed, and the representation is slightly more compact since this class
 * does not implement {@link S2Region}.
 *
 * <p>Polylines may have any number of vertices, but note that polylines with fewer than 2 vertices
 * do not define any edges. To create a polyline consisting of a single degenerate edge, repeat the
 * same vertex twice.
 */
@GwtIncompatible("S2ShapeAspect incompatible")
public interface S2LaxPolylineShape extends S2ShapeAspect.EdgeAspect.Open {
  /** A polyline with no edges. */
  static final S2LaxPolylineShape EMPTY = new SimpleArray(ImmutableList.of());

  /** Creates a lax polyline from the {@code line} by copying its data. */
  public static S2LaxPolylineShape create(S2Polyline line) {
    return create(line.vertices());
  }

  /** Creates a new lax polyline from the given vertices. */
  public static S2LaxPolylineShape create(Iterable<S2Point> vertices) {
    vertices = filterLine(vertices);
    return Iterables.isEmpty(vertices) ? EMPTY : new SimpleArray(vertices);
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
    } else if (Iterables.size(lines) == 1) {
      return new SimpleArray(Iterables.getOnlyElement(lines));
    } else {
      return new MultiArray(lines);
    }
  }
  
  /** As {@link #create}, but with coordinates packed into a double[]. */
  public static S2LaxPolylineShape createMultiPacked(Iterable<? extends Iterable<S2Point>> lines) {
    lines = filterLines(lines);
    if (Iterables.isEmpty(lines)) {
      return EMPTY;
    } else if (Iterables.size(lines) == 1) {
      return new SimplePacked(Iterables.getOnlyElement(lines));
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
    } else if (Iterables.size(lines) == 1) {
      return new SimpleSnapped(Iterables.getOnlyElement(lines));
    } else {
      return new MultiSnapped(lines);
    }
  }
  
  /** Returns 'input' or an empty iterable if 'input' has only one vertex. */
  static <T> Iterable<T> filterLine(Iterable<T> input) {
    return Iterables.size(input) < 2 ? ImmutableList.of() : input;
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
  
  @Override default int dimension() {
    return 1;
  }
  
  @Override default boolean hasInterior() {
    return false;
  }
  
  @Override default boolean containsOrigin() {
    return false;
  }

  @Override
  default int numEdges() {
    return numVertices() == 0 ? 0 : numVertices() - numChains();
  }

  /** Returns true unless there is at least one edge in this line. */
  default boolean isEmpty() {
    return numEdges() == 0;
  }
  
  /** Returns false in all cases since a polyline may never cover the entire sphere. */
  default boolean isFull() {
    return false;
  }
  
  /** A polyline storing references to previously allocated S2Point instances. */
  static class SimpleArray extends ChainAspect.Simple.Array implements S2LaxPolylineShape {
    private SimpleArray(Iterable<S2Point> vertices) {
      super(vertices);
    }
  }

  /** A polyline storing {@link S2Point}s in a {@link List<S2Point>}. */
  static class SimpleList extends ChainAspect.Simple implements S2LaxPolylineShape {
    private final List<S2Point> vertices;

    private SimpleList(List<S2Point> vertices) {
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

  /** An encoder/decoder of {@link S2LaxPolylineShape}s. */
  @GwtIncompatible("Uses EncodedS2PointVector")
  class Coder implements S2Coder<S2LaxPolylineShape> {

    /**
     * An instance of {@link Coder} which encodes/decodes {@link S2LaxPolylineShape}s in the {@code
     * FAST} format.
     */
    static final Coder FAST = new Coder(S2PointVectorCoder.FAST);
    /**
     * An instance of {@link Coder} which encodes/decodes {@link S2LaxPolylineShape}s in the {@code
     * COMPACT} format.
     */
    static final Coder COMPACT = new Coder(S2PointVectorCoder.COMPACT);

    private final S2Coder<List<S2Point>> coder;

    private Coder(S2Coder<List<S2Point>> coder) {
      this.coder = coder;
    }

    @Override
    public void encode(S2LaxPolylineShape shape, OutputStream output) throws IOException {
      coder.encode(shape.vertices(), output);
    }

    @Override
    public S2LaxPolylineShape decode(Bytes data, Cursor cursor) {
      return new SimpleList(coder.decode(data, cursor));
    }
  }
}
