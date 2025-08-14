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

import com.google.common.collect.ImmutableList;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.PrimitiveArrays.Longs;
import com.google.common.geometry.S2ShapeAspect.ChainAspect;
import com.google.common.primitives.ImmutableLongArray;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

/**
 * A region defined by a collection of zero or more closed loops. The interior is the region to the
 * left of all loops. Loops are not closed, that is, the last edge is an implicit path from the last
 * vertex back to the first vertex.
 *
 * <p>This is similar to {@link S2Polygon#shape}, except that this class supports polygons with two
 * types of degeneracy:
 *
 * <ol>
 *   <li>Degenerate edges (from a vertex to itself)
 *   <li>Sibling edge pairs (consisting of two oppositely oriented edges)
 * </ol>
 *
 * <p>Degeneracies can represent either "shells" or "holes" depending on the loop they are contained
 * by. For example, a degenerate edge or sibling pair contained by a "shell" would be interpreted as
 * a degenerate hole. Such edges form part of the boundary of the polygon.
 *
 * <p>Loops with fewer than three vertices are interpreted as follows:
 *
 * <ul>
 *   <li>A loop with two vertices defines two edges (in opposite directions).
 *   <li>A loop with one vertex defines a single degenerate edge.
 *   <li>A loop with no vertices is interpreted as the "full loop" containing all points on the
 *       sphere. If this loop is present, then all other loops must form degeneracies (i.e.,
 *       degenerate edges or sibling pairs). For example, two loops {} and {X} would be interpreted
 *       as the full polygon with a degenerate single-point hole at X.
 * </ul>
 *
 * <p>No error checking is performed during construction. It is perfectly fine to create objects
 * that do not meet the requirements below (e.g., in order to analyze or fix those problems).
 * However, some additional conditions must be satisfied in order to perform certain operations:
 *
 * <ul>
 *   <li>In order to be valid for point containment tests, the polygon must satisfy the "interior is
 *       on the left" rule. This means that there must not be any crossing edges, and if there are
 *       duplicate edges then all but at most one of them must belong to a sibling pair (i.e., the
 *       number of edges in opposite directions must differ by at most one).
 *   <li>To be valid for boolean operations, degenerate edges and sibling pairs cannot coincide with
 *       any other edges. For example, the following situations are not allowed:
 *       <ul>
 *         <li>{AA, AA} // degenerate edge coincides with another edge
 *         <li>{AA, AB} // degenerate edge coincides with another edge
 *         <li>{AB, BA, AB} // sibling pair coincides with another edge
 *       </ul>
 * </ul>
 *
 * <p>Note that this class is faster to initialize and is more compact than {@link S2Polygon}, but
 * it does not have any built-in operations, as those are by design provided by other classes. All
 * the design considerations here are focused on the meaning and storage of the model itself. All
 * implementations store a single list of vertices, and if there are multiple loops, an int[] that
 * provides the offset into the list where each loop's vertices begin. This scales at a rate of one
 * int per loop. This compares favorably to {@link S2Polygon}, which requires 4 objects and an int
 * per loop. Heap size of the vertex data scales at different rates, depending on the storage:
 *
 * <ol>
 *   <li>{@link #create Standard} polygons copy points into one {@link ImmutableList immutable
 *       list}, which requires 48 bytes per vertex on a standard 64-bit JVM.
 *   <li>{@link #createPacked Packed} polygons copy coordinates into a double[], and convert these
 *       coordinates back to {@link S2Point point} instances on demand. This consumes 24
 *       bytes/vertex, but construction and other operations are slower and drive the garbage
 *       collector harder.
 *   <li>{@link #createSnapped Snapped} polygons copy {@link S2CellId cells} into a long[], and
 *       convert the cell centers to {@link S2Point point} instances on demand. This consumes just 8
 *       bytes/vertex, but construction and especially operations are even slower and drive the
 *       garbage collector even harder.
 * </ol>
 *
 * <p>S2LaxPolygonShape implementations are expected to be immutable.
 */
public interface S2LaxPolygonShape extends S2ShapeAspect.EdgeAspect.Closed {
  // When adding a new encoding, be aware that old binaries will not be able to decode it.
  static final byte CURRENT_ENCODING_VERSION = 1;

  /** A singleton for the empty polygon. */
  public static S2LaxPolygonShape EMPTY = new MultiArray(ImmutableList.of());

  /** A singleton for the full polygon. */
  public static S2LaxPolygonShape FULL = new SimpleArray(ImmutableList.of());

  /**
   * An instance of {@link Coder} which encodes/decodes {@link S2LaxPolygonShape}s in the {@code
   * FAST} format.
   */
  public S2Coder<S2LaxPolygonShape> FAST_CODER = new Coder(S2PointVectorCoder.FAST);

  /**
   * An instance of {@link Coder} which encodes/decodes {@link S2LaxPolygonShape}s in the {@code
   * COMPACT} format.
   */
  public S2Coder<S2LaxPolygonShape> COMPACT_CODER = new Coder(S2PointVectorCoder.COMPACT);

  /** Creates a polygon from the given {@link S2Polygon} by copying its data. */
  public static S2LaxPolygonShape create(S2Polygon polygon) {
    if (polygon.isEmpty()) {
      return EMPTY;
    } else if (polygon.isFull()) {
      return FULL;
    } else {
      // S2Polygon filters out empty loops already. Convert full loops to empty lists.
      // Other loops must simply be oriented.
      return create(
          Lists.transform(
              polygon.getLoops(), x -> x.isFull() ? ImmutableList.of() : x.orientedVertices()));
    }
  }

  /**
   * Creates a polygon from a single loop, copying the vertices to ensure the resulting polygon is
   * deeply immutable. If the given loop is empty, the full polygon is the result. Otherwise the
   * resulting polygon's interior is on the left of the loop when walking the vertices in the given
   * order.
   *
   * <p>The loop should not be closed, that is, the last vertex in each inner iterable should differ
   * from the first vertex, since an implicit edge from the last vertex back to the first is
   * assumed.
   */
  public static S2LaxPolygonShape fromLoop(Iterable<S2Point> loop) {
    if (!loop.iterator().hasNext()) {
      return FULL;
    }
    return new SimpleArray(loop);
  }

  /**
   * Creates a polygon from the given loops, defensively copying any loop's Iterable except an
   * {@link ImmutableList}, to ensure the polygon is deeply immutable.
   *
   * <p>If given no loops, the empty polygon is the result. If given only empty loops, the full
   * polygon is the result. Otherwise the resulting polygon's interior is on the left of the loops
   * when walking the vertices in the given order.
   *
   * <p>Each loop should not be closed, that is, the last vertex in each inner iterable should
   * differ from the first vertex, since an implicit edge from the last vertex back to the first is
   * assumed.
   */
  public static S2LaxPolygonShape create(Iterable<? extends Iterable<S2Point>> loops) {
    if (Iterables.isEmpty(loops)) {
      return EMPTY;
    } else if (Iterables.all(loops, Iterables::isEmpty)) {
      return FULL;
    } else if (Iterables.size(loops) == 1) {
      return new SimpleArray(getOnlyElement(loops));
    } else {
      return new MultiArray(loops);
    }
  }

  /**
   * As {@link #create}, but packs coordinates into a double[] array. Operations are slower since
   * S2Points are constructed on each access, but this representation has vastly fewer objects, and
   * so can be a better choice if polygons may be held in RAM for a long time.
   */
  public static S2LaxPolygonShape createPacked(Iterable<? extends Iterable<S2Point>> loops) {
    if (Iterables.isEmpty(loops)) {
      return EMPTY;
    } else if (Iterables.all(loops, Iterables::isEmpty)) {
      return FULL;
    } else if (Iterables.size(loops) == 1) {
      return new SimplePacked(getOnlyElement(loops));
    } else {
      return new MultiPacked(loops);
    }
  }

  /**
   * As {@link #create}, but packs vertices into a long[] array. Operations may be much slower since
   * S2Points are constructed on each access, but this representation is the smallest, and so may be
   * far better if polygons may be held in RAM for a long time.
   */
  public static S2LaxPolygonShape createSnapped(Iterable<? extends Iterable<S2CellId>> loops) {
    if (Iterables.isEmpty(loops)) {
      return EMPTY;
    } else if (Iterables.all(loops, Iterables::isEmpty)) {
      return FULL;
    } else if (Iterables.size(loops) == 1) {
      return new SimpleSnapped(getOnlyElement(loops));
    } else {
      return new MultiSnapped(loops);
    }
  }

  /** Canonicalizes the empty/full instances on deserialization. */
  default Object readResolve() {
    int n = numChains();
    if (n == 0) {
      return EMPTY;
    }
    for (int i = 0; i < n; i++) {
      if (getChainLength(i) != 0) {
        return this;
      }
    }
    return FULL;
  }

  @Override
  default int dimension() {
    return 2;
  }

  /** Returns true if this polygon contains no area, i.e. has no loops. */
  @Override
  default boolean isEmpty() {
    return numChains() == 0;
  }

  /** Returns true if this polygon contains all points, i.e. there are loops, but all are empty. */
  @Override
  default boolean isFull() {
    int n = numChains();
    for (int i = 0; i < n; i++) {
      if (getChainLength(i) != 0) {
        return false;
      }
    }
    return n > 0;
  }

  @Override
  default boolean hasInterior() {
    return true;
  }

  @Override
  default boolean containsOrigin() {
    if (isFull()) {
      return true;
    } else if (isEmpty()) {
      return false;
    } else {
      return S2ShapeUtil.containsBruteForce(this, S2.origin());
    }
  }

  @Override
  default ReferencePoint getReferencePoint() {
    return S2ShapeUtil.getReferencePoint(this);
  }

  /** A simple polygon with points referenced from an array. */
  static class SimpleArray extends ChainAspect.Simple.Array implements S2LaxPolygonShape {
    private SimpleArray(Iterable<S2Point> vertices) {
      super(vertices);
    }
  }

  /** A simple polygon with vertices referenced from a {@link List} of {@link S2Point}. */
  static class SimpleList extends ChainAspect.Simple implements S2LaxPolygonShape {
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

  /** A simple polygon with vertex coordinates stored in a double[]. */
  static class SimplePacked extends ChainAspect.Simple.Packed implements S2LaxPolygonShape {
    SimplePacked(Iterable<S2Point> vertices) {
      super(vertices);
    }
  }

  /** A simple polygon with vertices at cell ID centers stored in a long[]. */
  static class SimpleSnapped extends ChainAspect.Simple.Snapped implements S2LaxPolygonShape {
    SimpleSnapped(Iterable<S2CellId> vertices) {
      super(vertices);
    }
  }

  /** A multi polygon with points referenced from an array. */
  static class MultiArray extends ChainAspect.Multi.Array implements S2LaxPolygonShape {
    private MultiArray(Iterable<? extends Iterable<S2Point>> loops) {
      super(loops);
    }
  }

  /**
   * A multi polygon with vertices referenced from a {@link List} of {@link S2Point}, and cumulative
   * edges referenced from an {@link Longs}.
   */
  static class MultiList extends ChainAspect.Multi implements S2LaxPolygonShape {
    private final List<S2Point> vertices;

    // Note this will throw an IllegalArgumentException if any of the Long cumulativeEdges cannot
    // be represented as an int.
    private MultiList(List<S2Point> vertices, Longs cumulativeEdges) {
      super(cumulativeEdges.toIntArray());
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

  /** A multi polygon with vertex coordinates stored in a double[]. */
  static class MultiPacked extends ChainAspect.Multi.Packed implements S2LaxPolygonShape {
    MultiPacked(Iterable<? extends Iterable<S2Point>> loops) {
      super(loops);
    }
  }

  /** A multi polygon with vertices at cell ID centers stored in a long[]. */
  static class MultiSnapped extends ChainAspect.Multi.Snapped implements S2LaxPolygonShape {
    MultiSnapped(Iterable<? extends Iterable<S2CellId>> loops) {
      super(loops);
    }
  }

  /** An encoder/decoder of {@link S2LaxPolygonShape}s. */
  class Coder implements S2Coder<S2LaxPolygonShape> {
    private final S2Coder<List<S2Point>> coder;

    private Coder(S2Coder<List<S2Point>> coder) {
      this.coder = coder;
    }

    @Override
    public void encode(S2LaxPolygonShape shape, OutputStream output) throws IOException {
      output.write(CURRENT_ENCODING_VERSION);
      // Write the number of loops.
      EncodedInts.writeVarint64(output, shape.numChains());
      coder.encode(shape.vertices(), output);
      if (shape.numChains() > 1) {
        ImmutableLongArray.Builder builder = ImmutableLongArray.builder();
        for (int i = 0; i < shape.numChains(); i++) {
          builder.add(shape.getChainStart(i));
        }
        builder.add(shape.numVertices());
        UintVectorCoder.UINT32.encode(Longs.fromImmutableLongArray(builder.build()), output);
      }
    }

    @Override
    public S2LaxPolygonShape decode(Bytes data, Cursor cursor) throws IOException {
      long numChains;
      try {
        // Bytes.get throws IndexOutOfBoundsException if the data is short.
        byte version = data.get(cursor.position++);
        if (version != CURRENT_ENCODING_VERSION) {
          throw new IOException(
              Platform.formatString(
                  "Expected encoding version %s, got %s.", CURRENT_ENCODING_VERSION, version));
        }
        // Bytes.readVarint64 throws IllegalArgumentException if the varint is malformed.
        numChains = Math.toIntExact(data.readVarint64(cursor));
      } catch (IndexOutOfBoundsException | IllegalArgumentException e) {
        throw new IOException("Insufficient or invalid input bytes: ", e);
      }

      // Both FAST and COMPACT coders are capable of decoding any encoding format. Both decode
      // points lazily (on-demand).
      List<S2Point> vertices = S2PointVectorCoder.FAST.decode(data, cursor);
      if (numChains == 0) {
        return S2LaxPolygonShape.EMPTY;
      } else if (numChains == 1) {
        return new SimpleList(vertices);
      }
      Longs cumulativeEdges = UintVectorCoder.UINT32.decode(data, cursor);
      if (cumulativeEdges.length() < 0) {
        throw new IOException("Invalid cumulative edges length: " + cumulativeEdges.length());
      }
      try {
        return new MultiList(vertices, cumulativeEdges);
      } catch (IllegalArgumentException e) {
        throw new IOException("Invalid input data: ", e);
      }
    }

    @Override
    public boolean isLazy() {
      return true;
    }
  }
}
