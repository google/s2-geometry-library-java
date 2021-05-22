/*
 * Copyright 2018 Google Inc.
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
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Multimap;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import com.google.common.geometry.S2ShapeUtil.S2EdgeVectorShape;
import com.google.common.primitives.Ints;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import javax.annotation.Nullable;

/**
 * An encoder/decoder of {@link S2ShapeIndex}s.
 *
 * <p>Values from the {@link S2ShapeIndex} returned by {@link #decode(Bytes, Cursor)} are decoded
 * only when they are accessed. This allows for very fast initialization and no additional memory
 * use beyond the encoded data, and a cache of the clipped shapes that have been accessed. When
 * accessing the entire index, this uses slightly more memory than {@link S2ShapeIndex}, but uses
 * dramatically less memory when accessing only a few cells of the index.
 */
@GwtIncompatible("S2LaxPolylineShape and S2LaxPolygonShape")
public class S2ShapeIndexCoder implements S2Coder<S2ShapeIndex> {

  /**
   * An instance of a {@code S2ShapeIndexCoder} which can encode an {@link S2ShapeIndex} but will
   * throw an {@link IllegalArgumentException} if used to decode an {@link S2ShapeIndex}.
   */
  public static final S2ShapeIndexCoder INSTANCE = new S2ShapeIndexCoder(null);

  private final List<S2Shape> shapes;

  /**
   * Constructs a {@code S2ShapeIndexCoder}.
   *
   * @param shapes the list of shapes, used only by {@link #decode}, commonly the result of {@link
   *     VectorCoder#FAST_SHAPE#decode(Bytes, Cursor)}.
   */
  public S2ShapeIndexCoder(@Nullable List<S2Shape> shapes) {
    this.shapes = shapes;
  }

  @Override
  public void encode(S2ShapeIndex value, OutputStream output) throws IOException {
    // The version number is encoded in 2 bits, under the assumption that by the time we need 5
    // versions the first version can be permanently retired. This only saves 1 byte, but that's
    // significant for very small indexes.
    long maxEdges = value.options().getMaxEdgesPerCell();
    EncodedInts.writeVarint64(output, maxEdges << 2 | EncodedS2ShapeIndex.CURRENT_ENCODING_VERSION);

    List<S2CellId> cellIds = new ArrayList<>();
    List<byte[]> encodedCells = new ArrayList<>();
    Multimap<S2Shape, Integer> shapeIds = S2ShapeUtil.shapeToShapeId(value);
    ByteArrayOutputStream baos = new ByteArrayOutputStream();
    for (S2Iterator<Cell> it = value.iterator(); !it.done(); it.next()) {
      cellIds.add(it.id());
      encodeCell(it.entry(), shapeIds, baos);
      encodedCells.add(baos.toByteArray());
      baos.reset();
    }
    S2CellIdVectorCoder.INSTANCE.encode(cellIds, output);
    VectorCoder.BYTE_ARRAY.encode(encodedCells, output);
  }

  @Override
  public S2ShapeIndex decode(Bytes data, Cursor cursor) {
    Preconditions.checkNotNull(shapes);
    return new EncodedS2ShapeIndex(data, cursor, shapes);
  }

  /**
   * Represents an encoded {@link S2ShapeIndex}.
   *
   * <p>This class is thread-safe.
   */
  private static final class EncodedS2ShapeIndex extends S2ShapeIndex {
    /**
     * Internal representation of an undecoded shape, which must be distinguished from a null shape.
     */
    private static final S2Shape UNDECODED_SHAPE = new S2EdgeVectorShape();

    /** The decoded options of this index. */
    private final Options options;

    /**
     * The array of not-yet-decoded and decoded shapes. The default value is {@link
     * #UNDECODED_SHAPE}. A value of {@code null} represents a null shape.
     */
    private final S2Shape[] cachedShapes;

    /** The encoded vector of cell IDs of this index. */
    private final S2CellIdVector encodedCellIds;

    /** The encoded cells of this index. */
    private final List<S2ClippedShape[]> encodedCells;

    /** The list of {@link S2ShapeIndex.Cell}s. */
    private final List<S2ShapeIndex.Cell> decodedCells;

    /** A coder of {@code S2ClippedShape[]}s. */
    private final S2Coder<S2ClippedShape[]> clippedShapeArrayCoder =
        new S2Coder<S2ClippedShape[]>() {
          @Override
          public void encode(S2ClippedShape[] values, OutputStream output) {
            throw new UnsupportedOperationException();
          }

          @Override
          public S2ClippedShape[] decode(Bytes data, Cursor cursor) {
            return decodeClippedShapes(shapes, data, cursor);
          }
        };

    /**
     * Initializes an {@link EncodedS2ShapeIndex} backed by {@code data} at {@code offset}.
     *
     * <p>Values are decoded only when they are accessed. This allows for very fast initialization
     * and little additional memory use beyond the encoded data.
     */
    EncodedS2ShapeIndex(Bytes data, Cursor cursor, List<S2Shape> shapeFactory) {
      long maxEdgesVersion = data.readVarint64(cursor);
      int version = (int) maxEdgesVersion & 3;
      Preconditions.checkArgument(
          version == S2ShapeIndex.CURRENT_ENCODING_VERSION, "Unknown encoding.");
      options = new Options();
      options.setMaxEdgesPerCell(Ints.checkedCast(maxEdgesVersion >> 2));
      cachedShapes = new S2Shape[Ints.checkedCast(shapeFactory.size())];
      Arrays.fill(cachedShapes, UNDECODED_SHAPE);
      shapes =
          new AbstractList<S2Shape>() {
            @Override
            public synchronized S2Shape get(int i) {
              return cachedShapes[i] == UNDECODED_SHAPE
                  ? cachedShapes[i] = shapeFactory.get(i)
                  : cachedShapes[i];
            }

            @Override
            public int size() {
              return cachedShapes.length;
            }
          };
      encodedCellIds = S2CellIdVectorCoder.INSTANCE.decode(data, cursor);
      encodedCells = new VectorCoder<>(clippedShapeArrayCoder).decode(data, cursor);
      ImmutableList.Builder<S2ShapeIndex.Cell> cellsBuilder = ImmutableList.builder();
      for (int i = 0; i < encodedCellIds.size(); i++) {
        cellsBuilder.add(new LazyCell(i));
      }
      decodedCells = cellsBuilder.build();
    }

    @Override
    public Options options() {
      return options;
    }

    @Override
    public void add(S2Shape shape) {
      throw new UnsupportedOperationException();
    }

    @Override
    public void remove(S2Shape shape) {
      throw new UnsupportedOperationException();
    }

    @Override
    public void reset() {
      throw new UnsupportedOperationException();
    }

    @Override
    public S2Iterator<S2ShapeIndex.Cell> iterator() {
      return S2Iterator.create(decodedCells, encodedCellIds::lowerBound);
    }

    @Override
    public boolean isFresh() {
      return true;
    }

    @Override
    void applyUpdates() {
      throw new UnsupportedOperationException();
    }

    /** A lazy implementation of {@link S2ShapeIndex.Cell} which decodes members on demand. */
    private final class LazyCell extends S2ShapeIndex.Cell {

      /** The index of this cell. */
      private final int i;

      private S2CellId cachedCellId = null;
      private volatile S2ClippedShape[] cachedClippedShapes;

      LazyCell(int i) {
        this.i = i;
      }

      /**
       * Returns {@link #cachedClippedShapes} if it's already cached. Otherwise, loads the clipped
       * shapes from {@link #encodedCells} and stores it in {@link #cachedClippedShapes}.
       */
      private S2ClippedShape[] loadClippedShapesFromCache() {
        if (cachedClippedShapes == null) {
          cachedClippedShapes = encodedCells.get(i);
        }
        return cachedClippedShapes;
      }

      @Override
      public long id() {
        synchronized (EncodedS2ShapeIndex.this) {
          if (cachedCellId == null) {
            cachedCellId = encodedCellIds.get(i);
          }
          return cachedCellId.id();
        }
      }

      @Override
      public int numShapes() {
        return loadClippedShapesFromCache().length;
      }

      @Override
      public S2ClippedShape clipped(int i) {
        return loadClippedShapesFromCache()[i];
      }
    }
  }

  private static void encodeCell(
      Cell cell, Multimap<S2Shape, Integer> shapeIds, OutputStream output) throws IOException {
    // The encoding is designed to be especially compact in certain common situations:
    //
    // 1. The S2ShapeIndex contains exactly one shape.
    //
    // 2. The S2ShapeIndex contains more than one shape, but a particular index cell contains only
    //    one shape (numShapes() == 1).
    //
    // 3. The edge ids for a given shape in a cell form a contiguous range.
    //
    // The details were optimized by constructing an S2ShapeIndex for each feature in Google's
    // geographic repository and measuring their total encoded size. The MutableS2ShapeIndex
    // encoding (of which this function is just one part) uses an average of 1.88 bytes per vertex
    // for features consisting of polygons or polylines.
    //
    // Note that this code does not bother handling numShapes() >= 2**28 or numEdges >= 2**29.
    // This could be fixed using varint64 in a few more places, but if a single cell contains this
    // many shapes or edges then we have bigger problems than just the encoding format :)
    Preconditions.checkArgument(cell.numShapes() < (1 << 28), "Too many shapes.");

    if (shapeIds.size() == 1) {
      S2ClippedShape clipped = cell.clipped(0);
      int n = clipped.numEdges();
      Preconditions.checkArgument(n < (1 << 29), "Too many edges.");

      int containsCenter = clipped.containsCenter() ? 1 : 0;
      if (n >= 2 && n <= 17 && clipped.edge(n - 1) - clipped.edge(0) == n - 1) {
        // The cell contains a contiguous range of edges (*most common case*).
        // If the starting edge id is small then we can encode the cell in one byte. (The n == 0
        // and n == 1 cases are encoded compactly below). This encoding uses a 1-bit tag because
        // it is by far the most common.
        //
        // Encoding: bit 0: 0
        //           bit 1: containsCenter
        //           bits 2-5: (numEdges - 2)
        //           bits 6+: edgeId
        EncodedInts.writeVarint64(
            output, clipped.edge(0) << 6 | (n - 2) << 2 | containsCenter << 1);
      } else if (n == 1) {
        // The cell contains only one edge. For edge ids up to 15, we can encode the cell in a
        // single byte.
        //
        // Encoding: bits 0-1: 1
        //           bit 2: containsCenter
        //           bits 3+: edgeId
        EncodedInts.writeVarint64(output, clipped.edge(0) << 3 | containsCenter << 2 | 1);
      } else {
        // General case (including n == 0, which is encoded compactly here).
        //
        // Encoding: bits 0-1: 3
        //           bit 2: containsCenter
        //           bits 3+: numEdges
        EncodedInts.writeVarint64(output, n << 3 | containsCenter << 2 | 3);
        encodeEdges(clipped, output);
      }
    } else {
      // Note that there are exactly two possible values for the first encoded tag:
      // 1. numShapes() > 1:  a 3-bit tag with value 3.
      // 2. numShapes() == 0: a 3-bit tag with value 7.
      if (cell.numShapes() > 1) {
        // The cell contains more than one shape. The tag for this encoding must be
        // distinguishable from the cases encoded below. We can afford to use a 3-bit tag because
        // numShapes() is generally small.
        EncodedInts.writeVarint64(output, (cell.numShapes() << 3) | 3);
      }
      // The shape ids are delta-encoded.
      int shapeIdBase = 0;
      for (int i = 0; i < cell.numShapes(); i++) {
        S2ClippedShape clipped = cell.clipped(i);
        int containsCenter = clipped.containsCenter() ? 1 : 0;
        int clippedShapeId = -1;
        Iterable<Integer> clippedShapeIds = shapeIds.get(clipped.shape());
        for (int id : clippedShapeIds) {
          if (id >= shapeIdBase) {
            clippedShapeId = id;
            break;
          }
        }
        assert clippedShapeId >= shapeIdBase;
        int shapeDelta = clippedShapeId - shapeIdBase;
        shapeIdBase = clippedShapeId + 1;

        // Like the code above except that we also need to encode shapeId(s).
        // Because of this some choices are slightly different.
        int n = clipped.numEdges();
        Preconditions.checkArgument(n < (1 << 29), "Too many edges.");

        if (n >= 1 && n <= 16 && clipped.edge(n - 1) - clipped.edge(0) == n - 1) {
          // The clipped shape has a contiguous range of up to 16 edges. This encoding uses a
          // 1-bit tag because it is by far the most common.
          //
          // Encoding: bit 0: 0
          //           bit 1: containsCenter
          //           bits 2+: edgeId
          // Next value: bits 0-3: (numEdges - 1)
          //             bits 4+: shapeDelta
          EncodedInts.writeVarint64(output, (clipped.edge(0) << 2) | (containsCenter << 1));
          EncodedInts.writeVarint64(output, (shapeDelta << 4) | (n - 1));
        } else if (n == 0) {
          // Special encoding for clipped shapes with no edges. Such shapes are common in polygon
          // interiors. This encoding uses a 3-bit tag in order to leave more bits available for
          // the other encodings.
          //
          // NOTE(user): When numShapes() > 1, this tag could be 2 bits (because the tag used to
          // indicate numShapes() > 1 can't appear). Alternatively, that tag can be considered
          // reserved for future use.
          //
          // Encoding: bits 0-2: 7
          //           bit 3: containsCenter
          //           bits 4+: shapeDelta
          EncodedInts.writeVarint64(output, (shapeDelta << 4) | (containsCenter << 3) | 7);
        } else {
          // General case. This encoding uses a 2-bit tag, and the first value typically is
          // encoded into one byte.
          //
          // Encoding: bits 0-1: 1
          //           bit 2: containsCenter
          //           bits 3+: (numEdges - 1)
          // Next value: shapeDelta
          EncodedInts.writeVarint64(output, ((n - 1) << 3) | (containsCenter << 2) | 1);
          EncodedInts.writeVarint64(output, shapeDelta);
          encodeEdges(clipped, output);
        }
      }
    }
  }

  /** Encodes the edge IDs of the given {@link S2ClippedShape}. */
  private static void encodeEdges(S2ClippedShape clipped, OutputStream output) throws IOException {
    // Each entry is an (edgeId, count) pair representing a contiguous range of edges. The edge
    // ids are delta-encoded such that 0 represents the minimum valid next edge id.
    //
    // Encoding: if bits 0-2 < 7: encodes (count - 1)
    //            - bits 3+: edge delta
    //           if bits 0-2 == 7:
    //            - bits 3+ encode (count - 8)
    //            - Next value is edge delta
    //
    // No count is encoded for the last edge (saving 3 bits).
    int edgeIdBase = 0;
    int numEdges = clipped.numEdges();
    for (int i = 0; i < numEdges; i++) {
      int edgeId = clipped.edge(i);
      assert edgeId >= edgeIdBase;
      int delta = edgeId - edgeIdBase;
      if (i + 1 == numEdges) {
        // This is the last edge; no need to encode an edge count.
        EncodedInts.writeVarint64(output, delta);
      } else {
        // Count the edges in this contiguous range.
        int count = 1;
        for (; i + 1 < numEdges && clipped.edge(i + 1) == edgeId + count; i++) {
          count++;
        }
        if (count < 8) {
          // Count is encoded in low 3 bits of delta.
          EncodedInts.writeVarint64(output, delta << 3 | (count - 1));
        } else {
          // Count and delta are encoded separately.
          EncodedInts.writeVarint64(output, (count - 8) << 3 | 7);
          EncodedInts.writeVarint64(output, delta);
        }
        edgeIdBase = edgeId + count;
      }
    }
  }

  /** Decodes {@code numEdges} edge IDs of a {@link S2ClippedShape}. */
  private static int[] decodeEdges(int numEdges, Bytes data, Cursor cursor) {
    // This function inverts the encodings documented above.
    int[] edges = new int[numEdges];
    int edgeId = 0;
    for (int i = 0; i < numEdges; ) {
      long delta = data.readVarint64(cursor);
      if (i + 1 == numEdges) {
        // The last edge is encoded without an edge count.
        edges[i++] = Ints.checkedCast(edgeId + delta);
      } else {
        // Otherwise decode the count and edge delta.
        long count = (delta & 7) + 1;
        delta >>>= 3;
        if (count == 8) {
          count = delta + 8;
          delta = data.readVarint64(cursor);
        }
        edgeId += Ints.checkedCast(delta);
        for (; count > 0; count--, i++, edgeId++) {
          edges[i] = edgeId;
        }
      }
    }
    return edges;
  }

  /**
   * Decodes an array of {@link S2ClippedShape} from {@code input} from the given {@code shapes}.
   * The {@link S2ClippedShape} at index 0 will store {@code cellId}.
   */
  private static S2ClippedShape[] decodeClippedShapes(
      List<S2Shape> shapes, Bytes data, Cursor cursor) {
    // This function inverts the encodings documented above.
    if (shapes.size() == 1) {
      S2ClippedShape[] clippedShapes = new S2ClippedShape[1];

      // Entire S2ShapeIndex contains only one shape.
      long header = data.readVarint64(cursor);
      if ((header & 1) == 0) {
        // The cell contains a contiguous range of edges.
        int numEdges = Ints.checkedCast(((header >>> 2) & 15) + 2);
        clippedShapes[0] =
            S2ClippedShape.create(
                null, shapes.get(0), (header & 2) != 0, Ints.checkedCast(header >>> 6), numEdges);
      } else if ((header & 2) == 0) {
        // The cell contains a single edge.
        clippedShapes[0] =
            S2ClippedShape.create(
                null, shapes.get(0), (header & 4) != 0, Ints.checkedCast(header >>> 3), 1);
      } else {
        // The cell contains some other combination of edges.
        int numEdges = Ints.checkedCast(header >> 3);
        int[] edges = decodeEdges(numEdges, data, cursor);
        clippedShapes[0] = S2ClippedShape.create(null, shapes.get(0), (header & 4) != 0, edges);
      }
      return clippedShapes;
    }

    // S2ShapeIndex contains more than one shape.
    long header = data.readVarint64(cursor);
    int numClipped = 1;
    if ((header & 7) == 3) {
      // This cell contains more than one shape.
      numClipped = Ints.checkedCast(header >>> 3);
      header = data.readVarint64(cursor);
    }

    S2ClippedShape[] clippedShapes = new S2ClippedShape[numClipped];

    long shapeId = 0;
    for (int j = 0; j < numClipped; j++, shapeId++) {
      if (j > 0) {
        header = data.readVarint64(cursor);
      }
      if ((header & 1) == 0) {
        // The clipped shape contains a contiguous range of edges.
        long shapeIdCount = data.readVarint64(cursor);
        shapeId += shapeIdCount >> 4;
        int numEdges = Ints.checkedCast((shapeIdCount & 15) + 1);
        clippedShapes[j] =
            S2ClippedShape.create(
                null,
                shapes.get(Ints.checkedCast(shapeId)),
                (header & 2) != 0,
                Ints.checkedCast(header >>> 2),
                numEdges);
      } else if ((header & 7) == 7) {
        // The clipped shape has no edges.
        shapeId += header >> 4;
        clippedShapes[j] =
            S2ClippedShape.create(
                null, shapes.get(Ints.checkedCast(shapeId)), (header & 8) != 0, 0, 0);
      } else {
        // The clipped shape contains some other combination of edges.
        assert (header & 3) == 1;
        long shapeDelta = data.readVarint64(cursor);
        shapeId += shapeDelta;
        int numEdges = Ints.checkedCast((header >>> 3) + 1);
        int[] edges = decodeEdges(numEdges, data, cursor);
        clippedShapes[j] =
            S2ClippedShape.create(
                null, shapes.get(Ints.checkedCast(shapeId)), (header & 4) != 0, edges);
      }
    }
    return clippedShapes;
  }
}
