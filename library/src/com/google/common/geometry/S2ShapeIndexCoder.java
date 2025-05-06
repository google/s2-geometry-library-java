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

import static com.google.common.geometry.EncodedInts.writeVarint64;
import static java.lang.Math.addExact;

import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.S2Iterator.ListIterator;
import com.google.common.geometry.S2ShapeIndex.Cell;
import com.google.common.geometry.S2ShapeIndex.Options;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape;
import com.google.common.geometry.S2ShapeIndex.S2ClippedShape.ManyEdges;
import com.google.common.geometry.S2ShapeUtil.S2EdgeVectorShape;
import com.google.common.primitives.Ints;
import com.google.common.primitives.UnsignedLongs;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import jsinterop.annotations.JsConstructor;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;
import org.jspecify.annotations.Nullable;

/**
 * An encoder/decoder of {@link S2ShapeIndex}s. Subclasses may override the cache methods to decide
 * when to decode shapes or clipped shapes. By default, values from the {@link S2ShapeIndex}
 * returned by {@link #decode(Bytes, Cursor)} are decoded only when they are accessed. This allows
 * for very fast initialization and no additional memory use beyond the encoded data, and a cache of
 * the clipped shapes that have been accessed. When accessing the entire index, this uses slightly
 * more memory than {@link S2ShapeIndex}, but uses dramatically less memory when accessing only a
 * few cells of the index.
 */
@JsType
@SuppressWarnings("Assertion")
public class S2ShapeIndexCoder implements S2Coder<S2ShapeIndex> {

  /**
   * An instance of a {@code S2ShapeIndexCoder} which can encode an {@link S2ShapeIndex} but will
   * throw an {@link IllegalArgumentException} if used to decode an {@link S2ShapeIndex}.
   */
  public static final S2ShapeIndexCoder INSTANCE = new S2ShapeIndexCoder(null);

  private final List<S2Shape> shapes;

  /**
   * Constructs a {@code S2ShapeIndexCoder}. If the optional {@code shapes} list is non-null, it
   * must contain the shapes that are in the index to be decoded, and this S2ShapeIndexCoder will
   * support decoding. If {@code shapes} is null, the constructed S2ShapeIndexCoder will only
   * support encoding.
   *
   * @param shapes the list of shapes, used only by {@link #decode}, commonly the result of {@link
   *     VectorCoder#decode} for {@link VectorCoder#FAST_SHAPE}.
   */
  public S2ShapeIndexCoder(@Nullable List<S2Shape> shapes) {
    this.shapes = shapes;
  }

  /** Returns a coder that preloads all shapes. Largest heap usage but best performance. */
  public S2ShapeIndexCoder preloaded() {
    return new S2ShapeIndexCoder(shapes) {
      @Override
      public List<S2Shape> cacheShapes(List<S2Shape> shapes) {
        return ImmutableList.copyOf(shapes);
      }

      @Override
      public List<S2ClippedShape[]> cacheClippedShapes(List<S2ClippedShape[]> clippedShapes) {
        return ImmutableList.copyOf(clippedShapes);
      }
    };
  }

  /** Returns a coder that won't cache, for worst performance but the smallest heap usage. */
  public S2ShapeIndexCoder uncached() {
    return new S2ShapeIndexCoder(shapes) {
      @Override
      public List<S2Shape> cacheShapes(List<S2Shape> shapes) {
        return shapes;
      }

      @Override
      public List<S2ClippedShape[]> cacheClippedShapes(List<S2ClippedShape[]> clippedShapes) {
        return clippedShapes;
      }
    };
  }

  /** Encodes the given S2ShapeIndex into the given OutputStream. */
  @Override
  @JsIgnore // OutputStream is not available to J2CL.
  public void encode(S2ShapeIndex value, OutputStream output) throws IOException {
    try (Encoder encoder = new Encoder(value.options())) {
      encoder.init(value.getShapes(), output);
      for (S2Iterator<Cell> it = value.iterator(); !it.done(); it.next()) {
        encoder.addCell(it.entry());
      }
    }
  }

  /**
   * Encodes an index from the options, a list of all the shapes, and an iterator of the cells. This
   * allows creating much larger encoded indices than could fit in RAM in decoded form. To use this
   * class, call {@link #init} to set the shapes and output, call {@link #addCell} for each index
   * cell (hopefully as they're generated instead of collecting them first), and finally call {@link
   * #close} to write the bytes to output. This class does not flush or close the output.
   */
  public static class Encoder implements AutoCloseable {
    private final Options options;
    private final List<S2CellId> cellIds = new ArrayList<>();
    private final List<byte[]> encodedCells = new ArrayList<>();
    private final ByteArrayOutputStream baos = new ByteArrayOutputStream();
    private boolean oneShape;
    private Cell last;
    private OutputStream output;

    /** Create an encode for the given options. */
    public Encoder(Options options) {
      this.options = options;
    }

    /** Initializes this encoder for the given shapes, writing the version to 'output'. */
    @CanIgnoreReturnValue
    public Encoder init(List<S2Shape> shapes, OutputStream output) throws IOException {
      this.output = output;
      this.oneShape = shapes.size() == 1;
      // The version number is encoded in 2 bits, under the assumption that by the time we need 5
      // versions the first version can be permanently retired. This only saves 1 byte, but that's
      // significant for very small indexes.
      long maxEdges = options.getMaxEdgesPerCell();
      writeVarint64(output, maxEdges << 2 | EncodedS2ShapeIndex.CURRENT_ENCODING_VERSION);
      return this;
    }

    /**
     * Encodes 'cell' immediately and holds the bytes in memory, flushing into the output when
     * {@link #close} is called.
     *
     * <p>Cells must be added in strictly increasing S2Cellid order.
     */
    public void addCell(Cell cell) throws IOException {
      Preconditions.checkArgument(last == null || UnsignedLongs.compare(last.id(), cell.id()) < 0);
      last = cell;
      cellIds.add(new S2CellId(cell.id()));
      encodeCell(oneShape, cell, baos);
      encodedCells.add(baos.toByteArray());
      baos.reset();
    }

    @Override
    public void close() throws IOException {
      S2CellIdVectorCoder.INSTANCE.encode(cellIds, output);
      this.cellIds.clear();
      VectorCoder.BYTE_ARRAY.encode(encodedCells, output);
      this.encodedCells.clear();
      this.output = null;
      this.last = null;
    }
  }

  /**
   * Internal representation of an undecoded shape, which must be distinguished from a null shape.
   */
  private static final S2Shape UNDECODED_SHAPE = new S2EdgeVectorShape();

  /** Returns a cache of the given list of shapes. */
  public List<S2Shape> cacheShapes(List<S2Shape> shapes) {
    S2Shape[] cachedShapes = new S2Shape[Ints.checkedCast(shapes.size())];
    Arrays.fill(cachedShapes, UNDECODED_SHAPE);
    return new AbstractList<S2Shape>() {
      @Override
      public int size() {
        return cachedShapes.length;
      }

      @Override
      public synchronized S2Shape get(int i) {
        return cachedShapes[i] == UNDECODED_SHAPE
            ? cachedShapes[i] = shapes.get(i)
            : cachedShapes[i];
      }
    };
  }

  /** Returns a cache of cells given the IDs and clipped shapes in pairwise order. */
  public List<S2ClippedShape[]> cacheClippedShapes(List<S2ClippedShape[]> clippedShapes) {
    S2ClippedShape[][] cache = new S2ClippedShape[clippedShapes.size()][];
    return new AbstractList<>() {
      @Override
      public int size() {
        return cache.length;
      }

      @Override
      public synchronized S2ClippedShape[] get(int i) {
        if (cache[i] == null) {
          cache[i] = clippedShapes.get(i);
        }
        return cache[i];
      }
    };
  }

  /**
   * Returns an immutable {@link S2ShapeIndex} backed by {@code data} starting at {@code cursor}. A
   * list of the shapes in the index must have been previously provided to the constructor. The
   * {@code cursor.position} is updated to the position of the first byte in {@code data} following
   * the encoded value.
   *
   * <p>Values are decoded only when they are accessed. This allows for very fast initialization and
   * little additional memory use beyond the encoded data.
   */
  @Override
  public S2ShapeIndex decode(Bytes data, Cursor cursor) {
    Preconditions.checkNotNull(shapes);
    long maxEdgesVersion = data.readVarint64(cursor);
    int version = (int) maxEdgesVersion & 3;
    Preconditions.checkArgument(
        version == S2ShapeIndex.CURRENT_ENCODING_VERSION, "Unknown encoding.");

    S2ShapeIndex.Options options = new S2ShapeIndex.Options();
    options.setMaxEdgesPerCell(Ints.checkedCast(maxEdgesVersion >> 2));
    try {
      var encodedCellIds = S2CellIdVectorCoder.INSTANCE.decode(data, cursor);
      var encodedCells = new VectorCoder<>(new ClippedShapeCoder(shapes)).decode(data, cursor);
      return new EncodedS2ShapeIndex(
          options, cacheShapes(shapes), encodedCellIds, cacheClippedShapes(encodedCells));
    } catch (IOException e) {
      // TODO(user): Throw IOExceptions.
      throw new IllegalStateException("Underlying bad data / IO error ", e);
    }
  }

  @Override
  public boolean isLazy() {
    return true;
  }

  /**
   * Represents an encoded {@link S2ShapeIndex}.
   *
   * <p>This class is thread-safe.
   */
  private static final class EncodedS2ShapeIndex extends S2ShapeIndex {
    /** The encoded vector of cell IDs of this index. */
    private final S2CellIdVector encodedCellIds;

    /** The encoded cells of this index. */
    private final List<S2ClippedShape[]> encodedCells;

    /**
     * Initializes an {@link EncodedS2ShapeIndex} with the given {@code options}, backed by {@code
     * data} at {@code offset}.
     *
     * <p>Values are decoded only when they are accessed. This allows for very fast initialization
     * and little additional memory use beyond the encoded data.
     */
    @JsConstructor
    EncodedS2ShapeIndex(
        Options options,
        List<S2Shape> shapes,
        S2CellIdVector cellIds,
        List<S2ClippedShape[]> cells) {
      super(options);
      this.shapes = shapes;
      this.encodedCellIds = cellIds;
      this.encodedCells = cells;
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
    public ListIterator<S2ShapeIndex.Cell> iterator() {
      // Create a list of lazy cells that decodes the cell on access.
      List<S2ShapeIndex.Cell> cells = new AbstractList<>() {
        @Override
        public int size() {
          return encodedCells.size();
        }

        @Override
        public S2ShapeIndex.Cell get(int i) {
          return new LazyCell(encodedCellIds.get(i).id(), encodedCells.get(i));
        }
      };

      // Override the id() method to provide fast access to the ID without decoding the Cell.
      return new ListIterator<>(cells) {
        @Override
        public S2CellId id() {
          return encodedCellIds.get(pos);
        }
      };
    }

    @Override
    public boolean isFresh() {
      return true;
    }

    @Override
    public void applyUpdates() {
      // Encoded shape index never has updates to apply.
      return;
    }

    /** A lazy implementation of {@link S2ShapeIndex.Cell} which decodes members on demand. */
    public static final class LazyCell extends S2ShapeIndex.Cell {
      private final long cachedCellId;
      private final S2ClippedShape[] cachedClippedShapes;

      @JsConstructor
      LazyCell(long cellId, S2ClippedShape[] clippedShapes) {
        this.cachedCellId = cellId;
        this.cachedClippedShapes = clippedShapes;
      }

      @Override
      public long id() {
        return cachedCellId;
      }

      @Override
      public int numShapes() {
        return cachedClippedShapes.length;
      }

      @Override
      public S2ClippedShape clipped(int i) {
        return cachedClippedShapes[i];
      }
    }
  }

  /**
   * Encodes a cell into 'output'.
   *
   * @param oneShape true iff the index has a single shape
   * @param cell the cell to encode
   * @param output the output stream to receive the encoded bytes
   */
  private static void encodeCell(boolean oneShape, Cell cell, OutputStream output)
      throws IOException {
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

    if (oneShape) {
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
        writeVarint64(output, clipped.edge(0) << 6 | (n - 2) << 2 | containsCenter << 1);
      } else if (n == 1) {
        // The cell contains only one edge. For edge ids up to 15, we can encode the cell in a
        // single byte.
        //
        // Encoding: bits 0-1: 1
        //           bit 2: containsCenter
        //           bits 3+: edgeId
        writeVarint64(output, clipped.edge(0) << 3 | containsCenter << 2 | 1);
      } else {
        // General case (including n == 0, which is encoded compactly here).
        //
        // Encoding: bits 0-1: 3
        //           bit 2: containsCenter
        //           bits 3+: numEdges
        writeVarint64(output, n << 3 | containsCenter << 2 | 3);
        encodeEdges(clipped, output);
      }
    } else {
      // Note that there are exactly two possible values for the first encoded tag:
      // 1. numShapes() > 1:  a 3-bit tag with value 3.
      // 2. numShapes() == 0: a 3-bit tag with value 7.
      if (cell.numShapes() > 1) {
        // The cell contains more than one shape. The tag for this encoding must be distinguishable
        // from the cases encoded below. We can afford to use a 3-bit tag because numShapes() is
        // generally small.
        writeVarint64(output, (cell.numShapes() << 3) | 3);
      }
      // The shape ids are delta-encoded.
      int shapeIdBase = 0;
      for (int i = 0; i < cell.numShapes(); i++) {
        S2ClippedShape clipped = cell.clipped(i);
        int containsCenter = clipped.containsCenter() ? 1 : 0;
        int clippedShapeId = clipped.shapeId();
        assert clippedShapeId >= shapeIdBase;
        int shapeDelta = clippedShapeId - shapeIdBase;
        shapeIdBase = clippedShapeId + 1;

        // Like the code above except that we also need to encode shapeId(s). Because of this some
        // choices are slightly different.
        int n = clipped.numEdges();
        Preconditions.checkArgument(n < (1 << 29), "Too many edges.");

        if (n >= 1 && n <= 16 && clipped.edge(n - 1) - clipped.edge(0) == n - 1) {
          // The clipped shape has a contiguous range of up to 16 edges. This encoding uses a 1-bit
          // tag because it is by far the most common.
          //
          // Encoding: bit 0: 0
          //           bit 1: containsCenter
          //           bits 2+: edgeId
          // Next value: bits 0-3: (numEdges - 1)
          //             bits 4+: shapeDelta
          writeVarint64(output, (clipped.edge(0) << 2) | (containsCenter << 1));
          writeVarint64(output, (shapeDelta << 4) | (n - 1));
        } else if (n == 0) {
          // Special encoding for clipped shapes with no edges. Such shapes are common in polygon
          // interiors. This encoding uses a 3-bit tag in order to leave more bits available for the
          // other encodings.
          //
          // NOTE: When numShapes() > 1, this tag could be 2 bits (because the tag used to indicate
          // numShapes() > 1 can't appear). Alternatively, that tag can be considered reserved for
          // future use.
          //
          // Encoding: bits 0-2: 7
          //           bit 3: containsCenter
          //           bits 4+: shapeDelta
          writeVarint64(output, (shapeDelta << 4) | (containsCenter << 3) | 7);
        } else {
          // General case. This encoding uses a 2-bit tag, and the first value typically is encoded
          // into one byte.
          //
          // Encoding: bits 0-1: 1
          //           bit 2: containsCenter
          //           bits 3+: (numEdges - 1)
          // Next value: shapeDelta
          writeVarint64(output, ((n - 1) << 3) | (containsCenter << 2) | 1);
          writeVarint64(output, shapeDelta);
          encodeEdges(clipped, output);
        }
      }
    }
  }

  /** Encodes the edge IDs of the given {@link S2ClippedShape}. */
  private static void encodeEdges(S2ClippedShape clipped, OutputStream output) throws IOException {
    // Each entry is an (edgeId, count) pair representing a contiguous range of edges. The edge ids
    // are delta-encoded such that 0 represents the minimum valid next edge id.
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
        writeVarint64(output, delta);
      } else {
        // Count the edges in this contiguous range.
        int count = 1;
        for (; i + 1 < numEdges && clipped.edge(i + 1) == edgeId + count; i++) {
          count++;
        }
        if (count < 8) {
          // Count is encoded in low 3 bits of delta.
          writeVarint64(output, delta << 3 | (count - 1));
        } else {
          // Count and delta are encoded separately.
          writeVarint64(output, (count - 8) << 3 | 7);
          writeVarint64(output, delta);
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
   *
   * <p>Will throw IllegalArgumentException for inputs that are not valid encodings, but not
   * necessarily all invalid inputs will be detected.
   */
  private static S2ClippedShape[] decodeClippedShapes(
      List<S2Shape> shapes, Bytes data, Cursor cursor) {
    // This function inverts the encodings documented above.
    if (shapes.size() == 1) {
      S2ClippedShape[] clippedShapes = new S2ClippedShape[1];

      // Entire S2ShapeIndex contains only one shape.
      int header = data.readVarint32(cursor);
      if ((header & 1) == 0) {
        // The cell contains a contiguous range of edges.
        int numEdges = ((header >>> 2) & 15) + 2;
        boolean center = (header & 2) != 0;
        clippedShapes[0] = S2ClippedShape.create(null, 0, center, header >>> 6, numEdges);
      } else if ((header & 2) == 0) {
        // The cell contains a single edge.
        boolean center = (header & 4) != 0;
        clippedShapes[0] = S2ClippedShape.create(null, 0, center, header >>> 3, 1);
      } else {
        // The cell contains some other combination of edges.
        int numEdges = Ints.checkedCast(header >> 3);
        int[] edges = decodeEdges(numEdges, data, cursor);
        boolean center = (header & 4) != 0;
        clippedShapes[0] = ManyEdges.create(null, 0, center, edges);
      }
      return clippedShapes;
    }

    // S2ShapeIndex contains more than one shape.
    int header = data.readVarint32(cursor);
    int numClipped = 1;
    if ((header & 7) == 3) {
      // This cell contains more than one shape.
      numClipped = Ints.checkedCast(header >>> 3);
      header = data.readVarint32(cursor);
    }

    S2ClippedShape[] clippedShapes = new S2ClippedShape[numClipped];

    int shapeId = 0;
    for (int j = 0; j < numClipped; j++, shapeId++) {
      if (j > 0) {
        header = data.readVarint32(cursor);
      }
      if ((header & 1) == 0) {
        // The clipped shape contains a contiguous range of edges.
        int shapeIdCount = data.readVarint32(cursor);
        shapeId = addExact(shapeId, shapeIdCount >> 4);
        int numEdges = (shapeIdCount & 15) + 1;
        boolean center = (header & 2) != 0;
        clippedShapes[j] = S2ClippedShape.create(null, shapeId, center, header >>> 2, numEdges);
      } else if ((header & 7) == 7) {
        // The clipped shape has no edges.
        shapeId = addExact(shapeId, header >> 4);
        boolean center = (header & 8) != 0;
        clippedShapes[j] = S2ClippedShape.create(null, shapeId, center, 0, 0);
      } else {
        // The clipped shape contains some other combination of edges.
        assert (header & 3) == 1;
        int shapeDelta = data.readVarint32(cursor);
        shapeId = addExact(shapeId, shapeDelta);
        int numEdges = (header >>> 3) + 1;
        int[] edges = decodeEdges(numEdges, data, cursor);
        boolean center = (header & 4) != 0;
        clippedShapes[j] = ManyEdges.create(null, shapeId, center, edges);
      }
    }
    return clippedShapes;
  }

  /** A coder of {@code S2ClippedShape[]}s. */
  private static final class ClippedShapeCoder implements S2Coder<S2ClippedShape[]> {
    private final List<S2Shape> shapes;

    public ClippedShapeCoder(List<S2Shape> shapes) {
      this.shapes = shapes;
    }

    @Override
    @JsIgnore // OutputStream is not available to J2CL.
    public void encode(S2ClippedShape[] values, OutputStream output) {
      throw new UnsupportedOperationException();
    }

    @Override
    public S2ClippedShape[] decode(Bytes data, Cursor cursor) {
      return decodeClippedShapes(shapes, data, cursor);
    }

    @Override
    public boolean isLazy() {
      return true;
    }
  }
}
