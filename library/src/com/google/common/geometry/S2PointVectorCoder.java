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

import static com.google.common.geometry.S2Projections.PROJ;
import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.PrimitiveArrays.Longs;
import com.google.common.geometry.S2Projections.FaceSiTi;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.ImmutableLongArray;
import com.google.common.primitives.Ints;
import com.google.common.primitives.UnsignedInts;
import com.google.common.primitives.UnsignedLongs;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.AbstractList;
import java.util.ArrayList;
import java.util.List;

/**
 * An encoder/decoder of Lists of {@link S2Point}s.
 *
 * <p>Values from the List of {@link S2Point} returned by {@link #decode(Bytes, Cursor)} are decoded
 * only when they are accessed. This allows for very fast initialization and no additional memory
 * use beyond the encoded data. The encoded data is not owned by the {@code List}; typically it
 * points into a large contiguous buffer that contains other encoded data as well. This buffer may
 * not be modified after decode() is called: {@link S2Coder#isLazy()} is true.
 */
public class S2PointVectorCoder implements S2Coder<List<S2Point>> {
  /**
   * An instance of a {@code S2PointVectorCoder} which encodes/decodes {@link S2Point}s in the FAST
   * encoding format. The FAST format is optimized for fast encoding/decoding. Decoding is on-
   * demand, {@link S2Coder#isLazy()} is true.
   */
  public static final S2PointVectorCoder FAST = new S2PointVectorCoder(Format.FAST);

  /**
   * An instance of a {@code S2PointVectorCoder} which encodes/decodes {@link S2Point}s in the
   * COMPACT encoding format. The COMPACT format is optimized for disk usage and memory footprint.
   * Decoding is on-demand, {@link S2Coder#isLazy()} is true.
   */
  public static final S2PointVectorCoder COMPACT = new S2PointVectorCoder(Format.COMPACT);

  /**
   * Controls whether to optimize for speed or size when encoding points. (Note that encoding is
   * always lossless, and that currently compact encodings are only possible when points have been
   * snapped to S2CellId centers).
   */
  private enum Format {
    FAST,
    COMPACT
  }

  private final Format type;

  /** The value of the FAST encoding format. */
  private static final int FORMAT_FAST = 0;
  /** The value of the COMPACT encoding format. */
  private static final int FORMAT_COMPACT = 1;

  // To save space (especially for vectors of length 0, 1, and 2), the encoding format is encoded in
  // the low-order 3 bits of the vector size. Up to 7 encoding formats are supported (only 2 are
  // currently defined). Additional formats could be supported by using "7" as an overflow indicator
  // and encoding the actual format separately, but it seems unlikely we will ever need to do that.
  private static final int ENCODING_FORMAT_BITS = 3;
  private static final byte ENCODING_FORMAT_MASK = (1 << ENCODING_FORMAT_BITS) - 1;

  /** The size of an encoded {@link S2Point} in bytes (3 doubles * 8 bytes per double). */
  private static final int SIZEOF_S2POINT = 3 * 8;

  /** The left shift factor for {@link #BLOCK_SIZE}. */
  private static final int BLOCK_SHIFT = 4;

  /**
   * {@link S2CellId}s are represented in a special 64-bit format and are encoded in fixed-size
   * blocks. {@code BLOCK_SIZE} represents the number of values per block. Block sizes of 4, 8, 16,
   * and 32 were tested, and {@code BLOCK_SIZE} == 16 seems to offer the best compression. (Note
   * that {@code BLOCK_SIZE == 32} requires some code modifications which have since been removed).
   */
  private static final int BLOCK_SIZE = 1 << BLOCK_SHIFT;

  /** The exception value in the COMPACT encoding format. */
  private static final long EXCEPTION = S2CellId.sentinel().id();

  private S2PointVectorCoder(Format type) {
    this.type = type;
  }

  @Override
  public void encode(List<S2Point> values, OutputStream output) throws IOException {
    switch (type) {
      case FAST:
        encodeFast(values, output);
        break;
      case COMPACT:
        encodeCompact(values, output);
        break;
    }
  }

  @Override
  public List<S2Point> decode(Bytes data, Cursor cursor) throws IOException {
    // Peek at the format but don't advance the decoder.
    int format = data.get(cursor.position) & ENCODING_FORMAT_MASK;

    switch (format) {
      case FORMAT_FAST:
        return decodeFast(data, cursor);
      case FORMAT_COMPACT:
        return decodeCompact(data, cursor);
      default:
        throw new IllegalArgumentException("Unexpected format: " + format);
    }
  }

  @Override
  public boolean isLazy() {
    // Both FORMAT_FAST and FORMAT_COMPACT are lazy.
    return true;
  }

  private static void encodeFast(List<S2Point> values, OutputStream output) throws IOException {
    // The encoding format is as follows:
    //
    //   varint64, bits 0-2:  encoding format (UNCOMPRESSED)
    //             bits 3-63: vector size
    //   array of values.size() S2Points in little-endian order

    long sizeFormat = (((long) values.size()) << ENCODING_FORMAT_BITS) | FORMAT_FAST;
    EncodedInts.writeVarint64(output, sizeFormat);
    for (S2Point point : values) {
      point.encode(output);
    }
  }

  /** Decodes S2Points {@link S2Coder#isLazy on demand}. */
  private static List<S2Point> decodeFast(Bytes data, Cursor cursor) {
    long tmpSize = data.readVarint64(cursor);
    tmpSize >>= ENCODING_FORMAT_BITS;
    int size = Ints.checkedCast(tmpSize);
    long offset = cursor.position;
    cursor.position += (long) size * (long) SIZEOF_S2POINT;

    return new AbstractList<S2Point>() {
      @Override
      public S2Point get(int index) {
        long position = offset + (long) index * (long) SIZEOF_S2POINT;
        return new S2Point(
            data.readLittleEndianDouble(position),
            data.readLittleEndianDouble(position + Doubles.BYTES),
            data.readLittleEndianDouble(position + Doubles.BYTES * 2));
      }

      @Override
      public int size() {
        return size;
      }
    };
  }

  /**
   * Encodes a vector of {@link S2Point}s, optimizing for space.
   *
   * <p><strong>OVERVIEW</strong>
   *
   * <p>We attempt to represent each S2Point as the center of an S2CellId. All S2CellIds must be at
   * the same level. Any points that cannot be encoded exactly as S2CellId centers are stored as
   * exceptions using 24 bytes each. If there are so many exceptions that the COMPACT encoding does
   * not save significant space, we give up and use the uncompressed encoding.
   *
   * <p>The first step is to choose the best S2CellId level. This requires converting each point to
   * (face, si, ti) coordinates and checking whether the point can be represented exactly as an
   * S2CellId center at some level. We then build a histogram of S2CellId levels (just like the
   * similar code in S2Polygon.encode) and choose the best level (or give up, if there are not
   * enough S2CellId-encodable points).
   *
   * <p>The simplest approach would then be to take all the S2CellIds and right-shift them to remove
   * all the constant bits at the chosen level. This would give the best spatial locality and hence
   * the smallest deltas. However instead we give up some spatial locality and use the similar but
   * faster transformation described below.
   *
   * <p>Each encodable point is first converted to the (sj, tj) representation defined below:
   *
   * <pre>{@code
   * sj = (((face & 3) << 30) | (si >> 1)) >> (30 - level);
   * tj = (((face & 4) << 29) | ti) >> (31 - level);
   * }</pre>
   *
   * These two values encode the (face, si, ti) tuple using (2 * level + 3) bits. To see this,
   * recall that "si" and "ti" are 31-bit values that all share a common suffix consisting of a "1"
   * bit followed by (30 - level) "0" bits. The code above right-shifts these values to remove the
   * constant bits and then prepends the bits for "face", yielding a total of (level + 2) bits for
   * "sj" and (level + 1) bits for "tj".
   *
   * <p>We then combine (sj, tj) into one 64-bit value by interleaving bit pairs:
   *
   * <pre>{@code
   * v = interleaveBitPairs(sj, tj);
   * }</pre>
   *
   * (We could also interleave individual bits, but it is faster this way). The result is similar to
   * right-shifting an S2CellId by (61 - 2 * level), except that it is faster to decode and the
   * spatial locality is not quite as good.
   *
   * <p>The 64-bit values are divided into blocks of size 8, and then each value is encoded as the
   * sum of a base value, a per-block offset, and a per-value delta within that block:
   *
   * <pre>{@code
   * v[i,j] = base + offset[i] + delta[i, j]
   * }</pre>
   *
   * <p>where "i" represents a block and "j" represents an entry in that block.
   *
   * <p>The deltas for each block are encoded using a fixed number of 4-bit nibbles (1-16 nibbles
   * per delta). This allows any delta to be accessed in constant time.
   *
   * <p>The "offset" for each block is a 64-bit value encoded in 0-8 bytes. The offset is left-
   * shifted such that it overlaps the deltas by a configurable number of bits (either 0 or 4),
   * called the "overlap". The overlap and offset length (0-8 bytes) are specified per block. The
   * reason for the overlap is that it allows fewer delta bits to be used in some cases. For example
   * if base == 0 and the range within a block is 0xf0 to 0x110, then rather than using 12-bits
   * deltas with an offset of 0, the overlap lets us use 8-bits deltas with an offset of 0xf0
   * (saving 7 bytes per block).
   *
   * <p>The global minimum value "base" is encoded using 0-7 bytes starting with the
   * most-significant non-zero bit possible for the chosen level. For example, if (level == 7) then
   * the encoded values have at most 17 bits, so if "base" is encoded in 1 byte then it is shifted
   * to occupy bits 9-16.
   *
   * <p>Example: at level == 15, there are at most 33 non-zero value bits. The following shows the
   * bit positions covered by "base", "offset", and "delta" assuming that "base" and "offset" are
   * encoded in 2 bytes each, deltas are encoded in 2 nibbles (1 byte) each, and "overlap" is 4
   * bits:
   *
   * <pre>{@code
   * Base:             1111111100000000-----------------
   * Offset:           -------------1111111100000000----
   * Delta:            -------------------------00000000
   * Overlap:                                   ^^^^
   * }</pre>
   *
   * <p>The numbers (0 or 1) in this diagram denote the byte number of the encoded value. Notice
   * that "base" is shifted so that it starts at the leftmost possible bit, "delta" always starts at
   * the rightmost possible bit (bit 0), and "offset" is shifted so that it overlaps "delta" by the
   * chosen "overlap" (either 0 or 4 bits). Also note that all of these values are summed, and
   * therefore each value can affect higher-order bits due to carries.
   *
   * <p>NOTE: Encoding deltas in 4-bit rather than 8-bit length increments reduces encoded sizes by
   * about 7%. Allowing a 4-bit overlap between the offset and deltas reduces encoded sizes by about
   * 1%. Both optimizations make the code more complex but don't affect running times significantly.
   *
   * <p><strong>ENCODING DETAILS</strong></br>
   *
   * <p>Now we can move on to the actual encodings. First, there is a 2 byte header encoded as
   * follows:
   *
   * <pre>{@code
   * Byte 0, bits 0-2: encodingFormat (COMPACT)
   * Byte 0, bit  3:   haveExceptions
   * Byte 0, bits 4-7: (lastBlockSize - 1)
   * Byte 1, bits 0-2: baseBytes
   * Byte 1, bits 3-7: level (0-30)
   * }</pre>
   *
   * <p>This is followed by an EncodedStringVector containing the encoded blocks. Each block
   * contains BLOCK_SIZE (8) values. The total size of the EncodedS2PointVector is not stored
   * explicitly, but instead is calculated as
   *
   * <pre>
   * num_values == BLOCK_SIZE * (numBlocks - 1) + lastBlockSize
   * </pre>
   *
   * <p>(An empty vector has numBlocks == 0 and lastBlockSize == BLOCK_SIZE.)
   *
   * <p>Each block starts with a 1 byte header containing the following:
   *
   * <pre>{@code
   * Byte 0, bits 0-2: (offsetBytes - overlapNibbles)
   * Byte 0, bit  3:   overlapNibbles
   * Byte 0, bits 4-7: (deltaNibbles - 1)
   * }</pre>
   *
   * <p>"overlapNibbles" is either 0 or 1 (indicating an overlap of 0 or 4 bits), while
   * "offsetBytes" is in the range 0-8 (indicating the number of bytes used to encode the offset for
   * this block). Note that some combinations cannot be encoded: in particular, offsetBytes == 0 can
   * only be encoded with an overlap of 0 bits, and offsetBytes == 8 can only be encoded with an
   * overlap of 4 bits. This allows us to encode offset lengths of 0-8 rather than just 0-7 without
   * using an extra bit. (Note that the combinations that can't be encoded are not useful anyway).
   *
   * <p>The header is followed by "offsetBytes" bytes for the offset, and then (4 * deltaNibbles)
   * bytes for the deltas.
   *
   * <p>If there are any points that could not be represented as S2CellIds, then "haveExceptions" in
   * the header is true. In that case the delta values within each block are encoded as (delta + 8),
   * and values 0-7 are used to represent exceptions. If a block has exceptions, they are encoded
   * immediately following the array of deltas, and are referenced by encoding the corresponding
   * exception index (0-7) as the delta.
   *
   * <p>TODO(user): A vector containing a single leaf cell is currently encoded as 13 bytes
   * (2 byte header, 7 byte base, 1 byte block count, 1 byte block length, 1 byte block header, 1
   * byte delta). However if this case occurs often, a better solution would be implement a separate
   * format that encodes the leading k bytes of an S2CellId. It would have a one-byte header
   * consisting of the encoding format (3 bits) and the number of bytes encoded (3 bits), followed
   * by the S2CellId bytes. The extra 2 header bits could be used to store single points using other
   * encodings, e.g. E7.
   *
   * <p>If we wind up using 8-value blocks, we could also use the extra bit in the first byte of the
   * header to indicate that there is only one value, and then skip the 2nd byte of header and the
   * EncodedStringVector. But this would be messy because it also requires special cases while
   * decoding. Essentially this would be a sub-format within the COMPACT format.
   */
  private static void encodeCompact(List<S2Point> values, OutputStream output) throws IOException {
    // 1. Compute (level, face, si, ti) for each point, build a histogram of levels, and determine
    // the optimal level to use for encoding (if any).
    List<CellPoint> cellPoints = Lists.newArrayListWithCapacity(values.size());
    int level = chooseBestLevel(values, cellPoints);
    if (level < 0) {
      encodeFast(values, output);
      return;
    }

    // 2. Convert the points into encodable 64-bit values. We don't use the S2CellId itself
    // because it requires a somewhat more complicated bit interleaving operation.
    //
    // TODO(user): Benchmark using shifted S2CellIds instead.
    ImmutableLongArray cellPointValues = convertCellsToValues(cellPoints, level);
    boolean haveExceptions = cellPointValues.contains(EXCEPTION);

    // 3. Choose the global encoding parameter "base" (consisting of the bit prefix shared by all
    // values to be encoded).
    Base base = chooseBase(cellPointValues, level, haveExceptions);

    // Now encode the output, starting with the 2-byte header (see above).
    int numBlocks = (cellPointValues.length() + BLOCK_SIZE - 1) >> BLOCK_SHIFT;
    int baseBytes = base.baseBits >> 3;
    int lastBlockCount = cellPointValues.length() - BLOCK_SIZE * (numBlocks - 1);
    Preconditions.checkArgument(lastBlockCount >= 0);
    Preconditions.checkArgument(lastBlockCount <= BLOCK_SIZE);
    Preconditions.checkArgument(baseBytes <= 7);
    Preconditions.checkArgument(level <= 30);
    output.write(FORMAT_COMPACT | ((haveExceptions ? 1 : 0) << 3) | ((lastBlockCount - 1) << 4));
    output.write(baseBytes | (level << 3));

    // Next we encode 0-7 bytes of "base".
    int baseShift = baseShift(level, base.baseBits);
    EncodedInts.encodeUintWithLength(output, base.base >> baseShift, baseBytes);

    // Now we encode the contents of each block.
    List<byte[]> blocks = new ArrayList<>();
    List<S2Point> exceptions = new ArrayList<>();
    MutableBlockCode code = new MutableBlockCode();
    for (int i = 0; i < cellPointValues.length(); i += BLOCK_SIZE) {
      int blockSize = min(BLOCK_SIZE, cellPointValues.length() - i);
      getBlockCode(code, cellPointValues.subArray(i, i + blockSize), base.base, haveExceptions);

      // Encode the one-byte block header (see above).
      ByteArrayOutput block = new ByteArrayOutput();
      int offsetBytes = code.offsetBits >> 3;
      int deltaNibbles = code.deltaBits >> 2;
      int overlapNibbles = code.overlapBits >> 2;
      Preconditions.checkArgument((offsetBytes - overlapNibbles) <= 7);
      Preconditions.checkArgument(overlapNibbles <= 1);
      Preconditions.checkArgument(deltaNibbles <= 16);
      block.write((offsetBytes - overlapNibbles) | (overlapNibbles << 3) | (deltaNibbles - 1) << 4);

      // Determine the offset for this block, and whether there are exceptions.
      long offset = -1L;
      int numExceptions = 0;
      for (int j = 0; j < blockSize; j++) {
        if (cellPointValues.get(i + j) == EXCEPTION) {
          numExceptions += 1;
        } else {
          Preconditions.checkArgument(cellPointValues.get(i + j) >= base.base);
          offset = UnsignedLongs.min(offset, cellPointValues.get(i + j) - base.base);
        }
      }
      if (numExceptions == blockSize) {
        offset = 0;
      }

      // Encode the offset.
      int offsetShift = code.deltaBits - code.overlapBits;
      offset &= ~bitMask(offsetShift);
      Preconditions.checkArgument((offset == 0) == (offsetBytes == 0));
      if (offset > 0) {
        EncodedInts.encodeUintWithLength(block, offset >>> offsetShift, offsetBytes);
      }

      // Encode the deltas, and also gather any exceptions present.
      int deltaBytes = (deltaNibbles + 1) >> 1;
      exceptions.clear();
      for (int j = 0; j < blockSize; j++) {
        long delta;
        if (cellPointValues.get(i + j) == EXCEPTION) {
          delta = exceptions.size();
          exceptions.add(values.get(i + j));
        } else {
          Preconditions.checkArgument(
              UnsignedLongs.compare(cellPointValues.get(i + j), offset + base.base) >= 0);
          delta = cellPointValues.get(i + j) - (offset + base.base);
          if (haveExceptions) {
            Preconditions.checkArgument(UnsignedLongs.compare(delta, -1L - BLOCK_SIZE) <= 0);
            delta += BLOCK_SIZE;
          }
        }
        Preconditions.checkArgument(UnsignedLongs.compare(delta, bitMask(code.deltaBits)) <= 0);
        if (((deltaNibbles & 1) != 0) && ((j & 1) != 0)) {
          // Combine this delta with the high-order 4 bits of the previous delta.
          int lastByte = block.removeLast();
          delta = (delta << 4) | (lastByte & 0xf);
        }
        EncodedInts.encodeUintWithLength(block, delta, deltaBytes);
      }
      // Append any exceptions to the end of the block.
      if (numExceptions > 0) {
        for (S2Point p : exceptions) {
          p.encode(block);
        }
      }
      blocks.add(block.toByteArray());
    }
    VectorCoder.BYTE_ARRAY.encode(blocks, output);
  }

  private static List<S2Point> decodeCompact(Bytes data, Cursor cursor) throws IOException {
    // First we decode the two-byte header:
    //  Byte 0, bits 0-2: encodingFormat (COMPACT)
    //  Byte 0, bit  3:   haveExceptions
    //  Byte 0, bits 4-7: (lastBlockSize - 1)
    //  Byte 1, bits 0-2: baseBytes
    //  Byte 1, bits 3-7: level (0-30)
    int header1 = data.get(cursor.position++) & 0xFF;
    int header2 = data.get(cursor.position++) & 0xFF;
    Preconditions.checkArgument((header1 & 7) == FORMAT_COMPACT);

    boolean haveExceptions = (header1 & 8) != 0;
    int lastBlockCount = (header1 >> 4) + 1;
    int baseBytes = header2 & 7;
    int level = header2 >> 3;

    // Decode the base value (if any).
    long tmpBase = data.readUintWithLength(cursor, baseBytes);
    long base = tmpBase << baseShift(level, baseBytes << 3);

    // Initialize the vector of encoded blocks.
    Longs blockOffsets = UintVectorCoder.UINT64.decode(data, cursor);
    long offset = cursor.position;

    int size = BLOCK_SIZE * (blockOffsets.length() - 1) + lastBlockCount;
    cursor.position += blockOffsets.length() > 0 ? blockOffsets.get(blockOffsets.length() - 1) : 0;

    return new AbstractList<S2Point>() {
      @Override
      public S2Point get(int index) {
        // First we decode the block header.
        int iShifted = index >> BLOCK_SHIFT;
        long position = offset + ((iShifted == 0) ? 0 : blockOffsets.get(iShifted - 1));
        int header = data.get(position++) & 0xFF;
        int overlapNibbles = (header >> 3) & 1;
        int offsetBytes = (header & 7) + overlapNibbles;
        int deltaNibbles = (header >> 4) + 1;

        // Decode the offset for this block.
        int offsetShift = (deltaNibbles - overlapNibbles) << 2;
        long offset;
        long delta;
        try {
          offset = data.readUintWithLength(data.cursor(position), offsetBytes) << offsetShift;
          position += offsetBytes;
          long exceptionPosition = position;

          // Decode the delta for the requested value.
          int deltaNibbleOffset = (index & (BLOCK_SIZE - 1)) * deltaNibbles;
          int deltaBytes = (deltaNibbles + 1) >> 1;
          position += deltaNibbleOffset >> 1;
          delta = data.readUintWithLength(data.cursor(position), deltaBytes);
          delta >>>= (deltaNibbleOffset & 1) << 2;
          delta &= bitMask(deltaNibbles << 2);

          // Test whether this point is encoded as an exception.
          if (haveExceptions) {
            if (delta < BLOCK_SIZE) {
              int blockSize = min(BLOCK_SIZE, size - (index & -BLOCK_SIZE));
              exceptionPosition += (blockSize * deltaNibbles + 1) >> 1;
              exceptionPosition += delta * SIZEOF_S2POINT;
              return S2Point.decode(data.toInputStream(exceptionPosition));
            }
            delta -= BLOCK_SIZE;
          }
        } catch (IOException e) {
          // This should never happen because Bytes.get() does not throw an IOException.
          throw new RuntimeException(e);
        }

        // Otherwise convert the 64-bit value back to an S2Point.
        long value = base + offset + delta;
        int shift = S2CellId.MAX_LEVEL - level;

        // The S2CellId version of the following code is:
        //   return S2CellId(((value << 1) | 1) << (2 * shift)).ToPoint();
        int sj = EncodedInts.deinterleaveBitPairs1(value);
        int tj = EncodedInts.deinterleaveBitPairs2(value);
        int si = (((sj << 1) | 1) << shift) & 0x7fffffff;
        int ti = (((tj << 1) | 1) << shift) & 0x7fffffff;
        int face = ((sj << shift) >>> 30) | (((tj << (shift + 1)) >>> 29) & 4);
        return S2Projections.faceUvToXyz(
                face,
                PROJ.stToUV(S2Projections.siTiToSt(si)),
                PROJ.stToUV(S2Projections.siTiToSt(ti)))
            .normalize();
      }

      @Override
      public int size() {
        return size;
      }
    };
  }

  /**
   * Represents the encoding parameters to be used for a given block (consisting of {@link
   * #BLOCK_SIZE} encodable 64-bit values).
   */
  private static final class MutableBlockCode {
    /** Delta length in bits (multiple of 4). */
    int deltaBits;
    /** Offset length in bits (multiple of 8). */
    int offsetBits;
    /** {Delta, Offset} overlap in bits (0 or 4). */
    int overlapBits;

    MutableBlockCode() {}

    public void set(int deltaBits, int offsetBits, int overlapBits) {
      this.deltaBits = deltaBits;
      this.offsetBits = offsetBits;
      this.overlapBits = overlapBits;
    }
  }

  /** Return type of {@link #chooseBase}. */
  private static final class Base {
    long base;
    int baseBits;

    Base(long base, int baseBits) {
      this.base = base;
      this.baseBits = baseBits;
    }
  }

  /**
   * Represents a point that can be encoded as an {@link S2CellId} center.
   *
   * <p>(If such an encoding is not possible then level < 0).
   */
  private static final class CellPoint {
    short level;
    short face;
    int si;
    int ti;

    CellPoint(int level, FaceSiTi faceSiTi) {
      this.level = (short) level;
      // TODO(user): UnsignedBytes is marked GwtIncompatible. Use that here when we can.
      Preconditions.checkArgument(faceSiTi.face >> Byte.SIZE == 0);
      this.face = (byte) faceSiTi.face;
      this.si = UnsignedInts.checkedCast(faceSiTi.si);
      this.ti = UnsignedInts.checkedCast(faceSiTi.ti);
    }
  }

  /**
   * A thin wrapper over {@link ByteArrayOutputStream} which allows the last written byte to be
   * removed.
   */
  private static class ByteArrayOutput extends ByteArrayOutputStream {
    /** Removes and returns the last written byte. */
    int removeLast() {
      Preconditions.checkState(count > 0);
      return buf[--count];
    }
  }

  /** Returns a bit mask with {@code n} low-order 1 bits, for {@code 0 <= n <= 64}. */
  private static long bitMask(int n) {
    return (n == 0) ? 0 : (-1L >>> (64 - n));
  }

  /** Returns the maximum number of bits per value at the given {@link S2CellId} level. */
  private static int maxBitsForLevel(int level) {
    return 2 * level + 3;
  }

  /**
   * Returns the number of bits that {@code base} should be right-shifted in order to encode only
   * its leading {@code baseBits} bits, assuming that all points are encoded at the given {@link
   * S2CellId} level.
   */
  private static int baseShift(int level, int baseBits) {
    return max(0, maxBitsForLevel(level) - baseBits);
  }

  /**
   * Returns the {@code S2CellId} level for which the greatest number of the given points can be
   * represented as the center of an {@code S2CellId}, or -1 if there is no S2CellId that would
   * result in significant space savings.
   *
   * <p>Adds the {@code S2CellId} representation of each point (if any) to {@code cellPoints}.
   */
  private static int chooseBestLevel(List<S2Point> points, List<CellPoint> cellPoints) {
    // Count the number of points at each level.
    int[] levelCounts = new int[S2CellId.MAX_LEVEL + 1];
    for (S2Point point : points) {
      FaceSiTi faceSiTi = PROJ.xyzToFaceSiTi(point);
      int level = PROJ.levelIfCenter(faceSiTi, point);
      cellPoints.add(new CellPoint(level, faceSiTi));
      if (level >= 0) {
        levelCounts[level]++;
      }
    }
    // Choose the level for which the most points can be encoded.
    int bestLevel = 0;
    for (int level = 1; level <= S2CellId.MAX_LEVEL; level++) {
      if (levelCounts[level] > levelCounts[bestLevel]) {
        bestLevel = level;
      }
    }
    // The uncompressed encoding is smaller *and* faster when very few of the points are encodable
    // as S2CellIds. The COMPACT encoding uses about 1 extra byte per point in this case,
    // consisting of up to a 3 byte EncodedArray offset for each block, a 1 byte block header, and
    // 4 bits per delta (encoding an exception number from 0-7), for a total of 8 bytes per block.
    // This represents a space overhead of about 4%, so we require that at least 5% of the input
    // points should be encodable as S2CellIds in order for the COMPACT format to be worthwhile.
    double minEncodableFraction = 0.05;
    if (levelCounts[bestLevel] <= minEncodableFraction * points.size()) {
      return -1;
    }
    return bestLevel;
  }

  /**
   * Given a vector of points in {@link CellPoint} format and an {@link S2CellId} level that has
   * been chosen for encoding, returns a vector of 64-bit values that should be encoded in order to
   * represent these points. Points that cannot be represented losslessly as the center of an {@code
   * S2CellId} at the chosen level are indicated by the value {@link #EXCEPTION}.
   */
  private static ImmutableLongArray convertCellsToValues(List<CellPoint> cellPoints, int level) {
    ImmutableLongArray.Builder builder = ImmutableLongArray.builder(cellPoints.size());
    int shift = S2CellId.MAX_LEVEL - level;
    for (CellPoint cp : cellPoints) {
      if (cp.level != level) {
        builder.add(EXCEPTION);
      } else {
        // Note that bit 31 of tj is always zero, and that bits are interleaved in
        // such a way that bit 63 of the result is always zero.
        //
        // The S2CellId version of the following code is:
        // long v = S2CellId.fromFaceIJ(cp.face, cp.si >> 1, cp.ti >> 1)
        //     .parent(level).id() >> (2 * shift + 1);
        int sj = (((cp.face & 3) << 30) | (cp.si >>> 1)) >>> shift;
        int tj = (((cp.face & 4) << 29) | cp.ti) >>> (shift + 1);
        long v = EncodedInts.interleaveBitPairs(sj, tj);
        Preconditions.checkArgument(UnsignedLongs.compare(v, bitMask(maxBitsForLevel(level))) <= 0);
        builder.add(v);
      }
    }
    return builder.build();
  }

  /**
   * Returns the global minimum value {@link Base#base} and the number of bits that should be used
   * to encode it ({@link Base#baseBits}).
   */
  private static Base chooseBase(ImmutableLongArray values, int level, boolean haveExceptions) {
    // Find the minimum and maximum non-exception values to be represented.
    long vMin = EXCEPTION;
    long vMax = 0;
    for (int i = 0; i < values.length(); i++) {
      long v = values.get(i);
      if (v != EXCEPTION) {
        vMin = UnsignedLongs.min(vMin, v);
        vMax = UnsignedLongs.max(vMax, v);
      }
    }
    if (vMin == EXCEPTION) {
      return new Base(0, 0);
    }

    // Generally "base" is chosen as the bit prefix shared by vMin and vMax. However, there are a
    // few adjustments we need to make.
    //
    // 1. Encodings are usually smaller if the bits represented by "base" and "delta" do not
    // overlap. Usually the shared prefix rule does this automatically, but if vMin == vMax or
    // there are special circumstances that increase deltaBits (such as values.size() == 1) then
    // we need to make an adjustment.
    //
    // 2. The format only allows us to represent up to 7 bytes (56 bits) of "base", so we need to
    // ensure that "base" conforms to this requirement.
    int minDeltaBits = (haveExceptions || values.length() == 1) ? 8 : 4;
    int excludedBits =
        Ints.max(
            63 - Long.numberOfLeadingZeros(vMin ^ vMax) + 1, minDeltaBits, baseShift(level, 56));
    long base = vMin & ~bitMask(excludedBits);

    int baseBits = 0;
    // Determine how many bytes are needed to represent this prefix.
    if (base != 0) {
      int lowBit = Long.numberOfTrailingZeros(base);
      baseBits = (maxBitsForLevel(level) - lowBit + 7) & ~7;
    }

    // Since baseBits has been rounded up to a multiple of 8, we may now be able to represent
    // additional bits of vMin. In general this reduces the final encoded size.
    //
    // NOTE(ericv): A different strategy for choosing "base" is to encode all blocks under the
    // assumption that "base" equals vMin exactly, and then set base equal to the minimum-length
    // prefix of "vMin" that allows these encodings to be used. This strategy reduces the encoded
    // sizes by about 0.2% relative to the strategy here, but is more complicated.
    return new Base(vMin & ~bitMask(baseShift(level, baseBits)), baseBits);
  }

  /**
   * Returns true if the range of values {@code [dMin, dMax]} can be encoded using the specified
   * parameters ({@code deltaBits}, {@code overlapBits}, and {@code haveExceptions}).
   */
  private static boolean canEncode(
      long dMin, long dMax, int deltaBits, int overlapBits, boolean haveExceptions) {
    // "offset" can't represent the lowest (deltaBits - overlapBits) of dMin.
    dMin &= ~bitMask(deltaBits - overlapBits);

    // The maximum delta is reduced by BLOCK_SIZE if any exceptions exist, since deltas
    // 0..BLOCK_SIZE-1 are used to indicate exceptions.
    long maxDelta = bitMask(deltaBits);
    if (haveExceptions) {
      if (UnsignedLongs.compare(maxDelta, BLOCK_SIZE) < 0) {
        return false;
      }
      maxDelta -= BLOCK_SIZE;
    }
    // The first test below is necessary to avoid 64-bit overflow.
    return (UnsignedLongs.compare(dMin, ~maxDelta) > 0)
        || (UnsignedLongs.compare(dMin + maxDelta, dMax) >= 0);
  }

  /**
   * Given a vector of 64-bit values to be encoded and an {@link S2CellId} level, returns the
   * optimal encoding parameters that should be used to encode each block.
   */
  private static void getBlockCode(
      MutableBlockCode code, ImmutableLongArray values, long base, boolean haveExceptions) {
    // "bMin" and "bMax"n are the minimum and maximum values within this block.
    long bMin = EXCEPTION;
    long bMax = 0;
    for (int i = 0; i < values.length(); i++) {
      long v = values.get(i);
      if (v != EXCEPTION) {
        bMin = UnsignedLongs.min(bMin, v);
        bMax = UnsignedLongs.max(bMax, v);
      }
    }
    if (bMin == EXCEPTION) {
      // All values in this block are exceptions.
      code.set(4, 0, 0);
      return;
    }

    // Adjust the min/max values so that they are relative to "base".
    bMin -= base;
    bMax -= base;

    // Determine the minimum possible delta length and overlap that can be used to encode this
    // block. The block will usually be encodable using the number of bits in (bMax - bMin)
    // rounded up to a multiple of 4. If this is not possible, the preferred solution is to shift
    // "offset" so that the delta and offset values overlap by 4 bits (since this only costs an
    // average of 4 extra bits per block). Otherwise we increase the delta size by 4 bits. Certain
    // cases require that both of these techniques are used.
    //
    // Example 1: bMin = 0x72, bMax = 0x7e. The range is 0x0c. This can be encoded using
    // deltaBits = 4 and overlapBits = 0, which allows us to represent an offset of 0x70 and a
    // maximum delta of 0x0f, so that we can encode values up to 0x7f.
    //
    // Example 2: bMin = 0x78, bMax = 0x84. The range is 0x0c, but in this case it is not
    // sufficient to use deltaBits = 4 and overlapBits = 0 because we can again only represent an
    // offset of 0x70, so the maximum delta of 0x0f only lets us encode values up to 0x7f. However
    // if we increase the overlap to 4 bits then we can represent an offset of 0x78, which lets us
    // encode values up to 0x78 + 0x0f = 0x87.
    //
    // Example 3: bMin = 0x08, bMax = 0x104. The range is 0xfc, so we should be able to use 8-bit
    // deltas. But even with a 4-bit overlap, we can still only encode offset = 0 and a maximum
    // value of 0xff. (We don't allow bigger overlaps because statistically they are not
    // worthwhile). Instead we increase the delta size to 12 bits, which handles this case easily.
    //
    // Example 4: bMin = 0xf08, bMax = 0x1004. The range is 0xfc, so we should be able to use
    // 8-bit deltas. With 8-bit deltas and no overlap, we have offset = 0xf00 and a maximum
    // encodable value of 0xfff. With 8-bit deltas and a 4-bit overlap, we still have
    // offset = 0xf00 and a maximum encodable value of 0xfff. Even with 12-bit deltas, we have
    // offset = 0 and we can still only represent 0xfff. However with deltaBits = 12 and
    // overlapBits = 4, we can represent offset = 0xf00 and a maximum encodable value of
    // 0xf00 + 0xfff = 0x1eff.
    //
    // It is possible to show that this last example is the worst case, i.e. we do not need to
    // consider increasing deltaBits or overlapBits further.
    int deltaBits = (max(1, 63 - Long.numberOfLeadingZeros(bMax - bMin)) + 3) & ~3;
    int overlapBits = 0;
    if (!canEncode(bMin, bMax, deltaBits, 0, haveExceptions)) {
      if (canEncode(bMin, bMax, deltaBits, 4, haveExceptions)) {
        overlapBits = 4;
      } else {
        Preconditions.checkArgument(deltaBits <= 60);
        deltaBits += 4;
        if (!canEncode(bMin, bMax, deltaBits, 0, haveExceptions)) {
          Preconditions.checkArgument(canEncode(bMin, bMax, deltaBits, 4, haveExceptions));
          overlapBits = 4;
        }
      }
    }

    // Avoid wasting 4 bits of delta when the block size is 1. This reduces the encoding size for
    // single leaf cells by one byte.
    if (values.length() == 1) {
      Preconditions.checkArgument(deltaBits == 4 && overlapBits == 0);
      deltaBits = 8;
    }

    // Now determine the number of bytes needed to encode "offset", given the chosen delta length.
    long maxDelta = bitMask(deltaBits) - (haveExceptions ? BLOCK_SIZE : 0);
    int offsetBits = 0;
    if (UnsignedLongs.compare(bMax, maxDelta) > 0) {
      // At least one byte of offset is required. Round up the minimum offset to the next
      // encodable value, and determine how many bits it has.
      int offsetShift = deltaBits - overlapBits;
      long mask = bitMask(offsetShift);
      long minOffset = (bMax - maxDelta + mask) & ~mask;
      Preconditions.checkArgument(minOffset != 0);
      offsetBits = ((63 - Long.numberOfLeadingZeros(minOffset)) + 1 - offsetShift + 7) & ~7;
      // A 64-bit offset can only be encoded with an overlap of 4 bits.
      if (offsetBits == 64) {
        overlapBits = 4;
      }
    }
    code.set(deltaBits, offsetBits, overlapBits);
  }
}
