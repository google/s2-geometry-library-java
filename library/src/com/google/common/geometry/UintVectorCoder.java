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

import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.PrimitiveArrays.Longs;
import com.google.common.primitives.Ints;
import java.io.IOException;
import java.io.OutputStream;

/**
 * An encoder/decoder of {@link Longs}. Either uint64 or uint32 values are supported. Decoding is
 * on-demand, so {@link S2Coder#isLazy()} is true.
 */
public class UintVectorCoder implements S2Coder<Longs> {
  /** An instance of an {@code UintVectorCoder} which encodes/decodes {@code uint32}s. */
  public static final UintVectorCoder UINT32 = new UintVectorCoder(Ints.BYTES);

  /** An instance of an {@code UintVectorCoder} which encodes/decodes {@code uint64}s. */
  public static final UintVectorCoder UINT64 =
      new UintVectorCoder(com.google.common.primitives.Longs.BYTES);

  private final int typeBytes;

  private UintVectorCoder(int typeBytes) {
    this.typeBytes = typeBytes;
  }

  /** Encodes the given {@code Longs} into the given OutputStream. */
  @Override
  public void encode(Longs values, OutputStream output) throws IOException {
    // The encoding format is as follows:
    //
    //   totalBytes (varint64): (values.size() * typeBytes) | (bytesPerWord - 1)
    //   array of values.size() elements [bytesPerWord bytes each]
    //
    // bytesPerWord must be >= 0 so we can encode it in (log2(typeBytes) - 1) bits.

    // oneBits = 1 ensures that bytesPerWord is at least 1.
    long oneBits = 1;
    for (int i = 0; i < values.length(); i++) {
      oneBits |= values.get(i);
    }

    // bytesPerWord is the minimum number of bytes required to encode the largest value in values.
    // It is computed by dividing the minimum number of bits required to represent the largest
    // integer in values by 8 (the division by 8 is the unsigned right shift by 3 bits).
    //
    // In the expression below, (63 - Long.numberOfLeadingZeros(oneBits)) is equivalent to
    // floor(log2(oneBits)). oneBits must be at least 1, so the largest value of
    // Long.numberOfLeadingZeros(oneBits) that is possible is 63.
    //
    // Examples:
    // - oneBits = ~0L: The number of leading 0s in oneBits is 0.
    //   - ((63 - 0) >>> 3) + 1 == 8 bytes per word.
    // - oneBits = 4321L: The number of leading 0s in oneBits is 51.
    //   - ((63 - 51) >>> 3) + 1 == 2 bytes per word.
    // - oneBits = 1L: The number of leading 0s in oneBits is 63.
    //   - ((63 - 63) >>> 3) + 1 == 1 byte per word.
    int bytesPerWord = ((63 - Long.numberOfLeadingZeros(oneBits)) >>> 3) + 1;

    // Since totalBytes must be a multiple of typeBytes, and bytesPerWord must be <= totalBytes,
    // (bytesPerWord - 1) can be encoded in the last few bits of totalBytes (e.g., if this is a
    // uint64 vector, typeBytes is 8, and bytesPerWord can be at most 8).
    //
    // For example, if typeBytes were 4, then any value of (values.length() * typeBytes) leaves us
    // the last 2 bits of totalBytes to encode the number of bytes in each word. Since there are
    // 2 bits to work with, and the largest possible bytesPerWord (4) requires 3 bits to encode, we
    // subtract 1 from bytesPerWord so the data fits in 2 bits. Note that this only works because
    // bytesPerWord cannot be 0 and because typeBytes is a power of 2.
    long totalBytes = ((long) values.length() * typeBytes) | (bytesPerWord - 1);
    EncodedInts.writeVarint64(output, totalBytes);
    for (int i = 0; i < values.length(); i++) {
      EncodedInts.encodeUintWithLength(output, values.get(i), bytesPerWord);
    }
  }

  /**
   * Returns a Longs implementation that on-demand-decodes long values from the underlying data,
   * starting at {@code cursor.position}. {@code cursor.position} is updated to the position of the
   * first byte in {@code data} following the encoded values.
   */
  @Override
  public Longs decode(Bytes data, Cursor cursor) throws IOException {
    // See encode for documentation on the encoding format.
    long totalBytes;
    int size;
    int bytesPerWord;

    try {
      // readVarint64 throws ArrayIndexOutOfBoundsException if data is too short.
      totalBytes = data.readVarint64(cursor);
      if (totalBytes < 0) {
        throw new IOException("Invalid input data, totalBytes = " + totalBytes);
      }
      // checkedCast throws IllegalArgumentException if int would overflow.
      size = Ints.checkedCast(totalBytes / typeBytes);
      bytesPerWord = Ints.checkedCast((totalBytes & (typeBytes - 1)) + 1);
    } catch (ArrayIndexOutOfBoundsException | IllegalArgumentException e) {
      throw new IOException("Input data invalid or too short.", e);
    }
    long offset = cursor.position;

    // Update the position to after these Longs. Position calculations must be 64 bit.
    cursor.position += (long) size * (long) bytesPerWord;

    // Check that the Longs we're going to return won't read past the end of 'data'.
    if (cursor.position > data.length()) {
      throw new IOException(
          Platform.formatString(
              "Decoding from 'data' with length %s bytes, but %s bytes are required.",
              data.length(), cursor.position));
    }

    return new Longs() {
      @Override
      public long get(int position) {
        return data.readUintWithLength(
            offset + (long) position * (long) bytesPerWord, bytesPerWord);
      }

      @Override
      public int length() {
        return size;
      }
    };
  }

  @Override
  public boolean isLazy() {
    return true;
  }
}
