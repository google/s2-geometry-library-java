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

import static java.lang.Math.max;
import static java.lang.Math.min;

import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.PrimitiveArrays.Longs;
import com.google.common.primitives.UnsignedLongs;
import java.io.IOException;
import java.io.OutputStream;
import java.util.List;

/** An encoder/decoder of Lists of {@link S2CellId}s. */
class S2CellIdVectorCoder implements S2Coder<List<S2CellId>> {

  /** An instance of an {@code S2CellIdVectorCoder}. */
  static final S2CellIdVectorCoder INSTANCE = new S2CellIdVectorCoder();

  /**
   * Encodes the given list of S2CellId values into the provided OutputStream.
   *
   * <p>The encoding format is as follows:
   * <ul>
   * <li> byte 0, bits 0-2: baseBytes
   * <li> byte 0, bits 3-7: shift
   * <li> byte 1: extended shift (only written for odd shift {@code >= 5})
   * <li> 0-7 bytes: base
   * <li> values.size() bytes: encoded uint64s of deltas, encoded by
   *      {@link UintVectorCoder.UINT64#encode(Longs, OutputStream)}.
   * </ul>
   * Where {@code base} is the baseBytes most-significant bytes of the minimum S2CellId. It is
   * 0-7 bytes, and is always shifted so its bytes are the most-significant bytes of a uint64
   * (little-endian).
   *
   * <p>{@code shift} is in the range 0-56. shift is odd only if all S2CellIds are at the same
   * level, in which case the bit at position (shift - 1) in base is automatically set to 1.
   *
   * <p>{@code base} (3 bits) and {@code shift} (6 bits) are encoded in either one or two bytes as
   * follows: If (shift <= 4 or shift is even), then 1 byte, else 2 bytes.
   *
   * <p>(shift == 1) means that all S2CellIds are leaf cells, and (shift == 2) means that all
   * S2CellIds are at level 29.
   */
  @Override
  public void encode(List<S2CellId> values, OutputStream output) throws IOException {
    long valuesOr = 0L;
    long valuesAnd = ~0L;
    long valuesMin = ~0L;
    long valuesMax = 0L;

    for (S2CellId cellId : values) {
      valuesOr |= cellId.id();
      valuesAnd &= cellId.id();
      valuesMin = UnsignedLongs.min(valuesMin, cellId.id());
      valuesMax = UnsignedLongs.max(valuesMax, cellId.id());
    }

    long base = 0L;
    // The number of bytes required to encode base.
    int baseBytes = 0;
    int shift = 0;
    // The bit position of the most-significant bit of the largest delta.
    int maxDeltaMsb = 0;

    if (UnsignedLongs.compare(valuesOr, 0) > 0) {
      // We only allow even shift, unless all values have the same low bit (in which case shift is
      // odd and the preceding bit is implicitly on). There is no point in allowing shifts > 56
      // because deltas are encoded in at least 1 byte each.
      shift = min(56, Long.numberOfTrailingZeros(valuesOr) & ~1);
      if ((valuesAnd & (1L << shift)) != 0) {
        // All S2CellIds are at the same level.
        shift++;
      }

      // base consists of the baseBytes most-significant bytes of the minimum S2CellId. We consider
      // all possible values of baseBytes (0-7) and choose the one that minimizes the total encoding
      // size.
      // The best encoding size so far.
      int minBytes = -1;
      for (int tmpBaseBytes = 0; tmpBaseBytes < 8; tmpBaseBytes++) {
        // The base value being tested (first tmpBaseBytes of valuesMin).
        long tmpBase = valuesMin & ~(~0L >>> (8 * tmpBaseBytes));
        // The most-significant bit position of the largest delta (or zero if there are no deltas
        // [i.e., if values.size == 0]).
        int tmpMaxDeltaMsb =
            max(0, 63 - Long.numberOfLeadingZeros((valuesMax - tmpBase) >>> shift));
        // The total size of the variable portion of the encoding.
        int candidateBytes = tmpBaseBytes + values.size() * ((tmpMaxDeltaMsb >> 3) + 1);

        if (UnsignedLongs.compare(candidateBytes, minBytes) < 0) {
          base = tmpBase;
          baseBytes = tmpBaseBytes;
          maxDeltaMsb = tmpMaxDeltaMsb;
          minBytes = candidateBytes;
        }
      }
      // It takes one extra byte to encode odd shifts (i.e., the case where all S2CellIds are at the
      // same level), so we check whether we can get the same encoding size per delta using an even
      // shift.
      if (((shift & 1) != 0) && (maxDeltaMsb & 7) != 7) {
        shift--;
      }
    }
    assert shift <= 56;

    // shift and baseBytes are encoded in 1 or 2 bytes.
    // shiftCode is 5 bits, values:
    // - <= 28 represent even shifts in the range 0-56.
    // - 29, 30 represent odd shifts 1 and 3.
    // - 31 indicates that the shift is odd and encoded in the next byte.
    int shiftCode = shift >> 1;
    if ((shift & 1) != 0) {
      shiftCode = min(31, shiftCode + 29);
    }
    output.write((byte) ((shiftCode << 3) | baseBytes));
    if (shiftCode == 31) {
      // shift is always odd, so 3 bits unused.
      output.write((byte) (shift >> 1));
    }

    // Encode the baseBytes most-significant bytes of base.
    long baseCode = base >>> (64 - 8 * max(1, baseBytes));
    EncodedInts.encodeUintWithLength(output, baseCode, baseBytes);

    // Encode the vector of deltas.
    long tmpBase = base;
    long tmpShift = shift;
    UintVectorCoder.UINT64.encode(
        new Longs() {
          @Override
          public long get(int position) {
            return (values.get(position).id() - tmpBase) >>> tmpShift;
          }

          @Override
          public int length() {
            return values.size();
          }
        },
        output);
  }

  /**
   * Returns a {@link S2CellIdVector} which wraps the provided Bytes to decode S2CellIds on demand,
   * starting at the given {@code cursor#position}. {@code cursor#position} is updated to the
   * position of the first byte in {@code data} following the encoded list of S2CellIds. See {@link
   * #encode(List, OutputStream)} for documentation of the encoding format.
   *
   * <p>Throws IOException if the input is invalid.
   */
  @Override
  public S2CellIdVector decode(Bytes data, Cursor cursor) throws IOException {
    try {
      return decodeInternal(data, cursor);
    } catch (IndexOutOfBoundsException e) {
      throw new IOException("Insufficient or invalid input bytes: ", e);
    }
  }

  /**
   * As above, but may also throw IndexOutOfBoundsException if the input data is too short.
   */
  private S2CellIdVector decodeInternal(Bytes data, Cursor cursor) throws IOException {
    // Invert the encoding of (shiftCode, baseBytes).
    int shiftCodeBaseBytes = data.get(cursor.position++) & 0xff;
    int shiftCode = shiftCodeBaseBytes >> 3;
    if (shiftCode == 31) {
      shiftCode = 29 + (data.get(cursor.position++) & 0xff);
    }

    // Decode the baseBytes most-significant bytes of base.
    int baseBytes = shiftCodeBaseBytes & 7;
    long base = data.readUintWithLength(cursor, baseBytes);
    base <<= 64 - 8 * max(1, baseBytes);

    // Invert the encoding of shiftCode.
    long shift;
    if (shiftCode >= 29) {
      shift = 2L * (shiftCode - 29) + 1;
      base |= 1L << (shift - 1);
    } else {
      shift = 2L * shiftCode;
    }

    long tmpBase = base;
    Longs deltas = UintVectorCoder.UINT64.decode(data, cursor);
    return new S2CellIdVector() {
      @Override
      public int size() {
        return deltas.length();
      }

      @Override
      public S2CellId get(int index) {
        return new S2CellId((deltas.get(index) << shift) + tmpBase);
      }

      @Override
      int lowerBound(S2CellId target) {
        if (UnsignedLongs.compare(target.id(), tmpBase) <= 0) {
          return 0;
        }
        if (target.greaterOrEquals(S2CellId.end(S2CellId.MAX_LEVEL))) {
          return size();
        }
        int low = 0;
        int high = deltas.length();
        long needle = (target.id() - tmpBase + (1L << shift) - 1) >>> shift;

        // Binary search for the index of the first element in deltas that is >= needle.
        while (low < high) {
          int mid = (low + high) >> 1;
          long value = deltas.get(mid);
          if (UnsignedLongs.compare(value, needle) < 0) {
            low = mid + 1;
          } else {
            high = mid;
          }
        }
        return low;
      }
    };
  }

  @Override
  public boolean isLazy() {
    return true;
  }
}
