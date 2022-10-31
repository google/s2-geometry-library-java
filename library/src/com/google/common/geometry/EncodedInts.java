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

import com.google.common.primitives.UnsignedInts;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

/** Utilities for encoding and decoding integers. */
public class EncodedInts {
  private EncodedInts() {}

  /**
   * Reads a variable-encoded signed long.
   *
   * <p>Note that if you frequently read/write negative numbers, you should consider zigzag-encoding
   * your values before storing them as varints. See {@link EncodedInts#encodeZigZag32} and {@link
   * #decodeZigZag32(int)}.
   *
   * @throws IOException if {@code input.read()} throws an {@code IOException} or returns -1 (EOF),
   *     or if the variable-encoded signed long is malformed.
   */
  public static long readVarint64(InputStream input) throws IOException {
    long result = 0;
    for (int shift = 0; shift < 64; shift += 7) {
      final byte b = InputStreams.readByte(input);
      result |= (long) (b & 0x7F) << shift;
      if ((b & 0x80) == 0) {
        return result;
      }
    }
    throw new IOException("Malformed varint.");
  }

  /** Returns the size in bytes of {@code value} when encoded by {@link #writeVarint64}. */
  public static int varIntSize(long value) {
    int bytes = 0;
    do {
      bytes++;
      value >>>= 7;
    } while (value != 0);
    return bytes;
  }

  /**
   * Writes a signed long using variable encoding.
   *
   * <p>Note that if you frequently read/write negative numbers, you should consider zigzag-encoding
   * your values before storing them as varints. See {@link EncodedInts#encodeZigZag32} and {@link
   * #decodeZigZag32(int)}.
   *
   * @throws IOException if {@code output.write(int)} throws an {@link IOException}.
   */
  public static void writeVarint64(OutputStream output, long value) throws IOException {
    while (true) {
      if ((value & ~0x7FL) == 0) {
        output.write((byte) value);
        return;
      } else {
        output.write((byte) (((int) value & 0x7F) | 0x80));
        value >>>= 7;
      }
    }
  }

  /**
   * Decodes a unsigned integer consisting of {@code bytesPerWord} bytes from {@code supplier} in
   * little-endian format as an unsigned 64-bit integer.
   *
   * <p>This method is not compatible with {@link #readVarint64(InputStream)} or {@link
   * #writeVarint64(OutputStream, long)}.
   *
   * @throws IOException if {@code input.read()} throws an {@code IOException} or returns -1 (EOF).
   */
  public static long decodeUintWithLength(InputStream input, int bytesPerWord) throws IOException {
    long x = 0;
    for (int i = 0; i < bytesPerWord; i++) {
      x += (InputStreams.readByte(input) & 0xffL) << (8 * i);
    }
    return x;
  }

  /**
   * Encodes an unsigned integer to {@code consumer} in little-endian format using {@code
   * bytesPerWord} bytes. (The client must ensure that the encoder's buffer is large enough).
   *
   * <p>This method is not compatible with {@link #readVarint64(InputStream)} or {@link
   * #writeVarint64(OutputStream, long)}.
   *
   * @throws IOException if {@code output.write(int)} throws an {@link IOException}.
   */
  public static void encodeUintWithLength(OutputStream output, long value, int bytesPerWord)
      throws IOException {
    while (--bytesPerWord >= 0) {
      output.write((byte) value);
      value >>>= 8;
    }
    assert value == 0;
  }

  /**
   * Encode a ZigZag-encoded 32-bit value. ZigZag encodes signed integers into values that can be
   * efficiently encoded with varint. (Otherwise, negative values must be sign-extended to 64 bits
   * to be varint encoded, thus always taking 10 bytes on the wire.)
   *
   * @param n A signed 32-bit integer.
   * @return An unsigned 32-bit integer, stored in a signed int because Java has no explicit
   *     unsigned support.
   */
  public static int encodeZigZag32(final int n) {
    // Note:  the right-shift must be arithmetic
    return (n << 1) ^ (n >> 31);
  }

  /**
   * Encode a ZigZag-encoded 64-bit value. ZigZag encodes signed integers into values that can be
   * efficiently encoded with varint. (Otherwise, negative values must be sign-extended to 64 bits
   * to be varint encoded, thus always taking 10 bytes on the wire.)
   *
   * @param n A signed 64-bit integer.
   * @return An unsigned 64-bit integer, stored in a signed int because Java has no explicit
   *     unsigned support.
   */
  public static long encodeZigZag64(final long n) {
    // Note:  the right-shift must be arithmetic
    return (n << 1) ^ (n >> 63);
  }

  /**
   * Decode a ZigZag-encoded 32-bit signed value. ZigZag encodes signed integers into values that
   * can be efficiently encoded with varint. (Otherwise, negative values must be sign-extended to 64
   * bits to be varint encoded, thus always taking 10 bytes on the wire.)
   *
   * @param n A 32-bit integer, stored in a signed int because Java has no explicit unsigned
   *     support.
   * @return A signed 32-bit integer.
   */
  public static int decodeZigZag32(final int n) {
    return (n >>> 1) ^ -(n & 1);
  }

  /**
   * Decode a ZigZag-encoded 64-bit signed value. ZigZag encodes signed integers into values that
   * can be efficiently encoded with varint. (Otherwise, negative values must be sign-extended to 64
   * bits to be varint encoded, thus always taking 10 bytes on the wire.)
   *
   * @param n A 64-bit integer, stored in a signed long because Java has no explicit unsigned
   *     support.
   * @return A signed 64-bit integer.
   */
  public static long decodeZigZag64(final long n) {
    return (n >>> 1) ^ -(n & 1);
  }

  /**
   * Returns the interleaving of bits of val1 and val2, where the LSB of val1 is the LSB of the
   * result, and the MSB of val2 is the MSB of the result.
   */
  public static long interleaveBits(int val1, int val2) {
    return insertBlankBits(val1) | (insertBlankBits(val2) << 1);
  }

  /** Returns the first int de-interleaved from the result of {@link #interleaveBits}. */
  public static int deinterleaveBits1(long bits) {
    return removeBlankBits(bits);
  }

  /** Returns the second int de-interleaved from the result of {@link #interleaveBits}. */
  public static int deinterleaveBits2(long bits) {
    return removeBlankBits(bits >>> 1);
  }

  /**
   * Inserts blank bits between the bits of 'value' such that the MSB is blank and the LSB is
   * unchanged.
   */
  private static final long insertBlankBits(int value) {
    long bits = UnsignedInts.toLong(value);
    bits = (bits | (bits << 16)) & 0x0000ffff0000ffffL;
    bits = (bits | (bits << 8)) & 0x00ff00ff00ff00ffL;
    bits = (bits | (bits << 4)) & 0x0f0f0f0f0f0f0f0fL;
    bits = (bits | (bits << 2)) & 0x3333333333333333L;
    bits = (bits | (bits << 1)) & 0x5555555555555555L;
    return bits;
  }

  /** Reverses {@link #insertBlankBits} by extracting the even bits (bit 0, 2, ...). */
  private static int removeBlankBits(long bits) {
    bits &= 0x5555555555555555L;
    bits |= bits >>> 1;
    bits &= 0x3333333333333333L;
    bits |= bits >>> 2;
    bits &= 0x0f0f0f0f0f0f0f0fL;
    bits |= bits >>> 4;
    bits &= 0x00ff00ff00ff00ffL;
    bits |= bits >>> 8;
    bits &= 0x0000ffff0000ffffL;
    bits |= bits >>> 16;
    return (int) bits;
  }

  /**
   * Like {@link #interleaveBits} but interleaves bit pairs rather than individual bits. This format
   * is faster to decode than the fully interleaved format, and produces the same results for our
   * use case.
   *
   * <p>This code is about 10% faster than {@link #interleaveBits}.
   */
  public static long interleaveBitPairs(int val1, int val2) {
    return insertBlankPairs(val1) | (insertBlankPairs(val2) << 2);
  }

  /** Returns the first int de-interleaved from the result of {@link #interleaveBitPairs}. */
  public static int deinterleaveBitPairs1(long pairs) {
    return removeBlankPairs(pairs);
  }

  /** Returns the second int de-interleaved from the result of {@link #interleaveBitPairs}. */
  public static int deinterleaveBitPairs2(long pairs) {
    return removeBlankPairs(pairs >>> 2);
  }

  /** Inserts 00 pairs in between the pairs from 'value'. */
  private static final long insertBlankPairs(int value) {
    long bits = UnsignedInts.toLong(value);
    bits = (bits | (bits << 16)) & 0x0000ffff0000ffffL;
    bits = (bits | (bits << 8)) & 0x00ff00ff00ff00ffL;
    bits = (bits | (bits << 4)) & 0x0f0f0f0f0f0f0f0fL;
    bits = (bits | (bits << 2)) & 0x3333333333333333L;
    return bits;
  }

  /**
   * Reverses {#link #insertBitPairs} by selecting the two LSB bits, dropping the next two,
   * selecting the next two, etc.
   */
  private static int removeBlankPairs(long pairs) {
    pairs &= 0x3333333333333333L;
    pairs |= pairs >>> 2;
    pairs &= 0x0f0f0f0f0f0f0f0fL;
    pairs |= pairs >>> 4;
    pairs &= 0x00ff00ff00ff00ffL;
    pairs |= pairs >>> 8;
    pairs &= 0x0000ffff0000ffffL;
    pairs |= pairs >>> 16;
    return (int) pairs;
  }
}
