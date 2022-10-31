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
import com.google.common.base.Preconditions;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.ImmutableLongArray;
import com.google.common.primitives.Ints;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.nio.ByteBuffer;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;

/** A set of interfaces for describing primitive arrays. */
@JsType
public final class PrimitiveArrays {
  private PrimitiveArrays() {}

  /**
   * An array of {@code byte}s.
   *
   * <p>Implementations will be thread-safe if the underlying data is not mutated. Users should
   * ensure the underlying data is not mutated in order to get predictable behaviour. Any buffering
   * should be done internally.
   *
   * <p>Implementations may support arrays &gt; 2GB in size like so:
   *
   * <pre>{@code
   * new Bytes() {
   *   byte get(long position) {
   *     if (position < b1.length) {
   *       return b1[Ints.checkedCast(position)];
   *     }
   *     return b2[Ints.checkedCast(position - b2.length)];
   *   }
   *   long length() { return b1.length + b2.length; }
   * }
   * }</pre>
   */
  @JsType
  public interface Bytes {
    /**
     * Returns the {@code byte} at position {@code position}.
     *
     * <p>Throws an {@link IndexOutOfBoundsException} if the absolute get on the underlying
     * implementation fails.
     */
    byte get(long position);

    /** Returns the length of this array. */
    long length();

    /** Returns a {@link Cursor} with the given {@code position} and {@code limit}. */
    default Cursor cursor(long position, long limit) {
      Preconditions.checkArgument(position <= limit && position <= length());
      return new Cursor(position, limit);
    }

    /**
     * Returns a {@link Cursor} with the given {@code position}.
     *
     * <p>The {@code limit} of the returned cursor is the {@link #length()} of this array.
     */
    @JsIgnore // overload.
    default Cursor cursor(long position) {
      return cursor(position, length());
    }

    /**
     * Returns a {@link Cursor}.
     *
     * <p>The {@code position} of the returned cursor is 0, and the {@code limit} is the {@link
     * #length()} of this array.
     */
    @JsIgnore // overload.
    default Cursor cursor() {
      return cursor(0);
    }

    /**
     * Returns a {@link Bytes} wrapping {@code buffer}.
     *
     * <p>The returned array starts from index 0 of buffer, and its length is {@code
     * buffer.limit()}.
     */
    // TODO(torrey): Consider removing ByteBuffer from PrimitiveArrays to avoid J2CL problems.
    // Provide another implementation for J2CL?
    @GwtIncompatible("ByteBuffer") // JsIgnore is insufficient.
    static Bytes fromByteBuffer(ByteBuffer buffer) {
      // TODO(user): Buffer positions > 0 trip bugs in the various methods. Exclude this
      // case.
      Preconditions.checkState(buffer.position() == 0);
      return new Bytes() {
        @Override
        public byte get(long position) {
          return buffer.get(Ints.checkedCast(position));
        }

        @Override
        public long length() {
          return (long) buffer.limit();
        }
      };
    }

    /** Returns a {@link Bytes} wrapping {@code bytes}. */
    static Bytes fromByteArray(byte[] bytes) {
      return new Bytes() {
        @Override
        public byte get(long position) {
          return bytes[Ints.checkedCast(position)];
        }

        @Override
        public long length() {
          return bytes.length;
        }
      };
    }

    /** Returns an {@link InputStream} wrapping this array starting at {@code offset}. */
    @JsIgnore // InputStream is not usable by Javascript.
    default InputStream toInputStream(long offset) {
      Preconditions.checkArgument(offset >= 0 && offset <= length());
      return new InputStream() {
        long position = offset;

        @Override
        public int read() {
          if (position == length()) {
            return -1;
          }
          return get(position++) & 0xFF;
        }
      };
    }

    /**
     * Returns an {@link InputStream} wrapping this array starting at {@code cursor.position}.
     *
     * <p>{@code cursor.position} is incremented for each byte read from the returned {@link
     * InputStream}.
     */
    @JsIgnore // InputStream is not usable by Javascript.
    default InputStream toInputStream(Cursor cursor) {
      Preconditions.checkArgument(cursor.position >= 0 && cursor.position <= cursor.limit);
      Preconditions.checkArgument(cursor.remaining() <= length());
      return new InputStream() {
        @Override
        public int read() {
          if (cursor.position == cursor.limit) {
            return -1;
          }
          return get(cursor.position++) & 0xFF;
        }
      };
    }

    /** Returns an {@link InputStream} wrapping this array starting at the 0th byte. */
    @JsIgnore // InputStream is not usable by Javascript.
    default InputStream toInputStream() {
      return toInputStream(0);
    }

    /** Writes this array to {@code output}. */
    @JsIgnore // OutputStream is not usable by Javascript.
    default void writeTo(OutputStream output) throws IOException {
      for (long i = 0; i < length(); i++) {
        output.write(get(i));
      }
    }

    /** Returns a byte at {@code cursor.position} and updates the position to the next byte. */
    default byte readByte(Cursor cursor) {
      return get(cursor.position++);
    }

    /**
     * Returns a unsigned integer consisting of {@code numBytes} bytes read from this array at
     * {@code cursor.position} in little-endian format as an unsigned 64-bit integer.
     *
     * <p>{@code cursor.position} is updated to the index of the first byte following the varint64.
     */
    default long readVarint64(Cursor cursor) {
      long result = 0;
      for (int shift = 0; shift < 64; shift += 7) {
        byte b = get(cursor.position++);
        result |= (long) (b & 0x7F) << shift;
        if ((b & 0x80) == 0) {
          return result;
        }
      }
      throw new IllegalArgumentException("Malformed varint.");
    }

    /**
     * Same as {@link #readVarint64(Cursor)}, but throws an {@link IllegalArgumentException} if the
     * read varint64 is greater than {@link Integer#MAX_VALUE}.
     */
    default int readVarint32(Cursor cursor) {
      return Ints.checkedCast(readVarint64(cursor));
    }

    /**
     * Returns a unsigned integer consisting of {@code numBytes} bytes read from this array at
     * {@code cursor.position} in little-endian format as an unsigned 64-bit integer.
     *
     * <p>{@code cursor.position} is updated to the index of the first byte following the uint.
     *
     * <p>This method is not compatible with {@link #readVarint64(Cursor)}.
     */
    default long readUintWithLength(Cursor cursor, int numBytes) {
      long result = readUintWithLength(cursor.position, numBytes);
      cursor.position += numBytes;
      return result;
    }

    /** Same as {@link #readUintWithLength(Cursor, int)}, but does not require a {@link Cursor}. */
    @JsIgnore // No method overloading in J2CL. Use readUintWithLength(Cursor, int).
    default long readUintWithLength(long position, int numBytes) {
      long x = 0;
      for (int i = 0; i < numBytes; i++) {
        x += (get(position++) & 0xffL) << (8 * i);
      }
      return x;
    }

    /** Reads a little endian long from the current cursor position. */
    default long readLittleEndianLong(Cursor cursor) {
      return readUintWithLength(cursor, com.google.common.primitives.Longs.BYTES);
    }

    /** Returns a little-endian double read from this array at {@code position}. */
    default double readLittleEndianDouble(long position) {
      return Double.longBitsToDouble(readUintWithLength(position, Doubles.BYTES));
    }
  }

  /**
   * An array of {@code long}s.
   *
   * <p>Implementations will be thread-safe if the underlying data is not mutated. Users should
   * ensure the underlying data is not mutated in order to get predictable behaviour. Any buffering
   * should be done internally.
   */
  @JsType
  public interface Longs {
    /**
     * Returns the {@code long} at position {@code position}.
     *
     * <p>Throws an {@link IndexOutOfBoundsException} if the absolute get on the underlying
     * implementation fails.
     */
    long get(int position);

    /** Returns the length of this array. */
    int length();

    /** Returns a {@link Longs} wrapping {@code immutableLongArray}. */
    @JsIgnore
    static Longs fromImmutableLongArray(ImmutableLongArray immutableLongArray) {
      return new Longs() {
        @Override
        public long get(int position) {
          return immutableLongArray.get(position);
        }

        @Override
        public int length() {
          return immutableLongArray.length();
        }
      };
    }

    /**
     * Decodes and returns this array as an {@code int[]}.
     *
     * <p>Throws an {@link IllegalArgumentException} if any value in this array is < {@link
     * Integer#MIN_VALUE} or > {@link Integer#MAX_VALUE}.
     */
    default int[] toIntArray() {
      int[] result = new int[length()];
      for (int i = 0; i < result.length; i++) {
        result[i] = Ints.checkedCast(get(i));
      }
      return result;
    }
  }

  /** A cursor storing a position and a limit. */
  @JsType
  public static class Cursor {
    public long position;
    public long limit;

    Cursor(long position, long limit) {
      Preconditions.checkArgument(position >= 0);
      Preconditions.checkArgument(position <= limit);
      this.position = position;
      this.limit = limit;
    }

    /** Returns the number of remaining elements ({@code limit - position}). */
    public long remaining() {
      return limit - position;
    }
  }
}
