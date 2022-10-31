/*
 * Copyright 2016 Google Inc.
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

import java.io.IOException;
import java.io.InputStream;

/** Simple utility for reading little endian primitives from a stream. */
public final class LittleEndianInput {
  private final InputStream input;

  /** Constructs a little-endian input that reads from the given stream. */
  public LittleEndianInput(InputStream input) {
    this.input = input;
  }

  /**
   * Reads a byte.
   *
   * @throws IOException if {@code input.read()} throws an {@code IOException} or returns -1 (EOF).
   */
  public byte readByte() throws IOException {
    return InputStreams.readByte(input);
  }

  /**
   * Reads a fixed size of bytes from the input.
   *
   * @param size the number of bytes to read.
   * @throws IOException if past end of input or error in underlying stream
   */
  public byte[] readBytes(final int size) throws IOException {
    byte[] result = new byte[size];
    int numRead = input.read(result);
    if (numRead < size) {
      throw new IOException("EOF");
    }
    return result;
  }

  /**
   * Reads a little-endian signed integer.
   *
   * @throws IOException if past end of input or error in underlying stream
   */
  public int readInt() throws IOException {
    return (readByte() & 0xFF)
        | ((readByte() & 0xFF) << 8)
        | ((readByte() & 0xFF) << 16)
        | ((readByte() & 0xFF) << 24);
  }

  /**
   * Reads a little-endian signed long.
   *
   * @throws IOException if past end of input or error in underlying stream
   */
  public long readLong() throws IOException {
    return readLong(input);
  }

  /**
   * Reads a little-endian signed long.
   *
   * @throws IOException if past end of input
   */
  public static long readLong(InputStream input) throws IOException {
    return (InputStreams.readByte(input) & 0xFFL)
        | ((InputStreams.readByte(input) & 0xFFL) << 8)
        | ((InputStreams.readByte(input) & 0xFFL) << 16)
        | ((InputStreams.readByte(input) & 0xFFL) << 24)
        | ((InputStreams.readByte(input) & 0xFFL) << 32)
        | ((InputStreams.readByte(input) & 0xFFL) << 40)
        | ((InputStreams.readByte(input) & 0xFFL) << 48)
        | ((InputStreams.readByte(input) & 0xFFL) << 56);
  }

  /**
   * Reads a little-endian IEEE754 32-bit float.
   *
   * @throws IOException if past end of input or error in underlying stream
   */
  public float readFloat() throws IOException {
    return Float.intBitsToFloat(readInt());
  }

  /**
   * Reads a little-endian IEEE754 64-bit double.
   *
   * @throws IOException if past end of input or error in underlying stream
   */
  public double readDouble() throws IOException {
    return Double.longBitsToDouble(readLong());
  }

  /**
   * Reads a variable-encoded signed integer with {@link #readVarint64()}.
   *
   * @throws IOException if past end of input or error in underlying stream
   */
  public int readVarint32() throws IOException {
    return (int) readVarint64();
  }

  /**
   * Reads a variable-encoded signed long with {@link EncodedInts#readVarint64(InputStream)}
   *
   * @throws IOException if past end of input or error in underlying stream
   */
  public long readVarint64() throws IOException {
    return EncodedInts.readVarint64(input);
  }

  /** Closes the underlying stream. */
  public void close() throws IOException {
    input.close();
  }
}
