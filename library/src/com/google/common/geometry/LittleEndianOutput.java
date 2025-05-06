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
import java.io.OutputStream;

/** Simple utility for writing little endian primitives to a stream. */
public final class LittleEndianOutput {
  private final OutputStream output;

  /** Constructs a little-endian output that writes to the given stream. */
  public LittleEndianOutput(OutputStream output) {
    this.output = output;
  }

  /** Returns the output stream. */
  public OutputStream output() {
    return output;
  }

  /** Writes a byte. */
  public void writeByte(byte value) throws IOException {
    output.write((int) value);
  }

  public void writeBytes(byte[] bytes) throws IOException {
    output.write(bytes);
  }

  /** Writes a little-endian signed integer. */
  public void writeInt(int value) throws IOException {
    output.write(value & 0xFF);
    output.write((value >> 8) & 0xFF);
    output.write((value >> 16) & 0xFF);
    output.write((value >> 24) & 0xFF);
  }

  /** Writes a little-endian signed long. */
  public void writeLong(long value) throws IOException {
    writeLong(output, value);
  }

  /** Writes a little-endian signed long to 'output'. */
  public static void writeLong(OutputStream output, long value) throws IOException {
    output.write((int) (value & 0xFF));
    output.write((int) (value >> 8) & 0xFF);
    output.write((int) (value >> 16) & 0xFF);
    output.write((int) (value >> 24) & 0xFF);
    output.write((int) (value >> 32) & 0xFF);
    output.write((int) (value >> 40) & 0xFF);
    output.write((int) (value >> 48) & 0xFF);
    output.write((int) (value >> 56) & 0xFF);
  }

  /** Writes a little-endian IEEE754 32-bit float. */
  public void writeFloat(float value) throws IOException {
    writeInt(Float.floatToIntBits(value));
  }

  /** Writes a little-endian IEEE754 64-bit double. */
  public void writeDouble(double value) throws IOException {
    writeDouble(output, value);
  }

  /** Writes a little-endian IEEE754 64-bit double to 'output'. */
  public static void writeDouble(OutputStream output, double value) throws IOException {
    writeLong(output, Double.doubleToLongBits(value));
  }

  /**
   * Writes a signed integer using variable encoding with {@link #writeVarint64(long)}.
   *
   * @throws IOException if past end of input or error in underlying stream
   */
  public void writeVarint32(int value) throws IOException {
    writeVarint64(value);
  }

  /**
   * Writes a signed long using variable encoding with {@link
   * EncodedInts#writeVarint64(OutputStream, long)}.
   *
   * @throws IOException if past end of input or error in underlying stream
   */
  public void writeVarint64(long value) throws IOException {
    EncodedInts.writeVarint64(output, value);
  }

  /** Closes the underlying output stream. */
  public void close() throws IOException {
    output.close();
  }
}
