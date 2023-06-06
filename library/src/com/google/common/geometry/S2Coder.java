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

import com.google.common.base.Preconditions;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.S2.SerializableFunction;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import java.nio.charset.Charset;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;

/**
 * An interface for encoding and decoding values.
 *
 * <p>This is one of several helper classes that allow complex data structures to be initialized
 * from an encoded format in constant time and then decoded on demand. This can be a big performance
 * advantage when only a small part of the data structure is actually used.
 */
@JsType
public interface S2Coder<T> extends Serializable {
  /** A delimited string coder that converts to/from UTF8, writing a varint length before it. */
  S2Coder<String> STRING =
      new S2Coder<String>() {
        /** The charset to use. Note the StandardCharsets constant is not available in Android. */
        private final Charset charset = Charset.forName("UTF-8");

        @Override
        public void encode(String value, OutputStream output) throws IOException {
          byte[] bytes = value.getBytes(charset);
          EncodedInts.writeVarint64(output, bytes.length);
          output.write(bytes);
        }

        @Override
        public String decode(PrimitiveArrays.Bytes data, Cursor cursor) {
          int length = data.readVarint32(cursor);
          Preconditions.checkArgument(length <= cursor.remaining(), "String too long");
          byte[] b = new byte[length];
          for (int i = 0; i < b.length; i++) {
            b[i] = data.get(cursor.position++);
          }
          return new String(b, charset);
        }

        @Override
        public boolean isLazy() {
          return false;
        }
      };

  /** Encodes {@code value} to {@code output}. */
  @JsIgnore // OutputStream is not available to J2CL.
  void encode(T value, OutputStream output) throws IOException;

  /**
   * Decodes a value of type {@link T} from {@code data} starting at {@code cursor.position}. {@code
   * cursor.position} is updated to the position of the first byte in {@code data} following the
   * encoded value.
   *
   * <p>Warning: If {@link #isLazy()} is true, the S2Coder may keep a reference to the provided
   * 'data' after decode() has returned. Callers should ensure that the bytes given to decode()
   * are unaltered for as long as the returned {@code T} is in use.
   */
  T decode(Bytes data, Cursor cursor) throws IOException;

  /** As {@link #decode(Bytes, Cursor)} but reads from the current position in {@code data}. */
  @JsIgnore // No method overloading in J2CL. Use decode(data, data.cursor()).
  default T decode(Bytes data) throws IOException {
    return decode(data, data.cursor());
  }

  /**
   * Must return true if the implementation of {@link #decode(Bytes)} retains a reference to the
   * buffer it is given. Lazy decoding is usually cheaper at decode time, but may be more expensive
   * if the entire object will be inspected. Callers should ensure that the bytes given to decode()
   * are unaltered for as long as the resulting {@code T} is in use.
   */
  default boolean isLazy() {
    return true;
  }

  /** Returns a coder that delegates to this coder via the given encode/decode transforms. */
  default <U> S2Coder<U> delegating(
      SerializableFunction<U, T> encode,
      SerializableFunction<T, U> decode) {
    S2Coder<T> base = this;
    return new S2Coder<U>() {
      @JsIgnore // OutputStream is not available to J2CL.
      @Override public void encode(U value, OutputStream output) throws IOException {
        base.encode(encode.apply(value), output);
      }

      @Override public U decode(Bytes data, Cursor cursor) throws IOException {
        // Note: When delegation passes through an AbstractList wrapper, and the data is short,
        // AbstractList may convert an IndexOutOfBoundsException into a NoSuchElementException.
        return decode.apply(base.decode(data, cursor));
      }

      @Override public boolean isLazy() {
        return base.isLazy();
      }
    };
  }
}
