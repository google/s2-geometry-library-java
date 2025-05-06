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

import static java.nio.charset.StandardCharsets.UTF_8;

import com.google.common.base.Preconditions;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.function.Supplier;
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

  /** A serializable function from type A to type B that may throw an IOException. */
  @JsType
  public interface SerializableFunction<A, B> extends Serializable {
    B apply(A input) throws IOException;
  }

  /** A coder of an unboxed varint. */
  @SuppressWarnings("CheckedExceptionNotThrown")
  S2Coder<Long> UNBOXED_VARINT =
      new S2Coder<>() {
        @Override
        public void encode(Long value, OutputStream output) throws IOException {
          EncodedInts.writeVarint64(output, value);
        }

        @Override
        public Long decode(Bytes data, Cursor cursor) throws IOException {
          return data.readVarint64(cursor);
        }
      };

  /** An encoder of byte[] that writes the varint length of the byte array before it. */
  S2Coder<byte[]> BYTES =
      new S2Coder<>() {
        @Override
        public void encode(byte[] bytes, OutputStream output) throws IOException {
          EncodedInts.writeVarint64(output, bytes.length);
          output.write(bytes);
        }

        @Override
        public byte[] decode(Bytes data, Cursor cursor) {
          int length = data.readVarint32(cursor);
          Preconditions.checkArgument(length <= cursor.remaining(), "Length too long");
          byte[] b = new byte[length];
          for (int i = 0; i < b.length; i++) {
            b[i] = data.get(cursor.position++);
          }
          return b;
        }

        @Override
        public boolean isLazy() {
          return false;
        }
      };

  /** A delimited string coder that converts to/from UTF8 and delegates to {@link #BYTES}. */
  S2Coder<String> STRING = BYTES.delegating(
      value -> value.getBytes(UTF_8),
      b -> new String(b, UTF_8));

  /** Encodes {@code value} to {@code output}. */
  @JsIgnore // OutputStream is not available to J2CL.
  void encode(T value, OutputStream output) throws IOException;

  @JsIgnore
  default byte[] encode(T value) throws IOException {
    ByteArrayOutputStream bytes = new ByteArrayOutputStream();
    encode(value, bytes);
    return bytes.toByteArray();
  }

  /** As {@link #encode(T)} but wraps any exceptions in IllegalArgumentException. */
  default byte[] unsafeEncode(T value) {
    try {
      return encode(value);
    } catch (IOException e) {
      throw new IllegalArgumentException(e);
    }
  }

  /** As {@link #decode(Bytes)} but wraps any exceptions in IllegalArgumentException. */
  default T unsafeDecode(Bytes data) {
    try {
      return decode(data);
    } catch (IOException e) {
      throw new IllegalArgumentException(e);
    }
  }

  /**
   * Returns a supplier that wraps the given value and serializes itself with this coder, decoding
   * afterward on first access. The resulting supplier is serializable but not thread-safe.
   */
  default SerializableSupplier<T> memoize(T value) {
    return new MemoizeEncoded<>(this, value);
  }

  /** A memoized value that serializes itself in terms of the encoded representation. */
  final class MemoizeEncoded<T> implements SerializableSupplier<T> {
    private final S2Coder<T> coder;
    private transient T decoded;
    private final byte[] encoded;
    MemoizeEncoded(S2Coder<T> coder, T value) {
      this.coder = coder;
      this.decoded = value;
      this.encoded = coder.unsafeEncode(value);
    }
    @Override public T get() {
      if (decoded == null) {
        return decoded = coder.unsafeDecode(Bytes.fromByteArray(encoded));
      }
      return decoded;
    }
  }

  /** A serializable supplier that can be returned by memoize. */
  interface SerializableSupplier<T> extends Supplier<T>, Serializable {}

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

  /**
   * Returns a coder that delegates to this coder via the given encode/decode transforms, which may
   * throw IOExceptions on invalid input.
   */
  default <U> S2Coder<U> delegating(
      SerializableFunction<U, T> encode,
      SerializableFunction<T, U> decode) {
    S2Coder<T> base = this;
    return new S2Coder<U>() {
      @JsIgnore // OutputStream is not available to J2CL.
      @Override
      public void encode(U value, OutputStream output) throws IOException {
        base.encode(encode.apply(value), output);
      }

      @Override
      public U decode(Bytes data, Cursor cursor) throws IOException {
        // Note: When delegation passes through an AbstractList wrapper, and the data is short,
        // AbstractList may convert an IndexOutOfBoundsException into a NoSuchElementException.
        return decode.apply(base.decode(data, cursor));
      }

      @Override
      public boolean isLazy() {
        return base.isLazy();
      }
    };
  }
}
