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

import com.google.common.annotations.GwtCompatible;
import java.io.IOException;
import java.io.OutputStream;

/**
 * An interface for encoding and decoding values.
 *
 * <p>This is one of several helper classes that allow complex data structures to be initialized
 * from an encoded format in constant time and then decoded on demand. This can be a big performance
 * advantage when only a small part of the data structure is actually used.
 */
@GwtCompatible
public interface S2Coder<T> {

  /** Encodes {@code value} to {@code output}. */
  void encode(T value, OutputStream output) throws IOException;

  /**
   * Decodes a value of type {@link T} from {@code data} starting at {@code cursor.position}. {@code
   * cursor.position} is updated to the position of the first byte in {@code data} following the
   * encoded value.
   */
  T decode(PrimitiveArrays.Bytes data, PrimitiveArrays.Cursor cursor);
}
