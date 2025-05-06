/*
 * Copyright 2023 Google Inc.
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
import com.google.common.primitives.Ints;
import java.nio.ByteBuffer;

/** Static utility methods for handling Java ByteBuffers. */
public class BufferUtils {

  /**
   * Returns a {@link Bytes} wrapping {@code buffer}.
   *
   * <p>The returned array starts from index 0 of buffer, and its length is {@code buffer.limit()}.
   */
  public static Bytes createBytes(ByteBuffer buffer) {
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

  private BufferUtils() {}
}
