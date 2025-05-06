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

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;

/** An abstract class that extends an OutputStream but provides a size() and writeTo() method. */
public abstract class InMemoryOutputStream extends OutputStream {

  public abstract long size();

  public abstract void writeTo(OutputStream outputStream) throws IOException;

  /** ByteArrayScratchOutputStream wraps a ByteArrayOutputStream. */
  static class ByteArrayInMemoryOutputStream extends InMemoryOutputStream {
    ByteArrayOutputStream backingStream;

    ByteArrayInMemoryOutputStream(ByteArrayOutputStream backingStream) {
      this.backingStream = backingStream;
    }

    ByteArrayInMemoryOutputStream() {
      this.backingStream = new ByteArrayOutputStream();
    }

    @Override
    public long size() {
      return backingStream.size();
    }

    @Override
    public void writeTo(OutputStream outputStream) throws IOException {
      backingStream.writeTo(outputStream);
    }

    @Override
    public void write(int b) throws IOException {
      backingStream.write(b);
    }

    @Override
    public void write(byte[] b) throws IOException {
      backingStream.write(b);
    }

    @Override
    public void write(byte[] b, int off, int len) throws IOException {
      backingStream.write(b, off, len);
    }
  }
}
