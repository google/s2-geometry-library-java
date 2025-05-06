/*
 * Copyright 2024 Google Inc.
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

import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;

/** An R2Edge is an mutable edge in two-dimensional space. */
@JsType
public class R2Edge {
  final R2Vector v0 = new R2Vector();
  final R2Vector v1 = new R2Vector();

  /** Creates a new mutable edge with both endpoints initially at (0, 0). */
  public R2Edge() {}

  /** Initializes this edge endpoints to be copies of the given endpoints. */
  @JsIgnore
  public void init(R2Vector v0, R2Vector v1) {
    this.v0.set(v0);
    this.v1.set(v1);
  }

  /** Sets this edge endpoints to be copies of the current endpoints of the given edge. */
  public void init(R2Edge edge) {
    this.v0.set(edge.v0);
    this.v1.set(edge.v1);
  }

  /**
   * Returns true if the current endpoints of this edge have exactly the same values as the current
   * endpoints of the given other edge.
   */
  public boolean isEqualTo(R2Edge other) {
    return v0.isEqualTo(other.v0) && v1.isEqualTo(other.v1);
  }

  /**
   * Mutable objects should not implement equals(). Use {@link #isEqualTo(R2Edge)} to compare the
   * current values of two R2Edges.
   */
  @Override
  public final boolean equals(Object o) {
    throw new UnsupportedOperationException("R2Edge is mutable and does not support equals.");
  }

  /** Mutable objects should not implement hashCode(). */
  @Override
  public final int hashCode() {
    throw new UnsupportedOperationException("R2Edge is mutable and does not support hashCode.");
  }

  @Override
  public String toString() {
    return Platform.formatString("R2Edge(%s, %s)", v0, v1);
  }
}
