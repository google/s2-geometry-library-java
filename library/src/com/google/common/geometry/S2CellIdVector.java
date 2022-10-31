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

import java.util.AbstractList;
import jsinterop.annotations.JsType;

/**
 * A list of {@link S2CellId}s, and specialized methods for directly operating on the encoded form.
 */
@JsType
abstract class S2CellIdVector extends AbstractList<S2CellId> {

  /**
   * Returns the index of the first element {@code x} such that {@code (x >= target)}, or {@link
   * #size()} if no such element exists.
   *
   * <p>The list must be sorted into ascending order prior to making this call. If it is not sorted,
   * the results are undefined.
   */
  abstract int lowerBound(S2CellId target);
}
