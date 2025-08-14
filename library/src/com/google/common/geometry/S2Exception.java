/*
 * Copyright 2025 Google Inc.
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

/**
 * An unchecked exception thrown from the S2 library. S2 methods where an error can occur generally
 * accept an {@link S2Error} parameter. If an error occurs, an {@link S2Error.Code} is assigned,
 * along with some textual description. This matches the C++ implementation of the S2 library.
 *
 * <p>In Java, throwing and catching exceptions for error handling is preferred by some clients. So,
 * for some S2 library methods that accept S2Error parameters, an alternate form that throws
 * S2Exception instead is provided. An S2Exception wraps an S2Error, and provides a convenience
 * method to get the underlying {@link S2Error.Code}. Methods that throw S2Exception are generally
 * named with an "Unsafe" suffix as a warning, because S2Exception is unchecked. For example,
 * {@link S2BooleanOperation#buildUnsafe(S2ShapeIndex, S2ShapeIndex)} is a wrapper for
 * {@link S2BooleanOperation#build(S2ShapeIndex, S2ShapeIndex, S2Error)}.
 */
public class S2Exception extends RuntimeException {

  private final S2Error error;

  /** Creates a new S2Exception wrapping the given S2Error. */
  public S2Exception(S2Error error) {
    this.error = error;
  }

  /** Returns the code of the S2Error wrapped by this S2Exception. */
  public S2Error.Code code() {
    return error.code();
  }

  @Override
  public String getMessage() {
    return error.text();
  }
}
