/*
 * Copyright 2014 Google Inc.
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

/**
 * An error code and text string describing the first error encountered during a validation process.
 */
// Package private for now, pending further usage.
@GwtCompatible
class S2Error {
  public enum Code {
    // Error codes that apply to more than one type of geometry:
    /** No problems detected. */
    NO_ERROR(0),
    /** Vertex is not unit length. */
    NOT_UNIT_LENGTH(1),
    /** There are two identical vertices. */
    DUPLICATE_VERTICES(2),
    /** There are two antipodal vertices. */
    ANTIPODAL_VERTICES(3),

    // Error codes that only apply to certain geometric objects:
    /** Loop with fewer than 3 vertices. */
    LOOP_NOT_ENOUGH_VERTICES(100),
    /** Loop has a self-intersection. */
    LOOP_SELF_INTERSECTION(101),

    /** Two polygon loops share an edge. */
    POLYGON_LOOPS_SHARE_EDGE(200),
    /** Two polygon loops cross. */
    POLYGON_LOOPS_CROSS(201),
    /** Polygon has an empty loop. */
    POLYGON_EMPTY_LOOP(202),
    /** Non-full polygon has a full loop. */
    POLYGON_EXCESS_FULL_LOOP(203);

    private int code;

    private Code(int code) {
      this.code = code;
    }

    /** Returns the numeric value of this error code. */
    public int code() {
      return code;
    }
  }

  private Code code = Code.NO_ERROR;
  private String text = "";

  /**
   * Sets the error code and text description; the description is formatted according to the rules
   * defined in {@link String#format(String, Object...)}. This method may be called more than once,
   * so that various layers may surround 
   */
  public void init(Code code, String format, Object ... args) {
    this.code = code;
    this.text = Platform.formatString(format, args);
  }

  /** Returns the code of this error. */
  public Code code() {
    return code;
  }

  /** Returns the text string. */
  public String text() {
    return text;
  }
}
