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
import com.google.common.base.Strings;
import jsinterop.annotations.JsType;

/**
 * An error code and text string describing the first error encountered during a validation process.
 */
@JsType
public class S2Error {
  /** Numeric values for S2 errors. */
  @JsType
  public enum Code {
    /** No problems detected. */
    NO_ERROR(0),

    ////////////////////////////////////////////////////////////////////
    // Generic errors, not specific to geometric objects:

    /** Unknown error. */
    UNKNOWN(1000),
    /** Operation is not implemented. */
    UNIMPLEMENTED(1001),
    /** Argument is out of range. */
    OUT_OF_RANGE(1002),
    /** Invalid argument (other than a range error). */
    INVALID_ARGUMENT(1003),
    /** Object is not in the required state. */
    FAILED_PRECONDITION(1004),
    /** An internal invariant has failed. */
    INTERNAL(1005),
    /** Data loss or corruption. */
    DATA_LOSS(1006),
    /** A resource has been exhausted. */
    RESOURCE_EXHAUSTED(1007),
    /** Operation was cancelled. */
    CANCELLED(1008),

    ////////////////////////////////////////////////////////////////////
    // Error codes that apply to more than one type of geometry:

    /** Vertex is not unit length. */
    NOT_UNIT_LENGTH(1),
    /** There are two identical vertices. */
    DUPLICATE_VERTICES(2),
    /** There are two antipodal vertices. */
    ANTIPODAL_VERTICES(3),
    /** Edges of a chain aren't continuous. */
    NOT_CONTINUOUS(4),
    /** Vertex has value that's inf or NaN. */
    INVALID_VERTEX(5),

    ////////////////////////////////////////////////////////////////////
    // S2Loop errors:

    /** Loop with fewer than 3 vertices. */
    LOOP_NOT_ENOUGH_VERTICES(100),
    /** Loop has a self-intersection. */
    LOOP_SELF_INTERSECTION(101),

    ////////////////////////////////////////////////////////////////////
    // S2Polygon / S2Shape errors:

    /** Two polygon loops share an edge. */
    POLYGON_LOOPS_SHARE_EDGE(200),
    /** Two polygon loops cross. */
    POLYGON_LOOPS_CROSS(201),
    /** Polygon has an empty loop. */
    POLYGON_EMPTY_LOOP(202),
    /** Non-full polygon has a full loop. */
    POLYGON_EXCESS_FULL_LOOP(203),
    /**
     * Inconsistent loop orientations were detected, indicating that the interior is not on the left
     * of all edges.
     */
    POLYGON_INCONSISTENT_LOOP_ORIENTATIONS(204),
    /** Loop depths don't correspond to any valid nesting hierarchy. */
    POLYGON_INVALID_LOOP_DEPTH(205),
    /** Actual polygon nesting does not correspond to the nesting given in the loop depths. */
    POLYGON_INVALID_LOOP_NESTING(206),

    /** Shape dimension isn't valid. */
    INVALID_DIMENSION(207),
    /** Interior split by holes. */
    SPLIT_INTERIOR(208),
    /** Geometry overlaps where it shouldn't. */
    OVERLAPPING_GEOMETRY(209),

    ////////////////////////////////////////////////////////////////////
    // S2Builder errors:

    /**
     * The S2Builder snap function moved a vertex by more than the specified snap radius. Also used
     * if a Double.NaN is encountered on an input vertex.
     */
    BUILDER_SNAP_RADIUS_TOO_SMALL(300),
    /**
     * S2Builder expected all edges to have siblings (as specified by {@link
     * S2Builder.GraphOptions.SiblingPairs#REQUIRE}), but some were missing.
     */
    BUILDER_MISSING_EXPECTED_SIBLING_EDGES(301),
    /**
     * S2Builder found an unexpected degenerate edge. For example, S2Builder.getLeftTurnMap() does
     * not support degenerate edges.
     */
    BUILDER_UNEXPECTED_DEGENERATE_EDGE(302),
    /**
     * S2Builder found a vertex with (indegree != outdegree), which means that the given edges
     * cannot be assembled into loops.
     */
    BUILDER_EDGES_DO_NOT_FORM_LOOPS(303),
    /** The edges provided to S2Builder cannot be assembled into a polyline. */
    BUILDER_EDGES_DO_NOT_FORM_POLYLINE(304),
    /**
     * There was an attempt to assemble a polygon from degenerate geometry without having specified
     * a predicate to decide whether the output is the empty polygon (containing no points) or the
     * full polygon (containing all points).
     */
    BUILDER_IS_FULL_PREDICATE_NOT_SPECIFIED(305),

    /** Codes in the range USER_DEFINED_START .. USER_DEFINED_END can be defined by clients. */
    USER_DEFINED_START(1000000),
    USER_DEFINED_END(9999999);

    private final int code;

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

  /** Prepares an S2Error instance for reuse by resetting it to its original state. */
  public void clear() {
    code = Code.NO_ERROR;
    text = "";
  }

  /**
   * Sets the error code and text description; the description is formatted according to the rules
   * defined in {@link Strings#lenientFormat(String, Object...)}, except that '%d' positional
   * arguments are also handled.
   *
   * <p>Note that this does not use String.format() because J2CL does not support that, and we
   * want to produce the same (or at least similar, readable) error texts on J2CL as Java. To
   * including Double values in the error text, %s in the format string and Double.toString() on the
   * value works well, except that the J2CL implementation of Double.toString() does not behave
   * exactly the same as the Java implementation, so messages will still be slightly different in
   * some cases (fewer decimal digits in J2CL, typically).
   *
   * <p>This method may be called more than once, so that various layers may surround the error
   * message with additional context:
   *
   * {@snippet :
   *     error.init(error.code(), "Loop %s: %s", j, error.text());
   * }
   */
  public void init(Code code, String format, Object... args) {
    this.code = code;
    format = format.replace("%d", "%s");
    this.text = Strings.lenientFormat(format, args);
  }

  /** Returns the code of this error. */
  public Code code() {
    return code;
  }

  /** Returns true if this error's code is NO_ERROR. */
  public boolean ok() {
    return code == Code.NO_ERROR;
  }

  /** Returns the text string. */
  public String text() {
    return text;
  }

  @Override
  public String toString() {
    if (code == Code.NO_ERROR) {
      return "OK";
    }
    return Strings.lenientFormat("%s: %s", code, text);
  }
}
