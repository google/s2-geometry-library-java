/*
 * Copyright 2005 Google Inc.
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

import static com.google.common.geometry.TestDataGenerator.makePoint;
import static java.lang.Math.abs;
import static java.lang.Math.min;
import static junit.framework.TestCase.assertEquals;

import com.google.common.base.Splitter;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.primitives.UnsignedLong;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import junit.framework.TestCase;

/** Common code for geometry tests. */
public strictfp class GeometryTestCase extends TestCase {
  /**
   * How many ULP's (Units in the Last Place) we want to tolerate when comparing two numbers. The
   * gtest framework for C++ also uses 4, and documents why in gtest-internal.h.
   */
  public static final int MAX_ULPS = 4;

  /**
   * The TestDataGenerator contains the Random used in unit tests as well as utility methods for
   * producing test data.
   */
  protected TestDataGenerator data;

  /**
   * Initializes the TestDataGenerator, and in particular, the random number generator it contains.
   */
  @Override
  protected void setUp() {
    data = new TestDataGenerator();
  }

  /**
   * Returns the next (more positive) representable number after x.
   */
  public static double nextUp(double x) {
    // java.lang.Math provides nextUp and nextDown, but not for Javascript in JCL, and they're not
    // in Platform either.
    return Platform.nextAfter(x, Double.POSITIVE_INFINITY);
  }

  /*
   * Returns the adjacent, more negative, representable number before x.
   */
  public static double nextDown(double x) {
    return Platform.nextAfter(x, Double.NEGATIVE_INFINITY);
  }

  /**
   * Tests that two double values have the same sign and are within 'maxUlps' of each other.
   */
  public static void assertDoubleUlpsWithin(String message, double a, double b, int maxUlps) {
    // The IEEE standard says that any comparison operation involving a NAN must return false.
    if (Double.isNaN(a)) {
      fail("'a' is NaN. " + message);
    }
    if (Double.isNaN(b)) {
      fail("'b' is NaN. " + message);
    }

    // Handle the exact equality case fast, as well as special cases like +0 == -0, and infinity.
    if (a == b) {
      return;
    }

    // If the signs are different, don't compare by ULP.
    if (Math.copySign(1.0, a) != Math.copySign(1.0, b)) {
      fail(a + " and " + b + " are not equal and have different signs. " + message);
    }

    UnsignedLong uA = UnsignedLong.fromLongBits(Double.doubleToLongBits(a));
    UnsignedLong uB = UnsignedLong.fromLongBits(Double.doubleToLongBits(b));
    int ulpsDiff = uA.minus(uB).intValue();
    assertTrue(a + " and " + b + " differ by " + ulpsDiff + " units in the last place, expected <= "
            + maxUlps + ". " + message,
        abs(ulpsDiff) <= maxUlps);
  }

  /**
   * Tests that two double values are almost equal. Uses ULP-based comparison to automatically pick
   * an error bound that is appropriate for the operands. "Almost equal" here means "a is at most
   * MAX_ULPS (which is 4) ULP's away from b".
   *
   * <p>This is similar to the implementation of "AlmostEquals" which underlies the gtest floating
   * point comparison macros. So, our Java unit tests can use assertDoubleEquals(a, b) and have the
   * same behavior as our equivalent C++ unit tests using EXPECT_DOUBLE_EQ(a, b).
   *
   * <p>There is one exception for comparing by ULPs: If a and b are both non-zero, and have
   * different signs, they are never considered almost equals, as a ULP based comparison of floating
   * point numbers with different signs doesn't make sense. For details, read:
   * <a>https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/</a>
   */
  public static void assertDoubleEquals(String message, double a, double b) {
    assertDoubleUlpsWithin(message, a, b, MAX_ULPS);
  }

  /** Returns true if and only if 'x' is within 4 units-in-the-last-place of 'y'. */
  public static void assertDoubleEquals(double a, double b) {
    assertDoubleEquals("", a, b);
  }

  /** Returns true if and only if 'x' is less than or equal to 'y'. */
  public static <T extends Comparable<T>> void assertLessOrEqual(T x, T y) {
    assertTrue("Expected " + x + " <= " + y + " but it is not.", x.compareTo(y) <= 0);
  }

  /** Returns true if and only if 'x' is greater than or equal to 'y'. */
  public static <T extends Comparable<T>> void assertGreaterOrEqual(T x, T y) {
    assertTrue("Expected " + x + " >= " + y + " but it is not.", x.compareTo(y) >= 0);
  }

  /** Returns true if and only if 'x' is less than 'y'. */
  public static <T extends Comparable<T>> void assertLessThan(T x, T y) {
    assertTrue("Expected " + x + " < " + y + " but it is not.", x.compareTo(y) < 0);
  }

  /** Returns true if and only if 'x' is greater than 'y'. */
  public static <T extends Comparable<T>> void assertGreaterThan(T x, T y) {
    assertTrue("Expected " + x + " > " + y + " but it is not.", x.compareTo(y) > 0);
  }


  /**
   * Assert that {@code val1} and {@code val2} are within the given {@code absError} of each other.
   *
   * <p>This implementation matches that of DoubleNearPredFormat() in gtest.cc, so our Java unit
   * tests can use assertDoubleNear() and get the same behavior as C++ tests using EXPECT_NEAR().
   */
  public static void assertDoubleNear(double val1, double val2, double absError) {
    double diff = Math.abs(val1 - val2);
    if (diff <= absError) {
      return;
    }

    // Find the value which is closest to zero.
    double minAbs = min(Math.abs(val1), Math.abs(val2));
    // Find the distance to the next double from that value.
    double epsilon = nextUp(minAbs) - minAbs;

    // Detect the case where absError is so small that "near" is effectively the same as "equal",
    // and give an informative error message so that the situation can be more easily understood.
    // Don't do an epsilon check if absError is zero because that implies that an equality check
    // was actually intended.
    if (!Double.isNaN(val1) && !Double.isNaN(val2) && absError > 0 && absError < epsilon) {
      fail("The difference between val1 (" + val1 + ") and val2 (" + val2 + ") is " + diff
            + ".\nThe absError parameter (" + absError + ") is smaller than the minimum "
            + "distance between doubles for numbers of this magnitude, which is " + epsilon
            + ", thus making this 'Near' check equivalent to an exact equality check.");
    }
    fail("The difference between " + val1 + " and " + val2 + " is " + diff + ", which exceeds "
          + absError + " by " + (diff - absError) + ".");
  }

  /**
   * Assert that {@code a} and {@code b} are within 1e-9 of each other.
   */
  public static void assertDoubleNear(double a, double b) {
    assertDoubleNear(a, b, 1e-9);
  }

  /**
   * Checks that two doubles are exactly equal, although note that 0.0 exactly equals -0.0: use
   * assertIdentical to differentiate positive and negative zero.
   *
   * <p>(JUnit 3 allows leaving off the third parameter of assertEquals, with a default of zero, but
   * JUnit4 does not. In S2 we often want to check that two doubles are exactly equal, so this is a
   * bit cleaner.)
   */
  public static void assertExactly(double expected, double actual) {
    assertEquals(expected, actual, 0.0);
  }

  /**
   * Assert that the S1ChordAngle {@code actual} is exactly equal to {@code expected}.
   */
  public static void assertExactly(S1ChordAngle expected, S1ChordAngle actual) {
    assertEquals(expected.getLength2(), actual.getLength2(), 0.0);
  }

  /** Checks that two doubles are exactly equal and have the same sign. */
  public static void assertIdentical(double a, double b) {
    assertEquals(Double.doubleToRawLongBits(a), Double.doubleToRawLongBits(b));
  }

  /**
   * Checks that the 3D distance between {@code expected} and {@code actual} is at most {@code eps}
   * units.
   */
  public static void assertEquals(S2Point expected, S2Point actual, double eps) {
    assertTrue(
        "expected: " + expected + " but was: " + actual, expected.getDistance2(actual) < eps * eps);
  }

  /** Assets that the actual shape {@link S2ShapeUtil#equals equals} the expected shape. */
  public static void assertShapesEqual(S2Shape expected, S2Shape actual) {
    assertTrue(S2ShapeUtil.equals(expected, actual));
  }

  /**
   * As {@link #checkCovering(S2Region, S2CellUnion, boolean, S2CellId)}, but creates a default and
   * invalid S2CellId for the last argument.
   */
  void checkCovering(S2Region region, S2CellUnion covering, boolean checkTight) {
    checkCovering(region, covering, checkTight, new S2CellId());
  }

  /** A helper for testing {@link S2Coder} implementations. */
  static <T> T roundtrip(S2Coder<T> coder, T value) {
    ByteArrayOutputStream output = new ByteArrayOutputStream();
    try {
      coder.encode(value, output);
      return coder.decode(Bytes.fromByteArray(output.toByteArray()));
    } catch (IOException e) {
      throw new AssertionError(e);
    }
  }

  /**
   * Checks that "covering" completely covers the given region. If "checkTight" is true, also checks
   * that it does not contain any cells that do not intersect the given region. ("id" is only used
   * internally.)
   */
  void checkCovering(S2Region region, S2CellUnion covering, boolean checkTight, S2CellId id) {
    if (!id.isValid()) {
      for (int face = 0; face < 6; ++face) {
        checkCovering(region, covering, checkTight, S2CellId.fromFacePosLevel(face, 0, 0));
      }
      return;
    }

    if (!region.mayIntersect(new S2Cell(id))) {
      // If region does not intersect id, then neither should the covering.
      if (checkTight) {
        assertTrue(!covering.intersects(id));
      }

    } else if (!covering.contains(id)) {
      // The region may intersect id, but we can't assert that the covering
      // intersects id because we may discover that the region does not actually
      // intersect upon further subdivision. (MayIntersect is not exact.)
      assertTrue(!region.contains(new S2Cell(id)));
      assertTrue(!id.isLeaf());
      S2CellId end = id.childEnd();
      for (S2CellId child = id.childBegin(); !child.equals(end); child = child.next()) {
        checkCovering(region, covering, checkTight, child);
      }
    }
  }

  /**
   * Asserts that the first N calls to {@link S2Shape#getEdge} returns edges are equivalent to those
   * specified in the edge list, described by the following format:
   *
   * <pre>edge1, edge2, ..., edgeN</pre>
   *
   * where edges are:
   *
   * <pre>point1 | point2</pre>
   *
   * and points are:
   *
   * <pre>lat:lng</pre>
   *
   * <p>Example:
   *
   * <table>
   * <tr><td>Two edges</td><td><pre>1:2 | 2:3, 4:5 | 6:7</pre></td></tr>
   * </table>
   */
  static void checkFirstNEdges(S2Shape shape, String edgeList) {
    MutableEdge result = new MutableEdge();
    consumeIndexedEdges(
        edgeList,
        (int edgeId, MutableEdge expected) -> {
          shape.getEdge(edgeId, result);
          assertEquals(expected.getStart(), result.getStart());
          assertEquals(expected.getEnd(), result.getEnd());
        });
  }

  /**
   * Similar to {@link #checkFirstNEdges}, except that {@link S2Shape#getChainEdge} is used with
   * specified chain id to return the edges.
   */
  static void checkFirstNChainEdges(S2Shape shape, int chainId, String edgeList) {
    MutableEdge result = new MutableEdge();
    consumeIndexedEdges(
        edgeList,
        (int offset, MutableEdge expected) -> {
          shape.getChainEdge(chainId, offset, result);
          assertEquals(expected.getStart(), result.getStart());
          assertEquals(expected.getEnd(), result.getEnd());
        });
  }

  private interface IndexedEdgeConsumer {
    void apply(int index, MutableEdge edge);
  }

  private static void consumeIndexedEdges(String edgeStrList, IndexedEdgeConsumer consumer) {
    MutableEdge edge = new MutableEdge();
    int i = 0;
    for (String edgeStr : Splitter.on(',').split(edgeStrList)) {
      int j = edgeStr.indexOf('|');
      edge.set(
          makePoint(edgeStr.substring(0, j).trim()), makePoint(edgeStr.substring(j + 1).trim()));
      consumer.apply(i++, edge);
    }
  }

  /** Asserts that the first N chain starts are equivalent to those specified in the list. */
  static void checkFirstNChainStarts(S2Shape shape, int... starts) {
    int chainId = 0;
    for (int start : starts) {
      assertEquals(start, shape.getChainStart(chainId++));
    }
  }

  /** Asserts that the first N chain lengths are equivalent to those specified in the list. */
  static void checkFirstNChainLengths(S2Shape shape, int... lengths) {
    int chainId = 0;
    for (int length : lengths) {
      assertEquals(length, shape.getChainLength(chainId++));
    }
  }
}
