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

import static com.google.common.geometry.S2TextFormat.makePoint;
import static java.lang.Math.abs;
import static java.lang.Math.min;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;
import static org.junit.Assert.fail;

import com.google.common.annotations.GwtIncompatible;
import com.google.common.base.Preconditions;
import com.google.common.base.Splitter;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.primitives.IntPairVector;
import com.google.common.geometry.primitives.IntVector;
import com.google.common.primitives.UnsignedLong;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.function.Supplier;
import org.jspecify.annotations.Nullable;
import org.junit.Before;

/** Common code for geometry tests. */
public class GeometryTestCase {
  /**
   * How many ULP's (Units in the Last Place) we want to tolerate when comparing two numbers. The
   * gtest framework for C++ also uses 4, and documents why in gtest-internal.h.
   */
  public static final int MAX_ULPS = 4;

  /**
   * A not very accurate value for the radius of the Earth. For testing only, matches the value
   * in the C++ implementation S2Testing::kEarthRadiusKm.
   */
  public static final double APPROXIMATE_EARTH_RADIUS_KM = 6371.01;

  /**
   * The TestDataGenerator contains the Random used in unit tests as well as utility methods for
   * producing test data.
   */
  protected TestDataGenerator data;

  /** For convenience, provides direct access to the TestDataGenerator's Random. */
  protected Random rand() {
    return data.rand;
  }

  /**
   * Initializes the TestDataGenerator, and in particular, the random number generator it contains.
   */
  @Before
  public final void setUp() {
    data = new TestDataGenerator();
  }

  /**
   * Returns the result of calling 'supplier' with many (not all) S2 internal Java assertions
   * disabled. Normally, when Java assertions are enabled, as they typically are when running unit
   * tests, creating S2 objects with invalid input parameters will cause AssertionErrors to be
   * thrown. That is usually the desired behavior for unit testing.
   *
   * <p>However, it is also necessary to test code with invalid S2 objects, such as when testing
   * methods like {@link S2Polygon#findValidationError(S2Error)}. For example:
   *
   * {@snippet :
   * @Test
   * public testRectangleValidation() {
   *   // Test that creating an invalid S2LatLngRect with lo.lat() > hi.lat() will throw an
   *   // AssertionError.
   *   assertThrows(AssertionError.class, () ->
   *       new S2LatLngRect(S2LatLng.fromE7(30, 30), S2LatLng.fromE7(10, 10)));
   *
   *   // But unit tests may use uncheckedCreate to deliberately create an invalid S2LatLngRect.
   *   S2LatLngRect invalidRect = uncheckedCreate(() ->
   *       new S2LatLngRect(S2LatLng.fromE7(30, 30), S2LatLng.fromE7(10, 10)));
   *
   *   // Verify that isValid() correctly returns false for the invalid S2LatLngRect.
   *   assertFalse(invalidRect.isValid());
   * }
   * }
   *
   * If you need to create an S2 object with a method that may throw a checked exception, such as
   * {@code S2TaggedShapeCoder.FAST.decode()}, use {@link #unsafeCreate(Callable)} instead.
   */
  public static synchronized <T> T uncheckedCreate(Supplier<T> supplier) {
    boolean skip = S2.skipAssertions;
    try {
      S2.skipAssertions = true;
      return supplier.get();
    } finally {
      S2.skipAssertions = skip;
    }
  }

  /**
   * Returns the result of calling 'callable' with many (not all) S2 internal Java assertions
   * enabled, to allow the creation of invalid S2 objects. This is exactly like
   * {@link #uncheckedCreate(Supplier)} above, except that it can be used with a {@link Callable}
   * that may throw checked exceptions. This method will only throw a checked exception if the given
   * Callable does.
   *
   * <p>If you have a constructor or other method that cannot throw a checked exception, it is
   * recommended and more convenient to use {@link #uncheckedCreate(Supplier)} instead.
   */
  public static synchronized <T> T unsafeCreate(Callable<T> callable) throws Exception {
    boolean skip = S2.skipAssertions;
    try {
      S2.skipAssertions = true;
      return callable.call();
    } finally {
      S2.skipAssertions = skip;
    }
  }

  /**
   * Runs the given 'initializer' with many (not all) S2 internal Java assertions disabled. Like the
   * {@link #unsafeCreate(Callable)} method above, and used for the same purposes, but with
   * initializers like {@code S2Polygon#init(List<S2Loop>)} that return void. For example, in a unit
   * test:
   *
   * {@snippet :
   * @Test
   * public testPolygonIsValid() {
   *   List<S2Loop> crossingLoops = ImmutableList.of(...);
   *   S2Polygon p = new S2Polygon();
   *
   *   // Accidentally initializing an S2Polygon with crossing loops throws an AssertionError.
   *   assertThrows(AssertionError.class, () -> p.initNested(crossingLoops));
   *
   *   // But unsafeInitialize can deliberately initialize an invalid S2Polygon:
   *   unsafeInitialize(() -> p.initNested(crossingLoops));
   *
   *   // Now we can test isValid() on the invalid S2Polygon.
   *   assertFalse(p.isValid());
   * }
   * }
   */
  public static synchronized void uncheckedInitialize(Runnable initializer) {
    boolean skip = S2.skipAssertions;
    try {
      S2.skipAssertions = true;
      initializer.run();
    } finally {
      S2.skipAssertions = skip;
    }
  }

  /**
   * Returns the result of calling {@link S2TextFormat#makeLoop(String)} with assertions disabled.
   */
  public static S2Loop makeInvalidLoop(String str) {
    return uncheckedCreate(() -> S2TextFormat.makeLoop(str));
  }

  /**
   * Like {@link S2TextFormat#makePolygonOrDie(String)}, except that it does not normalize loops
   * (i.e., it gives you exactly what you asked for). Disables internal assertion of polygon
   * validity, allowing invalid polygons to be constructed.
   *
   * @throws IllegalArgumentException on unparsable input.
   */
  public static S2Polygon makeVerbatimPolygonOrDie(String str) {
    S2Polygon polygon = makeVerbatimPolygon(str);
    Preconditions.checkArgument(polygon != null, ": str == \"%s\"", str);
    return polygon;
  }

  /**
   * As {@link #makeVerbatimPolygonOrDie(String)} above, but does not throw IllegalArgumentException
   * on invalid input. Disables internal assertions of polygon validity, allowing invalid polygons
   * to be constructed. Returns null if conversion is unsuccessful.
   */
  public static @Nullable S2Polygon makeVerbatimPolygon(String str) {
    return uncheckedCreate(() -> S2TextFormat.internalMakePolygon(str, false, -1));
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
   * point comparison macros. So, our Java unit tests can use assertAlmostEquals(a, b) and have the
   * same behavior as our equivalent C++ unit tests using EXPECT_DOUBLE_EQ(a, b).
   *
   * <p>There is one exception for comparing by ULPs: If a and b are both non-zero, and have
   * different signs, they are never considered almost equals, as a ULP based comparison of floating
   * point numbers with different signs doesn't make sense. For details, read:
   * <a>https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/</a>
   */
  public static void assertAlmostEquals(String message, double a, double b) {
    assertDoubleUlpsWithin(message, a, b, MAX_ULPS);
  }

  /** Succeeds if and only if 'x' is within 4 units-in-the-last-place of 'y'. */
  public static void assertAlmostEquals(double a, double b) {
    assertAlmostEquals("", a, b);
  }

  /** Succeeds if and only if 'x' is less than or equal to 'y'. */
  public static <T extends Comparable<T>> void assertLessOrEqual(T x, T y) {
    assertTrue("Expected " + x + " <= " + y + " but it is not.", x.compareTo(y) <= 0);
  }

  /** Succeeds if and only if 'x' is greater than or equal to 'y'. */
  public static <T extends Comparable<T>> void assertGreaterOrEqual(T x, T y) {
    assertTrue("Expected " + x + " >= " + y + " but it is not.", x.compareTo(y) >= 0);
  }

  /** Succeeds if and only if 'x' is less than 'y'. */
  public static <T extends Comparable<T>> void assertLessThan(T x, T y) {
    assertTrue("Expected " + x + " < " + y + " but it is not.", x.compareTo(y) < 0);
  }

  /** Succeeds if and only if 'x' is greater than 'y'. */
  public static <T extends Comparable<T>> void assertGreaterThan(T x, T y) {
    assertTrue("Expected " + x + " > " + y + " but it is not.", x.compareTo(y) > 0);
  }

  /** Succeeds if and only if 'x' is between 'hi' and 'lo' inclusive. */
  public static <T extends Comparable<T>> void assertBetween(T x, T lo, T hi) {
    assertTrue("Expected " + x + " >= " + lo + " but it is not.", x.compareTo(lo) >= 0);
    assertTrue("Expected " + x + " <= " + hi + " but it is not.", x.compareTo(hi) <= 0);
  }

  /**
   * Assert that {@code val1} and {@code val2} are within the given {@code absError} of each other.
   *
   * <p>This implementation matches that of DoubleNearPredFormat() in gtest.cc, so our Java unit
   * tests can use assertDoubleNear() and get the same behavior as C++ tests using EXPECT_NEAR().
   */
  public static void assertDoubleNear(double val1, double val2, double absError) {
    assertDoubleNear("", val1, val2, absError);
  }

  /**
   * As above but with a custom error message that is prefixed to the report of the difference, if
   * the distance is larger than absError.
   */
  public static void assertDoubleNear(String message, double val1, double val2, double absError) {
    double diff = Math.abs(val1 - val2);
    if (diff <= absError) {
      return;
    }

    // Find the value which is closest to zero.
    double minAbs = min(Math.abs(val1), Math.abs(val2));
    // Find the distance to the next double from that value.
    double epsilon = nextUp(minAbs) - minAbs;

    if (!message.isEmpty()) {
      message = message + "\n";
    }
    // Detect the case where absError is so small that "near" is effectively the same as "equal",
    // and give an informative error message so that the situation can be more easily understood.
    // Don't do an epsilon check if absError is zero because that implies that an equality check
    // was actually intended.
    if (!Double.isNaN(val1) && !Double.isNaN(val2) && absError > 0 && absError < epsilon) {
      fail(
          Platform.formatString(
                  "%sThe difference between val1 (%s) and val2 (%s) is %s.\n"
                  + "The absError parameter (%s) is smaller than the minimum distance between"
                  + " doubles for numbers of this magnitude, which is %s, thus making this 'Near'"
                  + " check equivalent to an exact equality check.",
              message, val1, val2, diff, absError, epsilon));
    }
    fail(Platform.formatString(
        "%sThe difference between %s and %s is %s, which exceeds %s by %s.",
             message, val1, val2, diff, absError, (diff - absError)));
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

  /** Checks that two doubles are exactly equal as above, but with a custom error message. */
  public static void assertExactly(String message, double expected, double actual) {
    assertEquals(message, expected, actual, 0.0);
  }

  /**
   * Assert that the S1ChordAngle {@code actual} is exactly equal to {@code expected}.
   */
  public static void assertExactly(S1ChordAngle expected, S1ChordAngle actual) {
    assertEquals(expected.getLength2(), actual.getLength2(), 0.0);
  }

  /**
   * Assert that the S1ChordAngle distance {@code actual} is within {@code tolerance} of {@code
   * expected}.
   */
  public static void assertDistanceWithin(
      S1ChordAngle expected, S1ChordAngle actual, S1ChordAngle tolerance) {
    double diff = Math.abs(expected.getLength2() - actual.getLength2());
    assertTrue(
        "expected: "
            + expected
            + ", actual: "
            + actual
            + ", difference = "
            + S1ChordAngle.fromLength2(diff)
            + " which exceeds tolerance "
            + tolerance,
        diff < tolerance.getLength2());
  }

  /** Checks that two doubles are exactly equal and have the same sign. */
  public static void assertIdentical(double a, double b) {
    TestPlatform.assertIdentical(a, b);
  }

  /**
   * Checks that the 3D distance between {@code expected} and {@code actual} is at most {@code eps}
   * units.
   */
  public static void assertPointsWithinDistance(S2Point expected, S2Point actual, double eps) {
    assertTrue(
        "expected: " + expected + " but was: " + actual, expected.getDistance2(actual) < eps * eps);
  }

  /**
   * Checks that the spherical distance between {@code expected} and {@code actual} is at most
   * {@code degrees} degrees.
   */
  public static void assertPointsWithinDegrees(S2Point expected, S2Point actual, double degrees) {
    double actualDegrees = new S1Angle(actual, expected).degrees();
    assertTrue(
        "Expected "
            + actual.toDegreesString()
            + " to be within "
            + degrees
            + " degrees of "
            + expected.toDegreesString()
            + " but the actual distance is "
            + actualDegrees
            + " degrees.",
        actualDegrees <= degrees);
  }

  /**
   * Checks that each of the three components of the two given points are within the given {@code
   * eps} units.
   */
  public static void assertPointsNear(S2Point expected, S2Point actual, double eps) {
    assertDoubleNear(expected.getX(), actual.getX(), eps);
    assertDoubleNear(expected.getY(), actual.getY(), eps);
    assertDoubleNear(expected.getZ(), actual.getZ(), eps);
  }

  /**
   * Checks that the distance between {@code expected} and {@code actual} is at most
   * {@code maxErrorRadians}.
   */
  public static void assertApproxEquals(S2Point expected, S2Point actual, double maxErrorRadians) {
    assertTrue(
        "expected: " + expected + " but was: " + actual,
        S2.approxEquals(expected, actual, maxErrorRadians));
  }

  /** Assets that the actual shape {@link S2ShapeUtil#equals equals} the expected shape. */
  public static void assertShapesEqual(S2Shape expected, S2Shape actual) {
    assertTrue(S2ShapeUtil.equals(expected, actual));
  }

  /**
   * Assets that the methods of the two shape indexes return the same results, including all the
   * shapes in the indexes.
   */
  public static void assertShapeIndexesEqual(S2ShapeIndex expected, S2ShapeIndex actual) {
    assertTrue(S2ShapeUtil.equals(expected, actual));
  }

  /** Returns true if the two S2Errors have the same content. */
  public static boolean errorsEqual(S2Error a, S2Error b) {
    return a.toString().equals(b.toString());
  }

  /**
   * Returns true if the "expected" and "actual" IntVectors contain the same ints, in any order,
   * ignoring duplicates. Note that this implementation is not very efficient, and should only be
   * used in testing with small lists.
   */
  public void assertSameElements(IntVector expected, IntVector actual) {
    expected.forEach(
        i -> assertTrue("Expected contains " + i + " but actual does not.", actual.contains(i)));
    actual.forEach(
        i -> assertTrue("Actual contains " + i + " but expected does not.", expected.contains(i)));
  }

  /**
   * Returns true if the "expected" and "actual" IntPairVectors contain the same pairs, in any
   * order, ignoring duplicates. Note that the given IntPairVectors are sorted.
   */
  public void assertSameElements(IntPairVector expected, IntPairVector actual) {
    expected.sort();
    actual.sort();

    expected.forEach(
        (l, r) ->
            assertTrue(
                "Expected contains (" + l + ", " + r + ") but actual does not.",
                actual.contains(l, r)));
    actual.forEach(
        (l, r) ->
            assertTrue(
                "Actual contains (" + l + ", " + r + ") but expected does not.",
                expected.contains(l, r)));
  }

  /**
   * Returns an S1Angle approximately equal to the given distance in meters. Use S2Earth where
   * accuracy is important.
   */
  public static S1Angle metersToAngle(double meters) {
    return kmToAngle(0.001 * meters);
  }

  /**
   * Returns an S1Angle approximately equal to the given distance in kilometers. Use S2Earth where
   * accuracy is important.
   */
  public static S1Angle kmToAngle(double km) {
    return S1Angle.radians(km / APPROXIMATE_EARTH_RADIUS_KM);
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

  /** A concise expression to get an S2Point from a lat,lng in degrees. */
  protected S2Point ll(double lat, double lng) {
    return S2LatLng.fromDegrees(lat, lng).toPoint();
  }

  /** Rotates the given list of edges left by 'n' places. */
  public static void rotate(List<S2Edge> edges, int n) {
    Collections.rotate(edges, -1 * n);
  }

  /** Randomly permutes the given list of edges. */
  public static void shuffle(List<S2Edge> edges) {
    Collections.shuffle(edges);
  }

  /**
   * As {@link #checkCovering(S2Region, S2CellUnion, boolean, S2CellId)}, but creates a default and
   * invalid S2CellId for the last argument.
   */
  public void checkCovering(S2Region region, S2CellUnion covering, boolean checkTight) {
    checkCovering(region, covering, checkTight, new S2CellId());
  }

  /**
   * Checks that "covering" completely covers the given region. If "checkTight" is true, also checks
   * that it does not contain any cells that do not intersect the given region. ("id" is only used
   * internally.)
   */
  public void checkCovering(S2Region region, S2CellUnion covering, boolean checkTight, S2CellId id) {
    if (!id.isValid()) {
      for (int face = 0; face < 6; ++face) {
        checkCovering(region, covering, checkTight, S2CellId.fromFacePosLevel(face, 0, 0));
      }
      return;
    }

    if (!region.mayIntersect(new S2Cell(id))) {
      // If region does not intersect id, then neither should the covering.
      if (checkTight) {
        assertFalse(covering.intersects(id));
      }

    } else if (!covering.contains(id)) {
      // The region may intersect id, but we can't assert that the covering
      // intersects id because we may discover that the region does not actually
      // intersect upon further subdivision. (MayIntersect is not exact.)
      assertFalse(region.contains(new S2Cell(id)));
      assertFalse(id.isLeaf());
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

  /** Expects obj.equals(deserialize(serialize(obj))) to be true. */
  @GwtIncompatible("Javascript doesn't support Java serialization.")
  protected void doSerializationTest(Object obj) {
    doSerializationTest(obj, null);
  }

  /**
   * If a non-null comparator is provided, expects obj.compare(deserialize(serialize(obj))) to be
   * zero. Otherwise, expects obj.equals(deserialize(serialize(obj))) to be true.
   */
  @GwtIncompatible("Javascript doesn't support Java serialization.")
  protected <T> void doSerializationTest(T obj, Comparator<T> comparator) {
    try {
      TestPlatform.testSerialization(obj, comparator);
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  /**
   * Returns a set of strings that are the "toToken" representation of the cell ids in the given
   * S2CellUnion.
   */
  public static Set<String> toTokenSet(S2CellUnion cells) {
    HashSet<String> tokens = new HashSet<>();
    cells.forEach(cell -> tokens.add(cell.toToken()));
    return tokens;
  }

  /**
   * Returns a set of strings that are the "toToken" representation of the given set of cell ids.
   */
  public static Set<String> toTokenSet(Set<S2CellId> cells) {
    HashSet<String> tokens = new HashSet<>();
    cells.forEach(cell -> tokens.add(cell.toToken()));
    return tokens;
  }

  /**
   * Returns a set of Strings that are the "toToken" representations of the cell id of every node in
   * the given density tree.
   */
  public static Set<String> toTokenSet(S2DensityTree tree) {
    HashSet<String> tokens = new HashSet<>();
    tree.visitCells(
        (S2CellId cellId, S2DensityTree.Cell node) -> {
          tokens.add(cellId.toToken());
          return S2DensityTree.CellVisitor.Action.ENTER_CELL;
        });
    return tokens;
  }

  /** Returns a String representation of the given tree. */
  public static String toString(S2DensityTree tree) {
    return decodedTreeToString(tree.decode());
  }

  /**
   * Returns a multi-line string representation of a decoded S2DensityTree, using indenting to show
   * the tree structure. Useful for debugging and tests.
   */
  public static String decodedTreeToString(Map<S2CellId, Long> decodedTree) {
    StringBuilder sb = new StringBuilder();
    dumpTree(Arrays.asList(S2CellId.FACE_CELLS), decodedTree, 0, sb);
    return sb.toString();
  }

  /** Recursive helper for decodedTreeToString. Appends to the given StringBuilder. */
  private static void dumpTree(
      Iterable<S2CellId> start, Map<S2CellId, Long> decoded, int indent, StringBuilder sb) {
    for (S2CellId id : start) {
      if (decoded.containsKey(id)) {
        for (int i = 0; i < indent; i++) {
          sb.append("    ");
        }
        sb.append(Platform.formatString("Level %s: node[%s] %s with weight %s\n",
                      id.level(), id.toToken(), id, decoded.get(id)));
        dumpTree(id.children(), decoded, indent + 1, sb);
      }
    }
  }

  /**
   * Like S2Shape.MutableEdge, but immutable. Overrides hashCode() and equals() so it can be used
   * as a HashMap key. This is not very memory efficient, but is convenient in tests. Equality to
   * a MutableEdge is also defined and works as expected.
   */
  static class ImmutableEdge {
    final S2Point a;
    final S2Point b;

    public ImmutableEdge(S2Point a, S2Point b) {
      this.a = a;
      this.b = b;
    }
    /** Constructs an ImmutableEdge with the same endpoints as the given S2Shape.MutableEdge. */
    public ImmutableEdge(MutableEdge e) {
      this.a = e.a;
      this.b = e.b;
    }

    /** Returns true iff 'point' is either endpoint of this edge. */
    public boolean isEndpoint(S2Point point) {
      return a.equalsPoint(point) || b.equalsPoint(point);
    }

    /**
     * Returns true if this ImmutableEdge has the same endpoints as the given S2Shape.MutableEdge
     * currently does.
     */
    public boolean isEqualTo(MutableEdge other) {
      return a.equalsPoint(other.a) && b.equalsPoint(other.b);
    }

    /**
     * Returns true if this ImmutableEdge has the same endpoints as the other ImmutableEdge.
     */
    public boolean isEqualTo(ImmutableEdge other) {
      return a.equalsPoint(other.a) && b.equalsPoint(other.b);
    }

    @Override
    public boolean equals(Object other) {
      if (other instanceof ImmutableEdge) {
        return isEqualTo((ImmutableEdge) other);
      }
      if (other instanceof MutableEdge) {
        return isEqualTo((MutableEdge) other);
      }
      return false;
    }

    @Override
    public int hashCode() {
      return a.hashCode() * 3 + b.hashCode();
    }

    public String toDegreesString() {
      return a.toDegreesString() + "-" + b.toDegreesString();
    }

    @Override
    public String toString() {
      return toDegreesString();
    }
  }

  /** A simple implementation of Pair for unit tests. */
  public static final class Pair<A, B> { // TODO(user): Use a record when J2CL supports it.
    public A first;
    public B second;

    /** Constructs a new Pair with the given values. */
    public Pair(A a, B b) {
      this.first = a;
      this.second = b;
    }

    /** Returns the first element of this Pair. */
    public A getFirst() {
      return first;
    }

    /** Returns the second element of this Pair. */
    public B getSecond() {
      return second;
    }

    /** Returns a new Pair containing the given values. */
    public static <A, B> Pair<A, B> of(A a, B b) {
      return new Pair<A, B>(a, b);
    }

    @Override
    public boolean equals(Object other) {
      if (other instanceof Pair) {
        Pair<?, ?> otherPair = (Pair<?, ?>) other;
        return first.equals(otherPair.first) && second.equals(otherPair.second);
      }
      return false;
    }

    @Override
    public int hashCode() {
      return first.hashCode() * 3 + second.hashCode();
    }
  }
}
