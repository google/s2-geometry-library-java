/*
 * Copyright 2013 Google Inc.
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

import com.google.errorprone.annotations.FormatMethod;
import java.io.PrintStream;
import java.math.BigDecimal;
import java.util.Locale;
import java.util.logging.Logger;

/**
 * Contains utility methods which require different implementations for the Java and J2CL-generated
 * Javascript implementations. This contains the Java implementations. The J2CL versions are in
 * java/com/google/gwt/corp/emul/com/google/common/geometry/Platform.java.
 */
public final class Platform {

  private Platform() {}

  /**
   * @see Math#IEEEremainder(double, double)
   */
  public static double IEEEremainder(double f1, double f2) {
    return Math.IEEEremainder(f1, f2);
  }

  /**
   * @see Math#getExponent(double)
   */
  public static int getExponent(double d) {
    return Math.getExponent(d);
  }

  /**
   * Returns the {@link Logger} for the class.
   *
   * @see Logger#getLogger(String)
   */
  static Logger getLoggerForClass(Class<?> clazz) {
    return Logger.getLogger(clazz.getCanonicalName());
  }

  /**
   * Workaround for the difficulty of getting identical formatting of numbers between J2CL and Java.
   * Supports different "expected" strings for Java and J2CL. Those should normally be similar, but
   * J2CL behaves differently in ways that are hard to work around otherwise. For example, J2CL's
   * Double.toString() will use fewer decimal digits, in some cases.
   *
   * <p>TODO(user): Remove this method in favor of a consistent cross-platform method of
   * formatting doubles.
   */
  public static void assertStringEquals(String expectedJava, String expectedJ2cl, String actual) {
    if (!expectedJava.contentEquals(actual)) {
      throw new AssertionError("Expected '" + expectedJava + "' to equal '" + actual + "'");
    }
  }

  /**
   * Invokes {@code stream.printf} with the arguments. The J2CL client just prints the format string
   * and the arguments separately. Using this method is not recommended; you should instead
   * construct strings with normal string concatenation whenever possible, so it will work the same
   * way in normal Java and J2CL client versions.
   */
  public static void printf(PrintStream stream, String format, Object... params) {
    stream.printf(format, params);
  }

  /**
   * Returns {@code String.format} with the arguments. The J2CL client just returns a string
   * consisting of the format string with the parameters concatenated to the end of it. Using this
   * method is not recommended; you should instead construct strings with normal string
   * concatenation whenever possible, so it will work the same way in normal Java and J2CL client
   * versions.
   */
  @FormatMethod
  public static String formatString(String format, Object... params) {
    return String.format(format, params);
  }

  /**
   * A cross platform (Java / Javascript) method to format a double as a string, very similar to
   * the behavior of printf("%.15g",d) in C++. Provides (almost?) always consistent behavior between
   * Java and J2CL. This method should only be used where cross-platform formatting is required, as
   * it is slower than Double.toString(), and nonstandard.
   */
  public static String formatDouble(double d) {
    StringBuilder out = new StringBuilder();
    if (d == 0d) {
      return "0";
    }
    // Style 'g' uses either 'e' or 'f', depending on the magnitude of the number.
    out.append(String.format(Locale.US, "%.15g", d));

    // If formatted with style 'e', the 'e' is always in the same place relative to the length,
    // and the string will be at least five chars long, like "1e-20".
    if ((out.length() >= 5) && (out.charAt(out.length() - 4) == 'e')) {
      // Remove trailing zeros before the 'e'.
      while ((out.length() >= 5) && out.charAt(out.length() - 5) == '0') {
        out.deleteCharAt(out.length() - 5);
      }
      // Remove trailing decimal point.
      if (out.charAt(out.length() - 5) == '.') {
        out.deleteCharAt(out.length() - 5);
      }
    } else {
      // Otherwise, it was formatted with style 'f'. Remove trailing zeros.
      while (out.length() > 0 && out.charAt(out.length() - 1) == '0') {
        out.setLength(out.length() - 1);
      }
      // Remove trailing decimal point.
      if (out.length() > 0 && out.charAt(out.length() - 1) == '.') {
        out.setLength(out.length() - 1);
      }
    }
    return out.toString();
  }

  /** A portable way to hash a double value. */
  public static long doubleHash(double value) {
    return Double.doubleToLongBits(value);
  }

  /**
   * Returns the sign of the determinant of the matrix constructed from the three column vectors
   * {@code a}, {@code b}, and {@code c}. This operation is very robust for small determinants, but
   * is extremely slow and should only be used if performance is not a concern or all faster
   * techniques have been exhausted.
   */
  public static int sign(S2Point a, S2Point b, S2Point c) {
    try {
      Real bycz = Real.strictMul(b.y, c.z);
      Real bzcy = Real.strictMul(b.z, c.y);
      Real bzcx = Real.strictMul(b.z, c.x);

      Real bxcz = Real.strictMul(b.x, c.z);
      Real bxcy = Real.strictMul(b.x, c.y);
      Real bycx = Real.strictMul(b.y, c.x);

      Real bcx = bycz.sub(bzcy);
      Real bcy = bzcx.sub(bxcz);
      Real bcz = bxcy.sub(bycx);
      Real x = bcx.strictMul(a.x);
      Real y = bcy.strictMul(a.y);
      Real z = bcz.strictMul(a.z);
      return x.add(y).add(z).signum();
    } catch (ArithmeticException unused) {
      return 0;
    }
  }

  /**
   * Returns the size of an ulp of the argument. An ulp of a double value is the positive distance
   * between this floating-point value and the double next larger in magnitude.
   */
  public static double ulp(double x) {
    return Math.ulp(x);
  }

  /**
   * Returns the next representable value in the direction of 'dir' starting from 'x', emulating the
   * behavior of {@link Math#nextAfter}.
   */
  public static double nextAfter(double x, double dir) {
    return Math.nextAfter(x, dir);
  }

  /**
   * Returns a new {@code BigDecimal} instance whose value is the exact decimal representation of
   * {@code x}, emulating the behavior of {@link BigDecimal#BigDecimal(double)}.
   */
  static BigDecimal newBigDecimal(double x) {
    return new BigDecimal(x);
  }
}
