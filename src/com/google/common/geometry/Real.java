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

import com.google.common.annotations.GwtIncompatible;

import java.math.BigDecimal;

/**
 * This class provides portable support for several exact arithmetic operations on double values,
 * without loss of precision. It stores an array of double values, and operations that require
 * additional bits of precision return Reals with larger arrays.
 *
 * <p>Converting a sequence of a dozen strictfp arithmetic operations to use Real can take up to 20
 * times longer than the natural but imprecise approach of using built in double operators like +
 * and *. Compared to other approaches like BigDecimal that consume more memory and typically slow
 * operations down by a factor of 100, that's great, but use of this class should still be avoided
 * when imprecise results will suffice.
 *
 * <p>This class exists as a package private element of the geometry library for several simple
 * algorithms such as {@link S2#robustCCW(S2Point, S2Point, S2Point)}, that require arbitrary
 * precision arithmetic. It could be made suitable for general usage by adding robust
 * implementations of multiplication and division between two Reals, and a toString() implementation
 * that prints the exact summation of all the components.
 *
 * <p>Many of the algorithms in this class were adapted from the multiple components technique for
 * extended 64-bit IEEE 754 floating point precision, as described in:
 *
 * <pre>
 * Robust Adaptive Floating-Point Geometric Predicates
 * Jonathan Richard Shewchuk
 * School of Computer Science
 * Carnegie Mellon University
 * </pre>
 *
 * <p>Faster adaptive techniques are also presented in that paper, but are not implemented here.
 */
@GwtIncompatible(value = "No javascript support for strictfp.")
strictfp class Real extends Number {
  private static final long serialVersionUID = 1L;

  /**
   * Used to split doubles into two half-length values, for exact multiplication.
   * The value should be Math.pow(2, Math.ceil(mantissaBits / 2)) + 1.
   */
  private static final double SPLITTER;

  static {
    // Find half ulp(1). We could use Math.ulp but it's not supported on GWT.
    double epsilon = 1.0;
    do {
      epsilon *= 0.5;
    } while (1.0 + epsilon != 1.0);
    int mantissaBits = (int) Math.round(-Math.log(epsilon) / Math.log(2));
    SPLITTER = (1 << ((mantissaBits + 1) / 2)) + 1;
  }

  /** Returns the result of a + b, without loss of precision. */
  public static Real add(double a, double b) {
    double x = a + b;
    double error = twoSumError(a, b, x);
    return new Real(error, x);
  }

  /** Returns the result of a - b, without loss of precision. */
  public static Real sub(double a, double b) {
    double x = a - b;
    double error = twoDiffError(a, b, x);
    return new Real(error, x);
  }

  /** Returns the result of a * b, without loss of precision. */
  public static Real mul(double a, double b) {
    double x = a * b;
    double bhi = splitHigh(b);
    double blo = splitLow(b, bhi);
    double error = twoProductError(a, bhi, blo, x);
    return new Real(error, x);
  }

  /**
   * A sequence of ordinary double values, ordered by magnitude in ascending
   * order, containing no zeroes and with no overlapping base 2 digits.
   */
  private final double[] values;

  /** Creates a Real based on the given double value. */
  public Real(double value) {
    values = new double[] {value};
  }

  private Real(double ... values) {
    this.values = values;
  }

  /** Returns the result of a + b, without loss of precision. */
  public Real add(Real that) {
    return add(this, that, false);
  }

  /** Returns the result of a - b, without loss of precision. */
  public Real sub(Real that) {
    return add(this, that, true);
  }

  /**
   * Returns the result of adding together the components of a and b, inverting
   * each element of b if negateB is true.
   */
  private static Real add(Real a, Real b, boolean negateB) {
    double bSign = negateB ? -1 : 1;
    double[] result = new double[a.values.length + b.values.length];
    int aIndex = 0;
    int bIndex = 0;

    double sum, newSum, error;
    if (smallerMagnitude(a.values[aIndex], b.values[bIndex])) {
      sum = a.values[aIndex++];
    } else {
      sum = bSign * b.values[bIndex++];
    }

    int resultIndex = 0;
    double smaller;
    if ((aIndex < a.values.length) && (bIndex < b.values.length)) {
      if (smallerMagnitude(a.values[aIndex], b.values[bIndex])) {
        smaller = a.values[aIndex++];
      } else {
        smaller = bSign * b.values[bIndex++];
      }
      newSum = smaller + sum;
      error = fastTwoSumError(smaller, sum, newSum);
      sum = newSum;
      if (error != 0.0) {
        result[resultIndex++] = error;
      }
      while ((aIndex < a.values.length) && (bIndex < b.values.length)) {
        if (smallerMagnitude(a.values[aIndex], b.values[bIndex])) {
          smaller = a.values[aIndex++];
        } else {
          smaller = bSign * b.values[bIndex++];
        }
        newSum = sum + smaller;
        error = twoSumError(sum, smaller, newSum);
        sum = newSum;
        if (error != 0.0) {
          result[resultIndex++] = error;
        }
      }
    }
    while (aIndex < a.values.length) {
      smaller = a.values[aIndex++];
      newSum = sum + smaller;
      error = twoSumError(sum, smaller, newSum);
      sum = newSum;
      if (error != 0.0) {
        result[resultIndex++] = error;
      }
    }
    while (bIndex < b.values.length) {
      smaller = bSign * b.values[bIndex++];
      newSum = sum + smaller;
      error = twoSumError(sum, smaller, newSum);
      sum = newSum;
      if (error != 0.0) {
        result[resultIndex++] = error;
      }
    }
    if ((sum != 0.0) || (resultIndex == 0)) {
      result[resultIndex++] = sum;
    }

    if (result.length > resultIndex) {
      result = copyOf(result, resultIndex);
    }
    return new Real(result);
  }

  /** Returns true if the magnitude of a is less than the magnitude of b. */
  private static boolean smallerMagnitude(double a, double b) {
    return (b > a) == (b > -a);
  }

  /** Returns the result of this * scale, without loss of precision. */
  public Real mul(double scale) {
    double[] result = new double[values.length * 2];
    double scaleHigh = splitHigh(scale);
    double scaleLow = splitLow(scale, scaleHigh);
    double quotient = values[0] * scale;
    double error = twoProductError(values[0], scaleHigh, scaleLow, quotient);
    int resultIndex = 0;
    if (error != 0) {
      result[resultIndex++] = error;
    }
    for (int i = 1; i < values.length; i++) {
      double term = values[i] * scale;
      double termError = twoProductError(values[i], scaleHigh, scaleLow, term);
      double sum = quotient + termError;
      error = twoSumError(quotient, termError, sum);
      if (error != 0) {
        result[resultIndex++] = error;
      }
      quotient  = term + sum;
      error = fastTwoSumError(term, sum, quotient);
      if (error != 0) {
        result[resultIndex++] = error;
      }
    }
    if ((quotient != 0.0) || (resultIndex == 0)) {
      result[resultIndex++] = quotient;
    }
    if (result.length > resultIndex) {
      result = copyOf(result, resultIndex);
    }
    return new Real(result);
  }

  /** Returns the negative of this number. */
  public Real negate() {
    double[] copy = new double[values.length];
    for (int i = values.length - 1; i >= 0; i--) {
      copy[i] = -values[i];
    }
    return new Real(copy);
  }

  /** Returns the signum of this number more quickly than via Math.signum(doubleValue()). */
  public int signum() {
    double msb = values[values.length - 1];
    if (msb > 0) {
      return 1;
    } else if (msb < 0) {
      return -1;
    } else {
      return 0;
    }
  }

  /** Returns the string representation of the double value nearest this Real. */
  @Override public String toString() {
    return Double.toString(doubleValue());
  }

  @Override public int intValue() {
    return (int) longValue();
  }

  @Override public long longValue() {
    return Math.round(doubleValue());
  }

  @Override public float floatValue() {
    return (float) doubleValue();
  }

  @Override
  public double doubleValue() {
    // Since the components are guaranteed to have no overlapping digits, we
    // could simply sum them without loss of precision... but to return a double
    // we truncate to the 53 bits of the largest exponent.
    double sum = 0;
    for (double value : values) {
      sum += value;
    }
    return sum;
  }

  /** Returns a BigDecimal representation of this extended precision real value. */
  public BigDecimal bigValue() {
    BigDecimal sum = new BigDecimal(values[0]);
    for (int i = 1; i < values.length; i++) {
      sum = sum.add(new BigDecimal(values[i]));
    }
    return sum.stripTrailingZeros();
  }

  private static double[] copyOf(double[] array, int newLength) {
    double[] result = new double[newLength];
    for (int i = 0; i < newLength; i++) {
      result[i] = array[i];
    }
    return result;
  }

  /** Returns the error in the sum x=a+b, when |a|>=|b|. */
  private static double fastTwoSumError(double a, double b, double x) {
    return b - (x - a);
  }

  /**
   * Returns the error in the sum x=a+b, when the relative magnitudes of a and
   * b are not known in advance.
   */
  private static double twoSumError(double a, double b, double x) {
    double error = x - a;
    return (a - (x - error)) + (b - error);
  }

  /** Returns the error in the difference x=a-b. */
  private static double twoDiffError(double a, double b, double x) {
    double error = a - x;
    return (a - (x + error)) + (error - b);
  }

  /** Returns the high split for the given value. */
  private static double splitHigh(double a) {
    double c = SPLITTER * a;
    return c - (c - a);
  }

  /**
   * Returns the low split for the given value and previously-computed high
   * split as returned by {@link #splitHigh(double)}.
   */
  private static double splitLow(double a, double ahi) {
    return a - ahi;
  }

  /** Returns the error in the product x=a*b, with precomputed splits for b. */
  private static double twoProductError(double a, double bhi, double blo, double x) {
    double ahi = splitHigh(a);
    double alo = splitLow(a, ahi);
    double err1 = x - (ahi * bhi);
    double err2 = err1 - (alo * bhi);
    double err3 = err2 - (ahi * blo);
    return (alo * blo) - err3;
  }
}
