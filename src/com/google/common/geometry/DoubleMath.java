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

/**
 * Defines the Java equivalent of a couple of advanced floating point functions
 * that are available in C.
 *
 */
strictfp class DoubleMath {
  /** Number of significant digits in a double */
  private static final int DIGITS = 52;

  /** A class that represents a double and a magnitude */
  public static class MantissaExponent implements Comparable<MantissaExponent> {
    public double mantissa;
    public int exp;

    /** No instantiation allowed */
    private MantissaExponent() {
    }

    @Override
    public int compareTo(MantissaExponent dm) {
      return Math.signum(mantissa) * exp < Math.signum(dm.mantissa) * dm.exp ? -1 :
          Math.signum(mantissa) * exp > Math.signum(dm.mantissa) * dm.exp ? 1 : mantissa
          < dm.mantissa ? -1 : mantissa > dm.mantissa ? 1 : 0;
    }
  }

  /**
   * If v is non-zero return an integer exp, so that (0.5 <= |v|*2^(-exp) < 1).
   * this is analogous to the integer part of the return value frexp If v is
   * zero return 0.
   */
  public static int getExp(double v) {
    if (v == 0) {
      return 0;
    }
    return (int) ((0x7ff0000000000000L & Double.doubleToLongBits(v)) >> DIGITS) - 1022;
  }


  /**
   * As in C++'s <code>double frexp ( double x , int * exp )</code> from math.h,
   * this function separates the mantissa and exponent of a floating-point
   * value.
   *
   *  This code certainly does not handle java's non-numerical values (NaN and
   * the like).
   */
  public static MantissaExponent frexp(double v) {
    MantissaExponent dm = new MantissaExponent();
    if (v == 0) {
      dm.mantissa = 0;
      dm.exp = 0;
      return dm;
    }
    long bits = Double.doubleToLongBits(v);
    dm.mantissa = Double.longBitsToDouble((0x800fffffffffffffL & bits) | 0x3fe0000000000000L);
    dm.exp = (int) ((0x7ff0000000000000L & bits) >> DIGITS) - 1022;
    return dm;
  }

}
