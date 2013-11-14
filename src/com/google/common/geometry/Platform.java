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

import com.google.common.annotations.GwtCompatible;

import java.io.PrintStream;
import java.util.logging.Logger;

/**
 * Contains utility methods which require different GWT client and server implementations.
 * This contains the server side implementations.
 */
@GwtCompatible(emulated = true)
final class Platform {

  private Platform() {}

  /** @see Math#IEEEremainder(double, double) */
  static double IEEEremainder(double f1, double f2) {
    return Math.IEEEremainder(f1, f2);
  }

  /** @see Math#getExponent(double) */
  static int getExponent(double d) {
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
   * Invokes {@code stream.printf} with the arguments. The GWT client just prints the format string
   * and the arguments separately. Using this method is not recommended; you should instead
   * construct strings with normal string concatenation whenever possible, so it will work the same
   * way in normal Java and GWT client versions.
   */
  static void printf(PrintStream stream, String format, Object... params) {
    stream.printf(format, params);
  }

  /**
   * Returns {@code String.format} with the arguments.  The GWT client just returns a string
   * consisting of the format string with the parameters concatenated to the end of it.  Using
   * this method is not recommended; you should instead construct strings with normal string
   * concatenation whenever possible, so it will work the same way in normal Java and GWT client
   * versions.
   */
  static String formatString(String format, Object... params) {
    return String.format(format, params);
  }

  /**
   * A portable way to hash a double value.
   */
  public static long doubleHash(double value) {
    return Double.doubleToLongBits(value);
  }
}
