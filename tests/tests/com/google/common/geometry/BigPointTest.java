/*
 * Copyright 2018 Google Inc.
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

import static org.junit.Assert.assertEquals;

import java.math.BigDecimal;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

/** Unit tests for {@link BigPoint}. */
@RunWith(JUnit4.class)
public class BigPointTest {
  @Test
  public void testCrossProd() {
    BigPoint a = createBp("1.2", "3.4", "5.6");
    BigPoint b = createBp("-7.8", "9.0", "-1.2");
    assertEquals(createBp("-54.48", "-42.24", "37.32"), a.crossProd(b));
    assertEquals(createBp("54.48", "42.24", "-37.32"), b.crossProd(a));
  }

  @Test
  public void testDotProd() {
    BigPoint a = createBp("1.2", "3.4", "5.6");
    BigPoint b = createBp("-7.8", "9.0", "-1.2");
    assertEquals(new BigDecimal("14.52"), a.dotProd(b));
    assertEquals(new BigDecimal("14.52"), b.dotProd(a));
  }

  @Test
  public void testNormalizedPoint() {
    // Zero vector.
    BigPoint p = createBp("0", "0", "0");
    S2Point normalized = p.normalizedPoint();
    assertEquals(0.0, normalized.norm(), 1e-15);
    assertEquals(new S1Angle(normalized, S2Point.ZERO).radians(), 0.0, 1e-15);

    // Ordinary vector without underflow.
    p = createBp("1", "2", "3");
    normalized = p.normalizedPoint();
    assertEquals(1.0, normalized.norm(), 1e-15);
    assertEquals(new S1Angle(normalized, new S2Point(1, 2, 3)).radians(), 0.0, 1e-15);

    // All underflowing components; normalized vector is nearly Y_NEG.
    p = createBp("1e-400", "-2e-200", "3e-300");
    normalized = p.normalizedPoint();
    assertEquals(1.0, normalized.norm(), 1e-15);
    assertEquals(new S1Angle(normalized, S2Point.Y_NEG).radians(), 0.0, 1e-15);

    // One underflowing component; no scaling needed.
    p = createBp("1", "2e-200", "1");
    normalized = p.normalizedPoint();
    assertEquals(1.0, normalized.norm(), 1e-15);
    assertEquals(new S1Angle(normalized, new S2Point(1, 0, 1)).radians(), 0.0, 1e-15);

    // Two components underflow; scaling needed.
    p = createBp("0", "2e-200", "3e-200");
    normalized = p.normalizedPoint();
    assertEquals(1.0, normalized.norm(), 1e-15);
    assertEquals(new S1Angle(normalized, new S2Point(0, 2, 3)).radians(), 0.0, 1e-15);
  }

  private static BigPoint createBp(String x, String y, String z) {
    return new BigPoint(new BigDecimal(x), new BigDecimal(y), new BigDecimal(z));
  }
}
