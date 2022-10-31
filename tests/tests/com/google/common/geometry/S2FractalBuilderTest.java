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

import static java.lang.Math.PI;
import static java.lang.Math.log;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

import com.google.common.base.Preconditions;
import java.util.Random;

public class S2FractalBuilderTest extends GeometryTestCase {
  /**
   * Constructs a fractal and then computes various metrics (number of vertices, total length,
   * minimum and maximum radius) and verifies that they are within expected tolerances. Essentially
   * this directly verifies that the shape constructed *is* a fractal, i.e. the total length of the
   * curve increases exponentially with the level, while the area bounded by the fractal is more or
   * less constant.
   */
  private void assertFractal(int minLevel, int maxLevel, double dimension) {
    // The radius needs to be fairly small to avoid spherical distortions.
    double nominalRadius = 0.001; // radians, or about 6 km
    double distortionError = 1e-5;

    Random rand = new Random(0);
    S2FractalBuilder fractal = new S2FractalBuilder(rand);
    fractal.setMinLevel(minLevel);
    fractal.setMaxLevel(maxLevel);
    fractal.setFractalDimension(dimension);
    S2Point p = data.getRandomPoint();
    S2Loop loop = fractal.makeLoop(S2.getFrame(p), S1Angle.radians(nominalRadius));
    assertTrue(loop.isValid());

    // If minLevel and maxLevel are not equal, then the number of vertices and the total length of
    // the curve are subject to random variation.  Here we compute an approximation of the standard
    // deviation relative to the mean, noting that most of the variance is due to the random choices
    // about whether to stop subdividing at 'minLevel' or not.  (The random choices at higher levels
    // contribute progressively less and less to the variance.)  The 'relativeError' below
    // corresponds to *one* standard deviation of error; it can be increased to a higher multiple
    // if necessary.
    //
    // Details: Let n = 3 * (4 ^ minLevel) and k = (maxLevel - minLevel + 1).  Each of the 'n'
    // edges at minLevel stops subdividing at that level with probability (1/k). This gives a
    // binomial distribution with mean μ = (n/k) and standard deviation s = sqrt((n/k) * (1 - 1/k)).
    // The relative error (s/μ) can be simplified to sqrt((k - 1)/n).
    int numLevels = maxLevel - minLevel + 1;
    int minVertices = numVerticesAtLevel(minLevel);
    double relativeError = sqrt((numLevels - 1.0) / minVertices);

    // 'expansionFactor' is the total fractal length at level 'n + 1' divided by the total fractal
    // length at level 'n'.
    double expansionFactor = pow(4, 1 - 1 / dimension);
    double expectedNumVertices = 0;
    double expectedLengthSum = 0;

    // 'trianglePerim' is the perimeter of the original equilateral triangle before any
    // subdivision occurs.
    double trianglePerim = 3 * sqrt(3) * nominalRadius;
    double minLengthSum = trianglePerim * pow(expansionFactor, minLevel);
    for (int level = minLevel; level <= maxLevel; ++level) {
      expectedNumVertices += numVerticesAtLevel(level);
      expectedLengthSum += pow(expansionFactor, level);
    }
    expectedNumVertices /= numLevels;
    expectedLengthSum *= trianglePerim / numLevels;

    assertTrue(loop.numVertices() >= minVertices);
    assertTrue(loop.numVertices() <= numVerticesAtLevel(maxLevel));
    assertDoubleNear(
        expectedNumVertices,
        loop.numVertices(),
        relativeError * (expectedNumVertices - minVertices) + 1e-11);

    S2Point center = p;
    double minRadius = 2 * PI;
    double maxRadius = 0;
    double lengthSum = 0;
    for (int i = 0; i < loop.numVertices(); ++i) {
      double r = center.angle(loop.vertex(i));
      minRadius = min(minRadius, r);
      maxRadius = max(maxRadius, r);
      lengthSum += loop.vertex(i).angle(loop.vertex(i + 1));
    }
    // kVertexError is an approximate bound on the error when computing vertex positions of the
    // fractal (trig calculations, etc.)
    double vertexError = 1e-14;

    // Although minRadiusFactor() is only a lower bound in general, it happens to be exact (to
    // within numerical errors) unless the dimension is in the range (1.0, 1.09).
    if (dimension == 1.0 || dimension >= 1.09) {
      // Expect the min radius to match very closely.
      assertDoubleNear(minRadius, fractal.minRadiusFactor() * nominalRadius, vertexError);
    } else {
      // Expect the min radius to satisfy the lower bound.
      assertTrue(minRadius >= fractal.minRadiusFactor() * nominalRadius - vertexError);
    }
    // Check that maxRadiusFactor() is exact (modulo errors) for all dimensions.
    assertDoubleNear(maxRadius, fractal.maxRadiusFactor() * nominalRadius, vertexError);

    assertDoubleNear(
        expectedLengthSum,
        lengthSum,
        relativeError * (expectedLengthSum - minLengthSum) + distortionError * lengthSum);
  }

  private static int numVerticesAtLevel(int level) {
    // Sanity and overflow check.
    Preconditions.checkArgument(level >= 0 && level <= 14);
    // 3 * (4 ^ level)
    return 3 * (1 << (2 * level));
  }

  public void testTriangleFractal() {
    assertFractal(7, 7, 1.0);
  }

  public void testTriangleMultiFractal() {
    assertFractal(1, 6, 1.0);
  }

  public void testKochCurveFractal() {
    assertFractal(7, 7, log(4) / log(3));
  }

  public void testKochCurveMultiFractal() {
    assertFractal(4, 8, log(4) / log(3));
  }

  public void testCesaroFractal() {
    assertFractal(7, 7, 1.8);
  }

  public void testCesaroMultiFractal() {
    assertFractal(2, 6, 1.8);
  }
}
