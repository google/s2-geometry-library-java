/*
 * Copyright 2019 Google Inc.
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

import static java.lang.Math.max;


/**
 * Defines various angle and area measures for {@link S2ShapeIndex} objects. In general, these
 * methods return the sum of the corresponding measure for all {@link S2Shape} in the index.
 */
public final class S2ShapeIndexMeasures {

  private S2ShapeIndexMeasures() {}

  /**
   * Returns the maximum dimension of any shape in shapeIndex, or -1 if shapeIndex has no shapes.
   *
   * <p>The dimension does <b>not</b> depend on whether the shapes in shapeIndex contain any points.
   * For example, the dimension of an empty point set is 0, and the dimension of an empty polygon is
   * 2.
   */
  public static int dimension(S2ShapeIndex shapeIndex) {
    int dimension = -1;
    for (S2Shape shape : shapeIndex.getShapes()) {
      dimension = max(dimension, shape.dimension());
    }
    return dimension;
  }

  /**
   * Returns the total length of all polylines in shapeIndex, or {@link S1Angle#ZERO} if shapeIndex
   * contains no polylines.
   */
  public static S1Angle length(S2ShapeIndex shapeIndex) {
    S1Angle.Builder builder = new S1Angle.Builder();
    for (S2Shape shape : shapeIndex.getShapes()) {
      builder.add(S2ShapeMeasures.length(shape));
    }
    return builder.build();
  }

  /**
   * Returns the total perimeter of all polygons in shapeIndex (including both "shells" and
   * "holes"), or {@link S1Angle#ZERO} shapeIndex contains no polygons.
   */
  public static S1Angle perimeter(S2ShapeIndex shapeIndex) {
    S1Angle.Builder builder = new S1Angle.Builder();
    for (S2Shape shape : shapeIndex.getShapes()) {
      builder.add(S2ShapeMeasures.perimeter(shape));
    }
    return builder.build();
  }

  /**
   * Returns the total area of all polygons in shapeIndex. Returns 0 if no polygons are present.
   * This method has good relative accuracy for both very large and very small regions. Note that
   * the result may exceed 4*Pi if shapeIndex contains overlapping polygons.
   */
  public static double area(S2ShapeIndex shapeIndex) {
    double area = 0;
    for (S2Shape shape : shapeIndex.getShapes()) {
      area += S2ShapeMeasures.area(shape);
    }
    return area;
  }

  /**
   * Returns the centroid of all shapes whose dimension is maximal within shapeIndex, multiplied by
   * the measure of those shapes. For example, if shapeIndex contains points and polylines, then the
   * result is defined as the centroid of the polylines multiplied by the total length of those
   * polylines. The points would be ignored when computing the centroid.
   *
   * <p>The measure of a given shape is defined as follows:
   *
   * <ul>
   *   <li>For dimension 0 shapes, the measure is {@link S2Shape#numEdges()}.
   *   <li>For dimension 1 shapes, the measure is {@link #length(S2ShapeIndex)}.
   *   <li>For dimension 2 shapes, the measure is {@link #area(S2ShapeIndex)}.
   * </ul>
   *
   * <p>The returned centroid is not unit length, so {@link S2Point#normalize()} may need to be
   * called before passing it to other S2 functions. (0, 0, 0) is returned if the index contains no
   * geometry.
   *
   * <p>The centroid is scaled by the total measure of the shapes for two reasons:
   *
   * <ol>
   *   <li>It is cheaper to compute this way.
   *   <li>This makes it easier to compute the centroid of a collection of shapes (since the
   *       individual centroids can simply be summed)
   * </ol>
   */
  public static S2Point centroid(S2ShapeIndex shapeIndex) {
    int dimension = dimension(shapeIndex);
    S2Point.Builder builder = new S2Point.Builder();
    for (S2Shape shape : shapeIndex.getShapes()) {
      if (shape.dimension() == dimension) {
        builder.add(S2ShapeMeasures.centroid(shape));
      }
    }
    return builder.build();
  }
}
