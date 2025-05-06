/*
 * Copyright 2011 Google Inc.
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

import com.google.common.base.Objects;
import java.io.Serializable;
import jsinterop.annotations.JsType;
import org.jspecify.annotations.Nullable;

/**
 * The area of an interior, i.e. the region on the left side of an odd number of loops and
 * optionally a centroid. The area is between 0 and 4*Pi. If it has a centroid, it is the true
 * centroid of the interior multiplied by the area of the shape. Note that the centroid may not be
 * contained by the shape.
 *
 * @author dbentley@google.com (Daniel Bentley)
 */
@JsType
public final class S2AreaCentroid implements Serializable {

  private final double area;
  private final @Nullable S2Point centroid;

  /** Constructs a new S2AreaCentroid with an area and optional centroid. */
  public S2AreaCentroid(double area, @Nullable S2Point centroid) {
    this.area = area;
    this.centroid = centroid;
  }

  /**
   * Returns the area of a shape interior, i.e. the region on the left side of an odd number of
   * loops.
   */
  public double getArea() {
    return area;
  }

  /**
   * Returns the true centroid of a shape, scaled by the area of the shape. Note that this not a
   * unit-length vector. The centroid might not be contained by the shape.
   */
  public @Nullable S2Point getCentroid() {
    return centroid;
  }

  @Override
  public boolean equals(Object obj) {
    if (obj instanceof S2AreaCentroid) {
      S2AreaCentroid that = (S2AreaCentroid) obj;
      return this.area == that.area && Objects.equal(this.centroid, that.centroid);
    }
    return false;
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(area, centroid);
  }
}
