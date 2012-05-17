/*
 * Copyright 2006 Google Inc.
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
 * An S2Point that also has a parameter associated with it, which corresponds
 * to a time-like order on the points.
 *
 * @author kirilll@google.com (Kirill Levin)
 */
public final class ParametrizedS2Point implements Comparable<ParametrizedS2Point> {
  private final double time;
  private final S2Point point;

  public ParametrizedS2Point(double time, S2Point point) {
    this.time = time;
    this.point = point;
  }

  public double getTime() {
    return time;
  }

  public S2Point getPoint() {
    return point;
  }

  @Override
  public int compareTo(ParametrizedS2Point o) {
    int compareTime = Double.compare(time, o.time);
    if (compareTime != 0) {
      return compareTime;
    }
    return point.compareTo(o.point);
  }

  @Override
  public boolean equals(Object other) {
    if (other instanceof ParametrizedS2Point) {
      ParametrizedS2Point x = (ParametrizedS2Point) other;
      return time == x.time && point.equals(x.point);
    } else {
      return false;
    }
  }
}
