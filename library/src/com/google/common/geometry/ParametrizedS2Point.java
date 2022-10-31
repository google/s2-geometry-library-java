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

import java.io.Serializable;

/**
 * An S2Point that also has a parameter associated with it, which corresponds to a time-like order
 * on the points. Comparing ParametrizedS2Points uses the time parameter, with the point as a
 * tiebreaker.
 *
 * @author kirilll@google.com (Kirill Levin)
 */
public final class ParametrizedS2Point implements Comparable<ParametrizedS2Point>, Serializable {

  private final double time;
  private final S2Point point;

  @SuppressWarnings("GoodTime") // should accept a java.time.Duration (?)
  public ParametrizedS2Point(double time, S2Point point) {
    this.time = time;
    this.point = point;
  }

  @SuppressWarnings("GoodTime") // should return a java.time.Duration (?)
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
      return time == x.time && point.equalsPoint(x.point);
    } else {
      return false;
    }
  }

  @Override
  public int hashCode() {
    // TODO(user): When this class or the library overall moves to SDK version 19 or higher,
    // use Objects.hash() here. The Android SDK version for this class is currently the default
    // (14), which does not support Objects.hash. Double.hashCode requires an even higher API level,
    // so just hash the point.
    return point.hashCode();
  }
}
