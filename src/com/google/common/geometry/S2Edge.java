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

import com.google.common.annotations.GwtCompatible;

import java.io.Serializable;

/**
 * An abstract directed edge from one S2Point to another S2Point.
 *
 * @author kirilll@google.com (Kirill Levin)
 */
@GwtCompatible(serializable = true)
public final class S2Edge implements Serializable {

  private final S2Point start;
  private final S2Point end;

  public S2Edge(S2Point start, S2Point end) {
    this.start = start;
    this.end = end;
  }

  public S2Point getStart() {
    return start;
  }

  public S2Point getEnd() {
    return end;
  }

  @Override
  public String toString() {
    return "Edge: (" + start.toDegreesString() + " -> " + end.toDegreesString() + ")";
  }

  @Override
  public int hashCode() {
    return getStart().hashCode() - getEnd().hashCode();
  }

  @Override
  public boolean equals(Object o) {
    if (o == null || !(o instanceof S2Edge)) {
      return false;
    }
    S2Edge other = (S2Edge) o;
    return getStart().equalsPoint(other.getStart())
        && getEnd().equalsPoint(other.getEnd());
  }
}
