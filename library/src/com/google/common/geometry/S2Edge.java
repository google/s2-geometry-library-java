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

import com.google.common.base.Preconditions;
import java.io.Serializable;
import jsinterop.annotations.JsType;

/**
 * An abstract directed edge from one S2Point to another S2Point.
 *
 * @author kirilll@google.com (Kirill Levin)
 */
@JsType
public final class S2Edge implements Serializable, S2Shape {

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
    return getStart().equalsPoint(other.getStart()) && getEnd().equalsPoint(other.getEnd());
  }

  @Override
  public int numEdges() {
    return 1;
  }

  @Override
  public void getEdge(int index, MutableEdge result) {
    result.set(start, end);
  }

  @Override
  public boolean hasInterior() {
    return false;
  }

  @Override
  public boolean containsOrigin() {
    return false;
  }

  @Override
  public int numChains() {
    return 1;
  }

  @Override
  public int getChainStart(int chainId) {
    Preconditions.checkElementIndex(chainId, numChains());
    return 0;
  }

  @Override
  public int getChainLength(int chainId) {
    Preconditions.checkElementIndex(chainId, numChains());
    return 1;
  }

  @Override
  public void getChainEdge(int chainId, int offset, MutableEdge result) {
    Preconditions.checkElementIndex(offset, getChainLength(chainId));
    result.set(start, end);
  }

  @Override
  public S2Point getChainVertex(int chainId, int edgeOffset) {
    Preconditions.checkElementIndex(edgeOffset, getChainLength(chainId));
    return edgeOffset == 0 ? start : end;
  }

  @Override
  public int dimension() {
    return 1;
  }
}
