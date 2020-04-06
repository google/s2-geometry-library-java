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

import com.google.common.annotations.GwtCompatible;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * This class determines whether a polygon contains one of its vertices given the edges incident to
 * that vertex. The result is +1 if the vertex is contained, -1 if it is not contained, and 0 if the
 * incident edges consist of matched sibling pairs (in which case the result cannot be determined
 * locally).
 *
 * <p>The {@link S2ContainsPointQuery.S2VertexModel#SEMI_OPEN "semi-open" boundary model} is used to
 * define point containment. This means that if several polygons tile the region around a vertex,
 * then exactly one of those polygons contains that vertex.
 *
 * <p>This class is not thread-safe.
 */
@GwtCompatible
public class S2ContainsVertexQuery {
  private final S2Point target;
  private final List<S2Point> outgoing = new ArrayList<>();
  private final List<S2Point> incoming = new ArrayList<>();

  /** Creates a contains vertex query to determine containment of 'target'. */
  public S2ContainsVertexQuery(S2Point target) {
    this.target = target;
  }

  /** Adds an edge outgoing from 'target' to 'v'. */
  public void addOutgoing(S2Point v) {
    outgoing.add(v);
  }

  /** Adds an edge from 'v' incoming to 'target'. */
  public void addIncoming(S2Point v) {
    incoming.add(v);
  }

  /**
   * Returns +1 if the vertex is contained, -1 if it is not contained, and 0 if the incident edges
   * consisted of matched sibling pairs.
   */
  public int containsSign() {
    // Find the unmatched edge that is immediately clockwise from S2.ortho(target).
    S2Point referenceDir = S2.ortho(target);
    S2Point bestPoint = referenceDir;
    int bestSum = 0;
    // Merge outgoing and incoming lists together, computing a sum of each distinct vertex as the
    // count of outgoing occurrences minus the count of incoming occurrences.
    Collections.sort(outgoing);
    Collections.sort(incoming);
    for (int out = 0, in = 0; out < outgoing.size() || in < incoming.size(); ) {
      S2Point v;
      int direction;
      if (out == outgoing.size()) {
        v = incoming.get(in++);
        direction = -1;
      } else if (in == incoming.size()) {
        v = outgoing.get(out++);
        direction = 1;
      } else {
        S2Point outPoint = outgoing.get(out);
        S2Point inPoint = incoming.get(in);
        int diff = outPoint.compareTo(inPoint);
        if (diff < 0) {
          // The out point is smaller, so increase direction by each occurrence.
          v = outPoint;
          direction = count(outgoing, out);
          out += direction;
        } else if (diff > 0) {
          // The in point is smaller, so decrease direction by each occurrence.
          v = inPoint;
          direction = -count(incoming, in);
          in -= direction;
        } else {
          // The points are equal, so increase direction by the difference in counts.
          v = outPoint;
          int outSum = count(outgoing, out);
          int inSum = count(incoming, in);
          direction = outSum - inSum;
          out += outSum;
          in += inSum;
        }
        assert Math.abs(direction) <= 1;
      }
      if (direction == 0) {
        // This is a "matched" edge.
        continue;
      }
      if (S2Predicates.orderedCCW(referenceDir, bestPoint, v, target)) {
        bestPoint = v;
        bestSum = direction;
      }
    }
    return bestSum;
  }

  /** Returns the count of vertices equal to vertices[start]. */
  private static int count(List<S2Point> vertices, int start) {
    S2Point v = vertices.get(start);
    int sum = 1;
    for (int i = start + 1; i < vertices.size() && vertices.get(i).equalsPoint(v); i++) {
      sum++;
    }
    return sum;
  }
}
