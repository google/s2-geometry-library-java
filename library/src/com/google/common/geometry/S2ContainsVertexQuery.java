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

import static java.lang.Math.abs;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.OptionalInt;

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
@SuppressWarnings("Assertion")
public class S2ContainsVertexQuery {
  /** The target vertex for which we are determining containment. */
  private S2Point target;

  /** Endpoints of outgoing polygon edges from "target". */
  private final List<S2Point> outgoing = new ArrayList<>();

  /** Starting points of incoming polygon edges to "target". */
  private final List<S2Point> incoming = new ArrayList<>();

  /** Creates a contains vertex query. init() must be called before use. */
  public S2ContainsVertexQuery() {
    this.target = null;
  }

  /** Creates a contains vertex query to determine containment of 'target'. */
  public S2ContainsVertexQuery(S2Point target) {
    this.target = target;
  }

  /** Initializes the query to determine containment of 'target'. */
  public void init(S2Point target) {
    this.target = target;
    outgoing.clear();
    incoming.clear();
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
   * Returns +1 if the target vertex is contained, -1 if it is not contained, and 0 if the incident
   * edges consisted of matched sibling pairs.
   *
   * <p>Throws an assertion error if Java assertions are enabled and duplicate unmatched edges
   * incident on the target are found, which indicates an invalid configuration of edges around the
   * target vertex. Consider {@link #safeContainsSign()} if you are not sure the polygon is valid.
   */
  public int containsSign() {
    OptionalInt s = safeContainsSign();
    assert s.isPresent();
    return s.getAsInt();
  }

  /**
   * Returns true if there are any duplicate edges incident on "target", where matched incoming and
   * outgoing edge pairs cancel out. Duplicate edges indicate that the edges around the target
   * vertex are invalid.
   */
  public boolean duplicateEdges() {
    OptionalInt s = safeContainsSign();
    return s.isEmpty();
  }

  /**
   * Returns an OptionalInt which is empty if duplicate unmatched edges incident on the target are
   * found. Otherwise, the OptionalInt contains +1 if the target vertex is contained, -1 if it is
   * not contained, and 0 if the incident edges consisted of matched sibling pairs.
   *
   * <p>Duplicate unmatched edges incident on the target indicate an invalid configuration of edges
   * around the target vertex.
   */
  public OptionalInt safeContainsSign() {
    // Find the unmatched edge that is immediately clockwise from S2.ortho(target) but not equal to
    // it. The result is +1 iff this edge is outgoing.
    //
    // A loop with consecutive vertices A,B,C contains vertex B if and only if the fixed vector R =
    // S2.ortho(B) is contained by the wedge ABC. The wedge is closed at A and open at C, i.e. the
    // point B is inside the loop if A = R but not if C = R.  This convention is required for
    // compatibility with S2EdgeUtil.vertexCrossing().
    boolean duplicateEdges = false;
    S2Point referenceDir = S2.refDir(target);
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
        v = incoming.get(in);
        direction = -count(incoming, in);
        in -= direction;
      } else if (in == incoming.size()) {
        v = outgoing.get(out);
        direction = count(outgoing, out);
        out += direction;
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
      }
      duplicateEdges |= (abs(direction) > 1);
      if (direction == 0) {
        // This is a "matched" edge.
        continue;
      }
      if (S2Predicates.orderedCCW(referenceDir, bestPoint, v, target)) {
        bestPoint = v;
        bestSum = direction;
      }
    }
    return duplicateEdges ? OptionalInt.empty() : OptionalInt.of(bestSum);
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
