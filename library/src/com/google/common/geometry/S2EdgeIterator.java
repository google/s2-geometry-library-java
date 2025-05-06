/*
 * Copyright 2025 Google Inc.
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

import com.google.common.geometry.S2Shape.MutableEdge;

/**
 * An iterator that advances through all edges in an S2ShapeIndex. Note that empty and full shapes
 * have no edges and the iterator will not visit them.
 *
 * <p>Example usage:
 *
 * {@snippet :
 * MutableEdge edge = new MutableEdge();
 * for (S2EdgeIterator it = new S2EdgeIterator(index); !it.done(); it.next()) {
 *   it.getEdge(edge);
 *   // use edge
 * }
 * }
 */
public final class S2EdgeIterator {
  // The index being iterated over.
  private final S2ShapeIndex index;
  // The current shape id.
  private int shapeId = -1;
  // The current edge id.
  private int edgeId = -1;
  // The number of edges in the current shape. Loaded when advancing to the next shape.
  private int numEdges = 0;

  /** Constructs an S2EdgeIterator for the given index. */
  public S2EdgeIterator(S2ShapeIndex index) {
    this.index = index;
    // Start the first shape, or possibly detect that we are done.
    next();
  }

  /** Copy constructor, creates a new iterator on the same index and at the same position. */
  public S2EdgeIterator(S2EdgeIterator other) {
    this.index = other.index;
    this.shapeId = other.shapeId;
    this.edgeId = other.edgeId;
    this.numEdges = other.numEdges;
  }

  /** Returns the current shape id. */
  public int shapeId() {
    return shapeId;
  }

  /** Returns the current edge id. */
  public int edgeId() {
    return edgeId;
  }

  /** Fills in the current edge. */
  public void getEdge(MutableEdge edge) {
    index.getShapes().get(shapeId).getEdge(edgeId, edge);
  }

  /** Returns true if there are no more edges to visit. */
  public boolean done() {
    return shapeId >= index.getShapes().size();
  }

  /**
   * Returns true if this iterator is iterating the same index and is at the same position as the
   * other iterator.
   */
  public boolean isEqualTo(S2EdgeIterator other) {
    return index == other.index && shapeId == other.shapeId && edgeId == other.edgeId;
  }

  /** Advances to the next edge. */
  public void next() {
    while (++edgeId >= numEdges) {
      // Advance to the next shape.
      if (++shapeId >= index.getShapes().size()) {
        // There are no more shapes, done() is true.
        break;
      }
      // Get the number of edges in the next shape.
      S2Shape shape = index.getShapes().get(shapeId);
      numEdges = (shape == null) ? 0 : shape.numEdges();
      edgeId = -1;
    }
  }
}
