/*
 * Copyright 2015 Google Inc.
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

/**
 * An interface to specify which vertices must be kept when simplifying a line or loop.
 */
@GwtCompatible
public interface S2VertexFilter {
  /**
   * Returns true if the given vertex must be kept in the simplified loop (returning false does not
   * mean that the vertex *must* be removed, but that it *can* be removed if it satisfies the
   * simplification criteria).
   */
  boolean shouldKeepVertex(S2Point vertex);
}
