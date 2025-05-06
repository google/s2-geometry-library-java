/*
 * Copyright 2022 Google Inc.
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
 * An S2BuilderLayer is responsible for assembling an S2BuilderGraph of snapped edges into the
 * desired output format (e.g., an S2Polygon).
 *
 * <p>This interface is not needed by ordinary S2Builder clients. It is only necessary if you want
 * to implement a new S2BuilderLayer subtype.
 */
public interface S2BuilderLayer {
  /**
   * Returns the GraphOptions used for building the edge graph passed to {@link
   * #build(S2BuilderGraph, S2Error)}.
   */
  S2Builder.GraphOptions graphOptions();

  /**
   * Assembles a graph of snapped edges into the geometry type implemented by this layer. If an
   * error is encountered, sets "error" appropriately and returns false, otherwise returns true.
   *
   * <p>Note that when there are multiple layers, the Graph objects passed to all layers are
   * guaranteed to be valid until the last build() method returns. This makes it easier to write
   * algorithms that gather the output graphs from several layers and process them all at once, such
   * as {@link S2ClosedSetNormalizer}.
   */
  boolean build(S2BuilderGraph graph, S2Error error);
}
