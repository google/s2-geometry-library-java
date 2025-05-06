/*
 * Copyright 2023 Google Inc.
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
 * An S2BuilderShapesLayer is an S2BuilderLayer that assembles an S2BuilderGraph of snapped edges
 * into one or more S2Shapes. Implementations of S2BuilderShapesLayer can be used with {@link
 * S2BuilderUtil.IndexedLayer} to automatically add output shapes to an S2ShapeIndex when build()
 * completes.
 */
public interface S2BuilderShapesLayer extends S2BuilderLayer {

  /**
   * After build() has returned true, returns one or more S2Shapes containing the geometry assembled
   * by this layer. Most S2BuilderLayer implementations should also provide a method that returns
   * the specific type of geometry built by that particular layer.
   */
  Iterable<? extends S2Shape> shapes();
}
