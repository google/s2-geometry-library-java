/*
 * Copyright 2018 Google Inc.
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

import com.google.common.annotations.GwtIncompatible;
import com.google.common.base.Preconditions;
import com.google.common.collect.ImmutableList;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.primitives.Ints;
import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;

/** An encoder/decoder of tagged {@link S2Shape}s. */
@GwtIncompatible("S2LaxPolylineShape and S2LaxPolygonShape")
public class S2TaggedShapeCoder implements S2Coder<S2Shape> {

  private static final int POLYGON_TYPE_TAG = 1;
  private static final int POLYLINE_TYPE_TAG = 2;
  private static final int POINT_TYPE_TAG = 3;
  private static final int LAX_POLYLINE_TYPE_TAG = 4;
  private static final int LAX_POLYGON_TYPE_TAG = 5;

  private static final S2Coder<S2Polygon.Shape> FAST_POLYGON_SHAPE_CODER =
      new S2Coder<S2Polygon.Shape>() {
        @Override
        public void encode(S2Polygon.Shape value, OutputStream output) throws IOException {
          value.polygon().encodeUncompressed(new LittleEndianOutput(output));
        }

        @Override
        public S2Polygon.Shape decode(Bytes data, Cursor cursor) {
          try {
            return S2Polygon.decode(data.toInputStream(cursor)).shape();
          } catch (IOException e) {
            throw new RuntimeException(e);
          }
        }
      };

  private static final S2Coder<S2Polygon.Shape> COMPACT_POLYGON_SHAPE_CODER =
      new S2Coder<S2Polygon.Shape>() {
        @Override
        public void encode(S2Polygon.Shape value, OutputStream output) throws IOException {
          value.polygon().encode(output);
        }

        @Override
        public S2Polygon.Shape decode(Bytes data, Cursor cursor) {
          return FAST_POLYGON_SHAPE_CODER.decode(data, cursor);
        }
      };

  private static final S2Coder<S2Polyline> FAST_POLYLINE_SHAPE_CODER =
      new S2Coder<S2Polyline>() {
        @Override
        public void encode(S2Polyline value, OutputStream output) throws IOException {
          value.encode(output);
        }

        @Override
        public S2Polyline decode(Bytes data, Cursor cursor) {
          try {
            return S2Polyline.decode(data.toInputStream(cursor));
          } catch (IOException e) {
            throw new RuntimeException(e);
          }
        }
      };

  private static final ImmutableList<Class<? extends S2Polygon.Shape>> POLYGON_SHAPE_CLASSES =
      ImmutableList.of(
          new S2Polygon().binarySearchShape().getClass(),
          new S2Polygon().linearSearchShape().getClass());
  private static final ImmutableList<Class<? extends S2Point.Shape>> POINT_SHAPE_CLASSES =
      ImmutableList.of(
          S2Point.Shape.class,
          S2Point.Shape.singleton(S2Point.ORIGIN).getClass(),
          S2Point.Shape.fromList(ImmutableList.of(S2Point.ORIGIN, S2Point.ORIGIN)).getClass());
  private static final ImmutableList<Class<? extends S2LaxPolylineShape>>
      LAX_POLYLINE_SHAPE_CLASSES =
          ImmutableList.of(
              S2LaxPolylineShape.SimpleArray.class,
              S2LaxPolylineShape.SimpleList.class,
              S2LaxPolylineShape.SimplePacked.class,
              S2LaxPolylineShape.SimpleSnapped.class);
  private static final ImmutableList<Class<? extends S2LaxPolygonShape>> LAX_POLYGON_SHAPE_CLASSES =
      ImmutableList.of(
          S2LaxPolygonShape.SimpleArray.class,
          S2LaxPolygonShape.SimpleList.class,
          S2LaxPolygonShape.SimplePacked.class,
          S2LaxPolygonShape.SimpleSnapped.class,
          S2LaxPolygonShape.MultiArray.class,
          S2LaxPolygonShape.MultiList.class,
          S2LaxPolygonShape.MultiPacked.class,
          S2LaxPolygonShape.MultiSnapped.class);

  /**
   * An instance of a {@code S2TaggedShapeCoder} which encodes/decodes {@link S2Shape}s in the FAST
   * encoding format. The FAST format is optimized for fast encoding/decoding.
   */
  public static final S2TaggedShapeCoder FAST =
      new Builder(true)
          .add(POLYGON_SHAPE_CLASSES, FAST_POLYGON_SHAPE_CODER, POLYGON_TYPE_TAG)
          .add(S2Polyline.class, FAST_POLYLINE_SHAPE_CODER, POLYLINE_TYPE_TAG)
          .add(POINT_SHAPE_CLASSES, S2Point.Shape.Coder.FAST, POINT_TYPE_TAG)
          .add(LAX_POLYLINE_SHAPE_CLASSES, S2LaxPolylineShape.Coder.FAST, LAX_POLYLINE_TYPE_TAG)
          .add(LAX_POLYGON_SHAPE_CLASSES, S2LaxPolygonShape.Coder.FAST, LAX_POLYGON_TYPE_TAG)
          .build();

  /**
   * An instance of a {@code S2TaggedShapeCoder} which encodes/decodes {@link S2Shape}s in the
   * COMPACT encoding format. The COMPACT format is optimized for disk usage and memory footprint.
   */
  public static final S2TaggedShapeCoder COMPACT =
      new Builder(true)
          .add(POLYGON_SHAPE_CLASSES, COMPACT_POLYGON_SHAPE_CODER, POLYGON_TYPE_TAG)
          .add(S2Polyline.class, FAST_POLYLINE_SHAPE_CODER, POLYLINE_TYPE_TAG)
          .add(POINT_SHAPE_CLASSES, S2Point.Shape.Coder.COMPACT, POINT_TYPE_TAG)
          .add(LAX_POLYLINE_SHAPE_CLASSES, S2LaxPolylineShape.Coder.COMPACT, LAX_POLYLINE_TYPE_TAG)
          .add(LAX_POLYGON_SHAPE_CLASSES, S2LaxPolygonShape.Coder.COMPACT, LAX_POLYGON_TYPE_TAG)
          .build();

  private final IdentityHashMap<Class<? extends S2Shape>, Integer> classToTypeTag;
  private final Map<Integer, S2Coder<? extends S2Shape>> typeTagToCoder;

  private S2TaggedShapeCoder(
      IdentityHashMap<Class<? extends S2Shape>, Integer> classToTypeTag,
      Map<Integer, S2Coder<? extends S2Shape>> typeTagToCoder) {
    this.classToTypeTag = classToTypeTag;
    this.typeTagToCoder = typeTagToCoder;
  }

  @Override
  @SuppressWarnings("unchecked") // safe covariant cast
  public void encode(S2Shape value, OutputStream output) throws IOException {
    if (value == null) {
      // A null shape is encoded as 0 bytes.
      return;
    }
    Integer typeTag = classToTypeTag.get(value.getClass());
    Preconditions.checkArgument(
        typeTag != null, "No S2Coder matched S2Shape with type %s", value.getClass().getName());
    EncodedInts.writeVarint64(output, typeTag);
    ((S2Coder<S2Shape>) typeTagToCoder.get(typeTag)).encode(value, output);
  }

  @Override
  public S2Shape decode(Bytes data, Cursor cursor) {
    if (cursor.remaining() == 0) {
      // A null shape is encoded as 0 bytes.
      return null;
    }
    int typeTag = Ints.checkedCast(data.readVarint64(cursor));
    S2Coder<? extends S2Shape> coder = typeTagToCoder.get(typeTag);
    Preconditions.checkArgument(coder != null, "No S2Coder matched type tag %s", typeTag);
    return typeTagToCoder.get(typeTag).decode(data, cursor);
  }

  /** Returns a new {@link Builder}. */
  public static Builder builder() {
    return new Builder(false);
  }

  /** Returns a new {@link Builder} initialized with the current {@link S2TaggedShapeCoder}. */
  public Builder toBuilder() {
    return new Builder(classToTypeTag, typeTagToCoder);
  }

  /** A builder for creating {@link S2TaggedShapeCoder} instances. */
  public static class Builder {
    /** The minimum non-reserved type tag. */
    public static final int MIN_USER_TYPE_TAG = 8192;

    private final boolean allowReservedTags;
    private final IdentityHashMap<Class<? extends S2Shape>, Integer> classToTypeTag;
    private final Map<Integer, S2Coder<? extends S2Shape>> typeTagToCoder;

    private Builder(boolean allowReservedTags) {
      this.allowReservedTags = allowReservedTags;
      classToTypeTag = new IdentityHashMap<>();
      typeTagToCoder = new HashMap<>();
    }

    private Builder(
        IdentityHashMap<Class<? extends S2Shape>, Integer> classToTypeTag,
        Map<Integer, S2Coder<? extends S2Shape>> typeTagToCoder) {
      this.allowReservedTags = false;
      this.classToTypeTag = classToTypeTag;
      this.typeTagToCoder = typeTagToCoder;
    }

    /**
     * Associates {@code clazz} with a unique {@code coder} and {@code typeTag}.
     *
     * <p>If {@code clazz} or {@code typeTag} was already added, an {@link IllegalArgumentException}
     * is thrown.
     */
    <T extends S2Shape> Builder add(Class<? extends T> clazz, S2Coder<T> coder, int typeTag) {
      validateTypeTag(typeTag);
      validateClass(clazz);
      classToTypeTag.put(clazz, typeTag);
      typeTagToCoder.put(typeTag, coder);
      return this;
    }

    /**
     * Same as {@link #add(Class, S2Coder, int)}, but associates all elements of {@code clazzes}
     * with a unique {@code coder} and {@code typeTag}.
     */
    <T extends S2Shape> Builder add(
        List<Class<? extends T>> clazzes, S2Coder<T> coder, int typeTag) {
      validateTypeTag(typeTag);
      for (Class<? extends T> clazz : clazzes) {
        validateClass(clazz);
        classToTypeTag.put(clazz, typeTag);
      }
      typeTagToCoder.put(typeTag, coder);
      return this;
    }

    private void validateTypeTag(int typeTag) {
      Preconditions.checkArgument(
          allowReservedTags || typeTag >= MIN_USER_TYPE_TAG,
          "Type tag must be greater than %s, got: %s",
          MIN_USER_TYPE_TAG,
          typeTag);
      Preconditions.checkArgument(
          !typeTagToCoder.containsKey(typeTag), "Duplicate type tag: %s", typeTag);
    }

    private <T extends S2Shape> void validateClass(Class<? extends T> clazz) {
      Preconditions.checkArgument(
          !classToTypeTag.containsKey(clazz), "Duplicate class: %s", clazz.getName());
    }

    /** Returns a newly-created {@link S2TaggedShapeCoder}. */
    S2TaggedShapeCoder build() {
      return new S2TaggedShapeCoder(classToTypeTag, typeTagToCoder);
    }
  }
}
