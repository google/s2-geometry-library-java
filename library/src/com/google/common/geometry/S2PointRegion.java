/*
 * Copyright 2017 Google Inc.
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

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;

/**
 * An S2PointRegion is a region that contains a single point. It is more expensive than the raw
 * S2Point type and is useful mainly for completeness.
 */
public final class S2PointRegion implements S2Region, Comparable<S2PointRegion>, Serializable {
  /** The byte in a stream that signifies the lossless encoding of an S2PointRegion follows. */
  private static final byte POINT_REGION_LOSSLESS_ENCODING_VERSION = 1;

  private final S2Point point;

  public S2PointRegion(double x, double y, double z) {
    this.point = new S2Point(x, y, z);
  }

  public S2PointRegion(S2Point point) {
    this.point = point;
  }

  public S2Point getPoint() {
    return point;
  }

  public double getX() {
    return point.getX();
  }

  public double getY() {
    return point.getY();
  }

  public double getZ() {
    return point.getZ();
  }

  @Override
  public boolean equals(Object that) {
    if (!(that instanceof S2PointRegion)) {
      return false;
    }
    S2PointRegion thatPointRegion = (S2PointRegion) that;
    return getPoint().equalsPoint(thatPointRegion.getPoint());
  }

  public boolean lessThan(S2PointRegion vb) {
    return getPoint().lessThan(vb.getPoint());
  }

  // Required for Comparable
  @Override
  public int compareTo(S2PointRegion other) {
    return (lessThan(other) ? -1 : (equals(other) ? 0 : 1));
  }

  @Override
  public String toString() {
    return point.toString();
  }

  public String toDegreesString() {
    return point.toDegreesString();
  }

  /**
   * Calculates hashcode based on stored coordinates. Since we want +0.0 and -0.0 to be treated the
   * same, we ignore the sign of the coordinates.
   */
  @Override
  public int hashCode() {
    return point.hashCode();
  }

  // S2Region implementation
  @Override
  public boolean contains(S2Cell cell) {
    return false;
  }

  @Override
  public boolean contains(S2Point p) {
    return getPoint().contains(p);
  }

  @Override
  public S2Cap getCapBound() {
    return S2Cap.fromAxisHeight(getPoint(), 0);
  }

  @Override
  public S2LatLngRect getRectBound() {
    S2LatLng latLng = new S2LatLng(getPoint());
    return S2LatLngRect.fromPoint(latLng);
  }

  @Override
  public boolean mayIntersect(S2Cell cell) {
    return cell.contains(getPoint());
  }

  /** Writes this point region to the given output stream. */
  public void encode(OutputStream os) throws IOException {
    encode(new LittleEndianOutput(os));
  }

  /** Writes this point region to the given little endian output stream. */
  void encode(LittleEndianOutput os) throws IOException {
    os.writeByte(POINT_REGION_LOSSLESS_ENCODING_VERSION);
    point.encode(os);
  }

  /** Returns a new S2PointRegion decoded from the given input stream. */
  public static S2PointRegion decode(InputStream is) throws IOException {
    return decode(new LittleEndianInput(is));
  }

  /** Returns a new S2PointRegion decoded from the given little endian input stream. */
  static S2PointRegion decode(LittleEndianInput is) throws IOException {
    byte version = is.readByte();
    if (version != POINT_REGION_LOSSLESS_ENCODING_VERSION) {
      throw new IOException("Unsupported S2PointRegion encoding version " + version);
    }
    return new S2PointRegion(S2Point.decode(is));
  }
}
