/*
 * Copyright 2006 Google Inc.
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

import static java.lang.Math.PI;
import static java.lang.Math.asin;
import static java.lang.Math.atan2;
import static java.lang.Math.cos;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sin;

import com.google.common.base.Objects;
import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import com.google.common.geometry.S2Projections.FaceSiTi;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.logging.Logger;
import jsinterop.annotations.JsConstructor;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsType;

/**
 * An S2Polyline represents a sequence of zero or more vertices connected by straight edges
 * (geodesics). Edges of length 0 and 180 degrees are not allowed, i.e. adjacent vertices should not
 * be identical or antipodal.
 *
 * <p>Note: Polylines do not have a Contains(S2Point) method, because "containment" is not
 * numerically well-defined except at the polyline vertices.
 *
 * @author shakusa@google.com (Steven Hakusa) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@JsType
public final strictfp class S2Polyline implements S2Shape, S2Region, Serializable {
  private static final Logger log = Platform.getLoggerForClass(S2Polyline.class);

  private static final S2Point[] ARR_TEMPLATE = {};
  private static final byte LOSSLESS_ENCODING_VERSION = 1;
  private static final byte COMPRESSED_ENCODING_VERSION = 2;

  /** A fast {@link S2Coder} of polylines that uses {@link #encode} and {@link #decode}. */
  public static final S2Coder<S2Polyline> FAST_CODER = new S2Coder<S2Polyline>() {
    @Override public void encode(S2Polyline value, OutputStream output) throws IOException {
      value.encodeUncompressed(new LittleEndianOutput(output));
    }

    @Override public S2Polyline decode(Bytes data, Cursor cursor) throws IOException {
      return S2Polyline.decode(data.toInputStream(cursor));
    }

    @Override public boolean isLazy() {
      return false;
    }
  };

  /** A fast {@link S2Coder} of polylines that uses {@link #encode} and {@link #decode}. */
  public static final S2Coder<S2Polyline> COMPACT_CODER = new S2Coder<S2Polyline>() {
    @Override public void encode(S2Polyline value, OutputStream output) throws IOException {
      value.encodeCompact(output);
    }

    @Override public S2Polyline decode(Bytes data, Cursor cursor) throws IOException {
      return S2Polyline.decode(data.toInputStream(cursor));
    }

    @Override public boolean isLazy() {
      return false;
    }
  };

  private final int numVertices;
  private final S2Point[] vertices;

  /**
   * Create a polyline that connects the given vertices. Empty polylines are allowed. Adjacent
   * vertices should not be identical or antipodal. All vertices should be unit length.
   */
  @JsIgnore
  public S2Polyline(List<S2Point> vertices) {
    this(vertices.toArray(ARR_TEMPLATE));
  }

  @JsConstructor
  private S2Polyline(S2Point[] vertices) {
    // assert isValid(vertices);
    this.numVertices = vertices.length;
    this.vertices = vertices;
  }

  /** Returns an unmodifiable view of the vertices of this polyline. */
  public List<S2Point> vertices() {
    return Collections.unmodifiableList(Arrays.asList(vertices));
  }

  /**
   * Return true if the polyline is valid having all vertices be in unit length and having no
   * identical or antipodal adjacent vertices.
   */
  public boolean isValid() {
    return isValid(vertices());
  }

  /** Return true if the given vertices form a valid polyline. */
  @JsIgnore
  public boolean isValid(List<S2Point> vertices) {
    // All vertices must be unit length.
    int n = vertices.size();
    for (int i = 0; i < n; ++i) {
      if (!S2.isUnitLength(vertices.get(i))) {
        log.info("Vertex " + i + " is not unit length");
        return false;
      }
    }

    // Adjacent vertices must not be identical or antipodal.
    for (int i = 1; i < n; ++i) {
      if (vertices.get(i - 1).equalsPoint(vertices.get(i))
          || vertices.get(i - 1).equalsPoint(S2Point.neg(vertices.get(i)))) {
        log.info("Vertices " + (i - 1) + " and " + i + " are identical or antipodal");
        return false;
      }
    }

    return true;
  }

  public int numVertices() {
    return numVertices;
  }

  public S2Point vertex(int k) {
    // assert (k >= 0 && k < numVertices);
    return vertices[k];
  }

  /** Return the angle corresponding to the total arclength of the polyline on a unit sphere. */
  public S1Angle getArclengthAngle() {
    double lengthSum = 0;
    for (int i = 1; i < numVertices(); ++i) {
      lengthSum += vertex(i - 1).angle(vertex(i));
    }
    return S1Angle.radians(lengthSum);
  }

  /**
   * Return the true centroid of the polyline multiplied by the length of the polyline. (See the
   * "About centroids" comments of S2.java for details on centroids). The result is not unit length,
   * so you may want to normalize it.
   *
   * <p>Prescaling by the polyline length makes it easy to compute the centroid of several polylines
   * by simply adding up their centroids.
   */
  public S2Point getCentroid() {
    return S2ShapeMeasures.polylineCentroid(this, /* chainId= */ 0);
  }

  /**
   * Return the point whose distance from vertex 0 along the polyline is the given fraction of the
   * polyline's total length. Fractions less than zero or greater than one are clamped. The return
   * value is unit length. This cost of this function is currently linear in the number of vertices.
   */
  public S2Point interpolate(double fraction) {
    // We intentionally let the (fraction >= 1) case fall through, since we need to handle it in the
    // loop below in any case because of possible roundoff errors.
    if (fraction <= 0) {
      return vertex(0);
    }

    double lengthSum = 0;
    for (int i = 1; i < numVertices(); ++i) {
      lengthSum += vertex(i - 1).angle(vertex(i));
    }
    double target = fraction * lengthSum;
    for (int i = 1; i < numVertices(); ++i) {
      double length = vertex(i - 1).angle(vertex(i));
      if (target < length) {
        // This code interpolates with respect to arc length rather than straight-line distance, and
        // produces a unit-length result.
        double f = sin(target) / sin(length);
        return vertex(i - 1).mul(cos(target) - f * cos(length)).add(vertex(i).mul(f));
      }
      target -= length;
    }
    return vertex(numVertices() - 1);
  }

  /**
   * Projects the query point to the nearest part of the polyline, and returns the fraction of the
   * polyline's total length traveled along the polyline from vertex 0 to the projected point.
   *
   * <p>For any query point, the returned fraction is at least 0 (when the query point projects to
   * the first vertex of the line) and at most 1 (when the query point projects to the last vertex).
   *
   * <p>This method is essentially the inverse of {@link #interpolate(double)}, except that this
   * method accepts any normalized point, whereas interpolate() only produces points on the line.
   *
   * <p>In the unusual case of multiple equidistant points on the polyline, one of the nearest
   * points is selected in a deterministic but unpredictable manner, and the fraction is computed up
   * to that position. For example, all points of the S2 edge from (1,0,0) to (0,1,0) are
   * equidistant from (0,0,1), so any fraction from 0 to 1 is a correct answer!
   */
  public double uninterpolate(S2Point queryPoint) {
    int i = getNearestEdgeIndex(queryPoint);

    double arcLength =
        S2EdgeUtil.getClosestPoint(queryPoint, vertex(i), vertex(i + 1)).angle(vertex(i));
    for (; i > 0; i--) {
      arcLength += vertex(i - 1).angle(vertex(i));
    }

    return min(arcLength / getArclengthAngle().radians(), 1);
  }

  // S2Region interface (see S2Region.java for details):

  /** Return a bounding spherical cap. */
  @Override
  public S2Cap getCapBound() {
    return getRectBound().getCapBound();
  }

  /** Return a bounding latitude-longitude rectangle. */
  @Override
  public S2LatLngRect getRectBound() {
    S2EdgeUtil.RectBounder bounder = new S2EdgeUtil.RectBounder();
    for (int i = 0; i < numVertices(); ++i) {
      bounder.addPoint(vertex(i));
    }
    return bounder.getBound();
  }

  /**
   * If this method returns true, the region completely contains the given cell. Otherwise, either
   * the region does not contain the cell or the containment relationship could not be determined.
   */
  @Override
  public boolean contains(S2Cell cell) {
    return false;
  }

  @Override
  public boolean contains(S2Point point) {
    return false;
  }

  /**
   * If this method returns false, the region does not intersect the given cell. Otherwise, either
   * region intersects the cell, or the intersection relationship could not be determined.
   */
  @Override
  public boolean mayIntersect(S2Cell cell) {
    if (numVertices() == 0) {
      return false;
    }

    // We only need to check whether the cell contains vertex 0 for correctness, but these tests
    // are cheap compared to edge crossings so we might as well check all the vertices.
    for (int i = 0; i < numVertices(); ++i) {
      if (cell.contains(vertex(i))) {
        return true;
      }
    }
    S2Point[] cellVertices = new S2Point[4];
    for (int i = 0; i < 4; ++i) {
      cellVertices[i] = cell.getVertex(i);
    }
    for (int j = 0; j < 4; ++j) {
      S2EdgeUtil.EdgeCrosser crosser =
          new S2EdgeUtil.EdgeCrosser(cellVertices[j], cellVertices[(j + 1) & 3], vertex(0));
      for (int i = 1; i < numVertices(); ++i) {
        if (crosser.robustCrossing(vertex(i)) >= 0) {
          // There is a proper crossing, or two vertices were the same.
          return true;
        }
      }
    }
    return false;
  }

  /**
   * Returns a new polyline where the vertices of the given polyline have been snapped to the
   * centers of cells at the specified level.
   */
  public static S2Polyline fromSnapped(final S2Polyline a, int snapLevel) {
    // TODO(user): Use S2Builder when available.
    List<S2Point> snappedVertices = new ArrayList<>(a.numVertices());
    S2Point prev = snapPointToLevel(a.vertex(0), snapLevel);
    snappedVertices.add(prev);
    for (int i = 1; i < a.numVertices(); i++) {
      S2Point curr = snapPointToLevel(a.vertex(i), snapLevel);
      if (!curr.equalsPoint(prev)) {
        prev = curr;
        snappedVertices.add(curr);
      }
    }
    return new S2Polyline(snappedVertices);
  }

  /**
   * Returns a new point, snapped to the center of the cell containing the given point at the
   * specified level.
   */
  private static S2Point snapPointToLevel(final S2Point p, int level) {
    return S2CellId.fromPoint(p).parent(level).toPoint();
  }

  /**
   * Return a subsequence of vertex indices such that the polyline connecting these vertices is
   * never further than "tolerance" from the original polyline. Provided the first and last vertices
   * are distinct, they are always preserved; if they are not, the subsequence may contain only a
   * single index.
   *
   * <p>Some useful properties of the algorithm:
   *
   * <ul>
   *   <li>It runs in linear time.
   *   <li>The output is always a valid polyline. In particular, adjacent output vertices are never
   *       identical or antipodal.
   *   <li>The method is not optimal, but it tends to produce 2-3% fewer vertices than the
   *       Douglas-Peucker algorithm with the same tolerance.
   *   <li>The output is *parametrically* equivalent to the original polyline to within the given
   *       tolerance. For example, if a polyline backtracks on itself and then proceeds onwards, the
   *       backtracking will be preserved (to within the given tolerance). This is different than
   *       the Douglas-Peucker algorithm, which only guarantees geometric equivalence.
   * </ul>
   */
  public S2Polyline subsampleVertices(S1Angle tolerance) {
    if (vertices.length == 0) {
      return this;
    }
    List<S2Point> results = Lists.newArrayList();
    results.add(vertex(0));
    S1Angle clampedTolerance = S1Angle.max(tolerance, S1Angle.ZERO);
    for (int i = 0; i < vertices.length - 1; ) {
      int nextIndex = findEndVertex(clampedTolerance, i);
      // Don't create duplicate adjacent vertices.
      if (!vertex(nextIndex).equalsPoint(vertex(i))) {
        results.add(vertex(nextIndex));
      }
      i = nextIndex;
    }
    return new S2Polyline(results);
  }

  /**
   * Given a polyline, a tolerance distance, and a start index, this function returns the maximal
   * end index such that the line segment between these two vertices passes within "tolerance" of
   * all interior vertices, in order.
   */
  private int findEndVertex(S1Angle tolerance, int index) {
    // assert tolerance.radians() >= 0;
    // assert index + 1 < polyline.num_vertices();

    // The basic idea is to keep track of the "pie wedge" of angles from the starting vertex such
    // that a ray from the starting vertex at that angle will pass through the discs of radius
    // "tolerance" centered around all vertices processed so far.

    // First we define a "coordinate frame" for the tangent and normal spaces at the starting
    // vertex.  Essentially this means picking three orthonormal vectors X,Y,Z such that X and Y
    // span the tangent plane at the starting vertex, and Z is "up".  We use the coordinate frame to
    // define a mapping from 3D direction vectors to a one-dimensional "ray angle" in the range
    // (-Pi, Pi].  The angle of a direction vector is computed by transforming it into the X,Y,Z
    // basis, and then calculating atan2(y,x).  This mapping allows us to represent a wedge of
    // angles as a 1D interval.  Since the interval wraps around, we represent it as an S1Interval,
    // i.e. an interval on the unit circle.
    S2Point origin = vertex(index);
    Matrix frame = S2.getFrame(origin);

    // As we go along, we keep track of the current wedge of angles and the distance to the last
    // vertex (which must be non-decreasing).
    S1Interval currentWedge = S1Interval.full();
    double lastDistance = 0;

    for (++index; index < vertices.length; ++index) {
      S2Point candidate = vertex(index);
      double distance = origin.angle(candidate);

      // We don't allow simplification to create edges longer than 90 degrees, to avoid numeric
      // instability as lengths approach 180 degrees.  (We do need to allow for original edges
      // longer than 90 degrees, though.)
      if (distance > PI / 2 && lastDistance > 0) {
        break;
      }

      // Vertices must be in increasing order along the ray, except for the initial disc around the
      // origin.
      if (distance < lastDistance && lastDistance > tolerance.radians()) {
        break;
      }
      lastDistance = distance;

      // Points that are within the tolerance distance of the origin do not constrain the ray
      // direction, so we can ignore them.
      if (distance <= tolerance.radians()) {
        continue;
      }

      // If the current wedge of angles does not contain the angle to this vertex, then stop right
      // now.  Note that the wedge of possible ray angles is not necessarily empty yet, but we can't
      // continue unless we are willing to backtrack to the last vertex that was contained within
      // the wedge (since we don't create new vertices).  This would be more complicated and also
      // make the worst-case running time more than linear.
      S2Point direction = S2.toFrame(frame, candidate);
      double center = atan2(direction.y, direction.x);
      if (!currentWedge.contains(center)) {
        break;
      }

      // To determine how this vertex constrains the possible ray angles, consider the triangle ABC
      // where A is the origin, B is the candidate vertex, and C is one of the two tangent points
      // between A and the spherical cap of radius "tolerance" centered at B.  Then from the
      // spherical law of sines, sin(a)/sin(A) = sin(c)/sin(C), where "a" and "c" are the lengths of
      // the edges opposite A and C.  In our case C is a 90 degree angle, therefore
      // A = asin(sin(a) / sin(c)).  Angle A is the half-angle of the allowable wedge.
      double halfAngle = asin(sin(tolerance.radians()) / sin(distance));
      S1Interval target = S1Interval.fromPoint(center).expanded(halfAngle);
      currentWedge = currentWedge.intersection(target);
      // assert !currentWedge.isEmpty();
    }

    // We break out of the loop when we reach a vertex index that can't be included in the line
    // segment, so back up by one vertex.
    return index - 1;
  }

  /**
   * Given a point, returns the index of the start point of the (first) edge on the polyline that is
   * closest to the given point. The polyline must have at least one vertex. Throws
   * IllegalStateException if this is not the case.
   */
  public int getNearestEdgeIndex(S2Point point) {
    Preconditions.checkState(numVertices() > 0, "Empty polyline");

    if (numVertices() == 1) {
      // If there is only one vertex, the "edge" is trivial, and it's the only one
      return 0;
    }

    // Initial value larger than any possible distance on the unit sphere.
    S1Angle minDistance = S1Angle.radians(10);
    int minIndex = -1;

    // Find the line segment in the polyline that is closest to the point given.
    for (int i = 0; i < numVertices() - 1; ++i) {
      S1Angle distanceToSegment = S2EdgeUtil.getDistance(point, vertex(i), vertex(i + 1));
      if (distanceToSegment.lessThan(minDistance)) {
        minDistance = distanceToSegment;
        minIndex = i;
      }
    }
    return minIndex;
  }

  /**
   * Given a point p and the index of the start point of an edge of this polyline, returns the point
   * on that edge that is closest to p.
   */
  public S2Point projectToEdge(S2Point point, int index) {
    Preconditions.checkState(numVertices() > 0, "Empty polyline");
    Preconditions.checkState(numVertices() == 1 || index < numVertices() - 1, "Invalid edge index");
    if (numVertices() == 1) {
      // If there is only one vertex, it is always closest to any given point.
      return vertex(0);
    }
    return S2EdgeUtil.getClosestPoint(point, vertex(index), vertex(index + 1));
  }

  /**
   * Returns the point on the polyline closest to {@code queryPoint}.
   *
   * <p>In the unusual case of a query point that is equidistant from multiple points on the line,
   * one is returned in a deterministic but otherwise unpredictable way.
   */
  public S2Point project(S2Point queryPoint) {
    Preconditions.checkState(numVertices() > 0, "Empty polyline");
    if (numVertices() == 1) {
      // If there is only one vertex, it is always closest to any given point.
      return vertex(0);
    }
    int i = getNearestEdgeIndex(queryPoint);
    return S2EdgeUtil.getClosestPoint(queryPoint, vertex(i), vertex(i + 1));
  }

  @Override
  public boolean equals(Object that) {
    if (!(that instanceof S2Polyline)) {
      return false;
    }

    S2Polyline thatPolyline = (S2Polyline) that;
    if (numVertices != thatPolyline.numVertices) {
      return false;
    }

    for (int i = 0; i < vertices.length; i++) {
      if (!vertices[i].equalsPoint(thatPolyline.vertices[i])) {
        return false;
      }
    }
    return true;
  }

  /**
   * Return true if this polyline intersects the given polyline. If the polylines share a vertex
   * they are considered to be intersecting. When a polyline endpoint is the only intersection with
   * the other polyline, the function may return true or false arbitrarily.
   *
   * <p>The running time is quadratic in the number of vertices.
   */
  public boolean intersects(S2Polyline line) {
    if (numVertices() <= 0 || line.numVertices() <= 0) {
      return false;
    }

    if (!getRectBound().intersects(line.getRectBound())) {
      return false;
    }

    // TODO(user): Use S2ShapeIndex here.
    for (int i = 1; i < numVertices(); ++i) {
      EdgeCrosser crosser = new EdgeCrosser(vertex(i - 1), vertex(i), line.vertex(0));
      for (int j = 1; j < line.numVertices(); ++j) {
        if (crosser.robustCrossing(line.vertex(j)) >= 0) {
          return true;
        }
      }
    }
    return false;
  }

  @Override
  public int hashCode() {
    return Objects.hashCode(numVertices, Arrays.deepHashCode(vertices));
  }

  @Override
  public String toString() {
    StringBuilder builder = new StringBuilder("S2Polyline, ");

    builder.append(vertices.length).append(" points. [");

    for (S2Point v : vertices) {
      builder.append(v.toDegreesString()).append(" ");
    }
    builder.append("]");

    return builder.toString();
  }

  // S2Shape interface (see S2Shape.java for details):

  @Override
  public int numEdges() {
    return max(0, numVertices - 1);
  }

  @Override
  public void getEdge(int index, MutableEdge result) {
    result.set(vertices[index], vertices[index + 1]);
  }

  @Override
  public boolean hasInterior() {
    return false;
  }

  @Override
  public boolean containsOrigin() {
    throw new IllegalStateException(
        "An S2Polyline has no interior, so " + "containsOrigin() should never be called on one.");
  }

  @Override
  public int numChains() {
    return (numVertices > 1) ? 1 : 0;
  }

  @Override
  public int getChainStart(int chainId) {
    Preconditions.checkElementIndex(chainId, numChains());
    return 0;
  }

  @Override
  public int getChainLength(int chainId) {
    Preconditions.checkElementIndex(chainId, numChains());
    return numEdges();
  }

  @Override
  public void getChainEdge(int chainId, int offset, MutableEdge result) {
    Preconditions.checkElementIndex(chainId, numChains());
    getEdge(offset, result);
  }

  @Override
  public S2Point getChainVertex(int chainId, int edgeOffset) {
    Preconditions.checkElementIndex(chainId, numChains());
    return vertex(edgeOffset);
  }

  @Override
  public int dimension() {
    return 1;
  }

  /** Encodes this polyline into the given output stream. */
  @JsIgnore // OutputStream is not available to J2CL.
  public void encode(OutputStream os) throws IOException {
    encodeUncompressed(new LittleEndianOutput(os));
  }

  /**
   * Encodes the polyline into an efficient, lossless binary representation, which can be decoded by
   * calling {@link S2Polyline#decode}. The encoding is byte-compatible with the C++ version of the
   * S2 library.
   *
   * @param output The output stream into which the encoding should be written.
   * @throws IOException if there was a problem writing into the output stream.
   */
  @JsIgnore // OutputStream is not available to J2CL.
  public void encodeCompact(OutputStream output) throws IOException {
    int level = numVertices == 0 ? S2CellId.MAX_LEVEL : getBestSnapLevel();
    LittleEndianOutput encoder = new LittleEndianOutput(output);
    if (level == -1) {
      encodeUncompressed(encoder);
    } else {
      encodeCompressed(level, encoder);
    }
  }

  /** Encodes this polyline into the given little endian output stream. */
  void encodeUncompressed(LittleEndianOutput os) throws IOException {
    os.writeByte(LOSSLESS_ENCODING_VERSION);
    os.writeInt(numVertices);
    for (S2Point p : vertices) {
      p.encode(os);
    }
  }

  /** Encodes a compressed polyline at requested snap level. */
  void encodeCompressed(int snapLevel, LittleEndianOutput encoder) throws IOException {
    encoder.writeByte(COMPRESSED_ENCODING_VERSION);
    encoder.writeByte((byte) snapLevel);
    encoder.writeVarint32(numVertices);
    S2PointCompression.encodePointsCompressed(vertices(), snapLevel, encoder);
  }

  @JsIgnore // InputStream is not available to J2CL.
  public static S2Polyline decode(InputStream is) throws IOException {
    return decode(new LittleEndianInput(is));
  }

  static S2Polyline decode(LittleEndianInput decoder) throws IOException {
    byte version = decoder.readByte();
    switch (version) {
      case LOSSLESS_ENCODING_VERSION:
        return decodeLossless(decoder);
      case COMPRESSED_ENCODING_VERSION:
        return decodeCompressed(decoder);
      default:
        throw new IOException("Unsupported S2Polyline encoding version " + version);
    }
  }

  private static S2Polyline decodeLossless(LittleEndianInput is) throws IOException {
    S2Point[] vertices = new S2Point[is.readInt()];
    for (int i = 0; i < vertices.length; i++) {
      vertices[i] = S2Point.decode(is);
    }
    return new S2Polyline(vertices);
  }

  private static S2Polyline decodeCompressed(LittleEndianInput decoder) throws IOException {
    int level = decoder.readByte();
    if (level > S2CellId.MAX_LEVEL) {
      throw new IOException("Invalid level " + level);
    }
    int numVertices = decoder.readVarint32();
    return new S2Polyline(S2PointCompression.decodePointsCompressed(numVertices, level, decoder));
  }

  /**
   * If all of the polyline's vertices happen to be the centers of S2Cells at some level, then
   * returns that level, otherwise returns -1. See also {@link #fromSnapped(S2Polyline, int)}.
   * Returns -1 if the polyline has no vertices.
   */
  public int getSnapLevel() {
    int snapLevel = -1;
    for (S2Point p : vertices) {
      FaceSiTi faceSiTi = S2Projections.PROJ.xyzToFaceSiTi(p);
      int level = S2Projections.PROJ.levelIfCenter(faceSiTi, p);
      if (level < 0) {
        // Vertex is not a cell center.
        return level;
      }
      if (level != snapLevel) {
        if (snapLevel < 0) {
          // First vertex.
          snapLevel = level;
        } else {
          // Vertices at more than one cell level.
          return -1;
        }
      }
    }

    return snapLevel;
  }

  /**
   * Computes the level at which most of the vertices are snapped. If multiple levels have the same
   * maximum number of vertices snapped to it, the first one (lowest level number / largest area /
   * smallest encoding length) will be chosen, so this is desired. Returns -1 for unsnapped
   * polylines.
   */
  int getBestSnapLevel() {
    int[] histogram = new int[S2CellId.MAX_LEVEL + 1];
    for (S2Point p : vertices) {
      FaceSiTi faceSiTi = S2Projections.PROJ.xyzToFaceSiTi(p);
      int level = S2Projections.PROJ.levelIfCenter(faceSiTi, p);
      // Level is -1 for unsnapped points.
      if (level >= 0) {
        histogram[level]++;
      }
    }
    int snapLevel = 0;
    for (int i = 1; i < histogram.length; i++) {
      if (histogram[i] > histogram[snapLevel]) {
        snapLevel = i;
      }
    }
    if (histogram[snapLevel] == 0 && numVertices > 0) {
      // This is an unsnapped polyline.
      return -1;
    }
    return snapLevel;
  }

  /***
   * Returns a new S2Polyline with reversed order of vertices to the original.
   */
  public S2Polyline reversed() {
    return new S2Polyline(Lists.reverse(vertices()));
  }
}
