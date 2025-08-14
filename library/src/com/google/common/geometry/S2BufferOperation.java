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

import static com.google.common.geometry.S2.M_PI_2;
import static com.google.common.geometry.S2EdgeUtil.getPointOnRay;
import static com.google.common.geometry.S2Predicates.DBL_ERR;
import static com.google.common.geometry.S2Predicates.angleContainsVertex;
import static com.google.common.geometry.S2RobustCrossProd.robustCrossProd;
import static java.lang.Math.PI;
import static java.lang.Math.acos;
import static java.lang.Math.ceil;
import static java.lang.Math.cos;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.signum;
import static java.lang.Math.tan;

import com.google.common.base.Preconditions;
import com.google.common.collect.Lists;
import com.google.common.geometry.S2Builder.SnapFunction;
import com.google.common.geometry.S2BuilderSnapFunctions.IdentitySnapFunction;
import com.google.common.geometry.S2EdgeUtil.EdgeCrosser;
import java.util.ArrayList;
import java.util.List;

/**
 * This class provides a way to expand an arbitrary collection of geometry by a fixed radius (an
 * operation variously known as "buffering", "offsetting", or "Minkowski sum with a disc"). The
 * output consists of a polygon (possibly with multiple shells) that contains all points within the
 * given radius of the original geometry.
 *
 * <p>The radius can also be negative, in which case the geometry is contracted. This causes the
 * boundaries of polygons to shrink or disappear, and removes all points and polylines.
 *
 * <p>The input consists of a sequence of layers. Each layer may consist of any combination of
 * points, polylines, and polygons, with the restriction that polygon interiors within each layer
 * may not intersect any other geometry (including other polygon interiors). The output is the union
 * of the buffered input layers. Note that only a single layer is allowed if the buffer radius is
 * negative.
 *
 * <p>This class may be used to compute polygon unions by setting the buffer radius to zero. The
 * union is computed using a single snapping operation.
 *
 * <p>Note that if you only want to compute an S2CellId covering of the buffered geometry, it is
 * much faster to use S2ShapeIndexBufferedRegion instead.
 *
 * <p>Keywords: buffer, buffering, expand, expanding, offset, offsetting, widen, contract, shrink,
 * Minkowski sum
 */
@SuppressWarnings("Assertion")
public class S2BufferOperation {

  /**
   * For polylines, specifies whether the end caps should be round or flat. See {@link
   * Options#setEndCapStyle(EndCapStyle)} below.
   */
  public static enum EndCapStyle {
    ROUND,
    FLAT
  }

  /**
   * Specifies whether polylines should be buffered only on the left, only on the right, or on both
   * sides.
   */
  public static enum PolylineSide {
    LEFT,
    RIGHT,
    BOTH
  }

  /** Options for an S2BufferOperation. */
  public static class Options {
    public static final double MIN_ERROR_FRACTION = 1e-6;
    public static final double MAX_CIRCLE_SEGMENTS = 1570.7968503979573;

    private S1Angle bufferRadius = S1Angle.ZERO;
    private double errorFraction = 0.02;
    private EndCapStyle endCapStyle = EndCapStyle.ROUND;
    private PolylineSide polylineSide = PolylineSide.BOTH;
    private SnapFunction snapFunction;

    /** Constructor that sets default options. */
    public Options() {
      this.snapFunction = new IdentitySnapFunction(S1Angle.ZERO);
    }

    /** Convenience constructor that calls setBufferRadius(). */
    public Options(S1Angle bufferRadius) {
      this();
      this.bufferRadius = bufferRadius;
    }

    /** Options copy constructor. */
    public Options(Options other) {
      this.bufferRadius = other.bufferRadius;
      this.errorFraction = other.errorFraction;
      this.endCapStyle = other.endCapStyle;
      this.polylineSide = other.polylineSide;
      // TODO(torrey): Ensure snap functions are immutable, or make a copy.
      this.snapFunction = other.snapFunction;
    }

    /** Returns the buffer radius. See {@link #setBufferRadius(S1Angle)}. */
    public S1Angle bufferRadius() {
      return bufferRadius;
    }

    /**
     * If positive, specifies that all points within the given radius of the input geometry should
     * be added to the output. If negative, specifies that all points within the given radius of the
     * complement of the input geometry should be subtracted from the output. If the buffer radius
     * is zero then the input geometry is passed through to the output layer after first converting
     * points and polylines into degenerate loops.
     *
     * <p>DEFAULT: S1Angle.ZERO
     */
    public void setBufferRadius(S1Angle bufferRadius) {
      this.bufferRadius = bufferRadius;
    }

    /**
     * Returns the error fraction, which is the allowable error when buffering as a fraction of the
     * bufferRadius. See {@link #setErrorFraction(double)}.
     */
    public double errorFraction() {
      return errorFraction;
    }

    /**
     * Specifies the allowable error when buffering, expressed as a fraction of bufferRadius(). The
     * actual buffer distance will be in the range {@code [(1-f) * r - C, (1 + f) * r + C] } where
     * "f" is the error fraction, "r" is the buffer radius, and "C" is S2BufferOperation.kAbsError.
     *
     * <p>Be aware that the number of output edges increases proportionally to {@code
     * (1 / sqrt(errorFraction))}, so setting a small value may increase the size of the output
     * considerably.
     *
     * <pre>
     * REQUIRES: {@code errorFraction() >= MIN_ERROR_FRACTION}
     * REQUIRES: {@code errorFraction() <= 1.0}
     * DEFAULT: 0.01  (i.e., maximum error of 1%)
     * </pre>
     */
    public void setErrorFraction(double errorFraction) {
      Preconditions.checkArgument(errorFraction >= MIN_ERROR_FRACTION);
      Preconditions.checkArgument(errorFraction <= 1.0);
      this.errorFraction = max(MIN_ERROR_FRACTION, min(1.0, errorFraction));
    }

    /**
     * Returns the maximum error in the buffered result for the current bufferRadius(),
     * errorFraction(), and snapFunction(). Note that the error due to buffering consists of both
     * relative errors (those proportional to the buffer radius) and absolute errors. The maximum
     * relative error is controlled by errorFraction(), while the maximum absolute error is about 10
     * nanometers on the Earth's surface and is defined internally. The error due to snapping is
     * defined by the specified snapFunction().
     */
    public S1Angle maxError() {
      // See comments for MIN_REQUESTED_ERROR above.
      S2Builder.Builder builderOptions = new S2Builder.Builder(snapFunction);
      builderOptions.setSplitCrossingEdges(true);
      return S1Angle.max(MIN_REQUESTED_ERROR, bufferRadius.abs().mul(errorFraction))
          .add(MAX_ABSOLUTE_INTERPOLATION_ERROR)
          .add(builderOptions.maxEdgeDeviation());
    }

    /**
     * Returns the number of polyline segments used to approximate a planar circle. See
     * {@link #setCircleSegments(double)}.
     */
    public double circleSegments() {
      // TODO(ericv):
      //   This formula assumes that vertices can be placed anywhere.
      //   return PI / acos((1 - errorFraction) / (1 + errorFraction));

      // This formula assumes that all vertices are placed on the midline.
      return PI / acos(1 - errorFraction);
    }

    /**
     * As an alternative to {@link #setErrorFraction(double)}, errorFraction() may be specified as
     * the number of polyline segments used to approximate a planar circle. These two values are
     * related according to the formula
     *
     * <pre>
     *    errorFraction = (1 - cos(theta)) / (1 + cos(theta))
     *                  ~= 0.25 * (theta ** 2)
     * </pre>
     *
     * where (theta == Pi / circleSegments), i.e. error decreases quadratically with the number of
     * circle segments.
     *
     * <pre>
     * REQUIRES: {@code circleSegments() >= 2.0}
     * REQUIRES: {@code circleSegments() <= MAX_CIRCLE_SEGMENTS} (about 1570; corresponds to
     * MIN_ERROR_FRACTION)
     * DEFAULT: about 15.76 (corresponding to errorFraction() default value)
     * </pre>
     */
    public void setCircleSegments(double circleSegments) {
      Preconditions.checkArgument(circleSegments >= 2.0);
      Preconditions.checkArgument(circleSegments <= MAX_CIRCLE_SEGMENTS);

      // We convert circleSegments to errorFraction using planar geometry, because the number of
      // segments required to approximate a circle on the sphere to within a given tolerance is not
      // constant. Unlike in the plane, the total curvature of a circle on the sphere decreases as
      // the area enclosed by the circle increases; great circles have no curvature at all. We round
      // up when converting to ensure that we won't generate any tiny extra edges.

      // Note that we take advantage of both positive and negative errors when approximating circles
      // (i.e., vertices are not necessarily on the midline) and thus the relationships between
      // circleSegments and errorFraction are
      //        e = (1 - cos(Pi/n)) / (1 + cos(Pi/n))
      //        n = Pi / acos((1 - e) / (1 + e))
      // double r = cos(PI / circleSegments);
      // setErrorFraction((1 - r) / (1 + r) + 1e-15);

      // When all vertices are on the midline, the relationships are
      //        e = 1 - cos(Pi/n)
      //        n = Pi / acos(1 - e)
      setErrorFraction(1.0 - cos(PI / circleSegments) + 1e-15);
    }

    /** Returns the end cap style. See {@link #setEndCapStyle(EndCapStyle)}. */
    public EndCapStyle endCapStyle() {
      return endCapStyle;
    }

    /**
     * For polylines, specifies whether the end caps should be round or flat. Note that with flat
     * end caps, there is no buffering beyond the polyline endpoints (unlike "square" end caps,
     * which are not implemented).
     *
     * <p>DEFAULT: EndCapStyle.ROUND
     */
    public void setEndCapStyle(EndCapStyle endCapStyle) {
      this.endCapStyle = endCapStyle;
    }

    /**
     * Returns the polyline side, indicating if polylines should be buffered on the right, the
     * left, or both. See {@link #setPolylineSide(PolylineSide)}.
     */
    public PolylineSide polylineSide() {
      return polylineSide;
    }

    /**
     * Specifies whether polylines should be buffered only on the left, only on the right, or on
     * both sides. For one-sided buffering please note the following:
     *
     * <ul>
     *   <li>EndCapStyle.ROUND yields two quarter-circles, one at each end.
     *   <li>To buffer by a different radius on each side of the polyline, you can use two
     *       S2BufferOperations and compute their union. (Note that round end caps will yield two
     *       quarter-circles at each end of the polyline with different radii.)
     *   <li>Polylines consisting of a single degenerate edge are always buffered identically to
     *       points, i.e. this option has no effect.
     *   <li>When the polyline turns right by more than 90 degrees, buffering may or may not extend
     *       to the non-buffered side of the polyline. For example if ABC makes a 170 degree right
     *       turn at B, it is unspecified whether the buffering of AB extends across edge BC and
     *       vice versa. Similarly if ABCD represents two right turns of 90 degrees where AB and CD
     *       are separated by less than the buffer radius, it is unspecified whether buffering of AB
     *       extends across CD and vice versa.
     * </ul>
     *
     * DEFAULT: PolylineSide.BOTH
     */
    public void setPolylineSide(PolylineSide polylineSide) {
      this.polylineSide = polylineSide;
    }

    /** Returns the snap function. See {@link #setSnapFunction(SnapFunction)}. */
    public SnapFunction snapFunction() {
      return snapFunction;
    }

    /**
     * Specifies the function used for snap rounding the output during the call to build(). Note
     * that any errors due to snapping are in addition to those specified by errorFraction().
     *
     * <p>DEFAULT: s2builderutil.IdentitySnapFunction(S1Angle.ZERO)
     */
    public void setSnapFunction(SnapFunction snapFunction) {
      this.snapFunction = snapFunction;
    }
  }

  // The errors due to buffering can be categorized as follows:
  //
  //  1. Requested error. This represents the error due to approximating the buffered boundary as a
  //     sequence of line segments rather than a sequence of circular arcs. It is largely controlled
  //     by options.errorFraction(), and can be bounded as
  //
  //       max(MIN_REQUESTED_ERROR, errorFraction * bufferRadius)
  //
  //     where MIN_REQUESTED_ERROR reflects the fact that S2Points do not have infinite precision.
  //     (For example, it makes no sense to attempt to buffer geometry by 1e-100 radians because the
  //     spacing between representable S2Points is only about 2e-16 radians in general.)
  //
  //  2. Relative interpolation errors. These are numerical errors whose magnitude is proportional
  //     to the buffer radius. For such errors the worst-case coefficient of proportionality turns
  //     out to be so tiny compared to the smallest allowable error fraction (MIN_ERROR_FRACTION)
  //     that we simply ignore such errors.
  //
  //  3. Absolute interpolation errors. These are numerical errors that are not proportional to the
  //     buffer radius. The main sources of such errors are (1) calling robustCrossProd() to compute
  //     edge normals, and (2) calls to getPointOnRay() to interpolate points along the buffered
  //     boundary. It is possible to show that this error is at most
  //     MAX_ABSOLUTE_INTERPOLATION_ERROR as defined below.
  //
  // Putting these together, the final error bound looks like this:
  //
  //   maxError = MAX_ABSOLUTE_INTERPOLATION_ERROR +
  //               max(MIN_REQUESTED_ERROR,
  //                   max(MIN_ERROR_FRACTION, options.errorFraction()) *
  //                   options.bufferRadius())

  /**
   * The maximum angular spacing between representable S2Points on the unit sphere is roughly 2 *
   * DBL_ERR. We require the requested absolute error to be at least this large because attempting
   * to achieve a smaller error does not increase the precision of the result and can increase the
   * running time and space requirements considerably.
   */
  private static final S1Angle MIN_REQUESTED_ERROR = S1Angle.radians(2 * DBL_ERR);

  /**
   * The maximum absolute error due to interpolating points on the buffered boundary. The following
   * constant bounds the maximum additional error perpendicular to the buffered boundary due to all
   * steps of the calculation (S2.RobustCrossProd, the two calls to GetPointOnRay, etc).
   *
   * <p>This distance represents about 10 nanometers on the Earth's surface. Note that this is a
   * conservative upper bound and that it is difficult to construct inputs where the error is
   * anywhere close to this large.
   */
  private static final S1Angle MAX_ABSOLUTE_INTERPOLATION_ERROR =
      S2EdgeUtil.GET_POINT_ON_LINE_ERROR.add(S2EdgeUtil.GET_POINT_ON_RAY_PERPENDICULAR_ERROR);

  /** The current Options for this S2BufferOperation. Set by the call to init(). */
  private Options options;

  /**
   * The number of layers containing two-dimension geometry that have been added so far. This is
   * used to enforce the requirement that negative buffer radii allow only a single such layer.
   */
  private int numPolygonLayers = 0;

  // Parameters for buffering vertices and edges.

  /** The sign of bufferRadius (-1, 0, or +1). */
  private int bufferSign;
  private S1ChordAngle absRadius;
  private S1ChordAngle vertexStep = S1ChordAngle.ZERO;
  private S1ChordAngle edgeStep = S1ChordAngle.ZERO;

  /**
   * We go to extra effort to ensure that points are transformed into regular polygons. (We don't do
   * this for arcs in general because we would rather use the allowable error to reduce the
   * complexity of the output rather than increase its symmetry.)
   */
  private S1ChordAngle pointStep; // set by init()

  /** Contains the buffered loops that have been accumulated so far. */
  private final S2WindingOperation op = new S2WindingOperation();

  /**
   * The current offset path. When each path is completed into a loop it is added to "op" (the
   * S2WindingOperation).
   */
  private final List<S2Point> path = new ArrayList<>();

  /**
   * As buffered loops are added we keep track of the winding number of a fixed reference point.
   * This is used to derive the winding numbers of every region in the spherical partition induced
   * by the buffered loops.
   */
  private S2Point refPoint; // set by init()

  /** The winding number associated with refPoint. */
  private int refWinding; // set by init()

  /** An endpoint of the current sweep edge. sweepA is a vertex of the original geometry. */
  private S2Point sweepA;
  /** An endpoint of the current sweep edge. sweepB is a vertex of the current offset path. */
  private S2Point sweepB;

  /**
   * The starting vertex of the current input loop. With offsetStart, used to close the buffer
   * region when a loop is completed.
   */
  private S2Point inputStart;
  /**
   * The starting vertex of the current offset curve. With inputSTart, used to close the buffer
   * region when a loop is completed.
   */
  private S2Point offsetStart;
  /** Is inputStart set? */
  private boolean haveInputStart; // set by init()
  /** Is offsetStart set? */
  private boolean haveOffsetStart; // set by init()

  /** Constructor that does not set options. Requires init() to be called. */
  public S2BufferOperation() {}

  /** Convenience constructor that calls init(). */
  public S2BufferOperation(S2BuilderLayer resultLayer, Options options) {
    init(resultLayer, options);
  }

  /**
   * Starts a buffer operation that sends the output polygon to the given S2Builder layer. This
   * method may be called more than once.
   *
   * <p>Note that buffering always yields a polygon, even if the input includes polylines and
   * points. If the buffer radius is zero, points and polylines will be converted into degenerate
   * polygon loops; if the buffer radius is negative, points and polylines will be removed.
   */
  public void init(S2BuilderLayer resultLayer, Options options) {
    this.options = options;
    this.refPoint = S2.origin();
    this.refWinding = 0;
    this.haveInputStart = false;
    this.haveOffsetStart = false;
    this.bufferSign = (int) signum(options.bufferRadius().radians());

    S1Angle absRadius = options.bufferRadius().abs();
    S1Angle requestedError = S1Angle.max(MIN_REQUESTED_ERROR, absRadius.mul(options.errorFraction()));
    S1Angle maxError = MAX_ABSOLUTE_INTERPOLATION_ERROR.add(requestedError);

    if (absRadius.lessOrEquals(maxError)) {
      // If the requested radius is smaller than the maximum error, buffering could yield points on
      // the wrong side of the original input boundary (e.g., shrinking geometry slightly rather
      // than expanding it). Rather than taking that risk, we set the buffer radius to zero when
      // this happens (which causes the original geometry to be returned).
      this.absRadius = S1ChordAngle.ZERO;
      this.bufferSign = 0;
    } else if (absRadius.add(maxError).greaterOrEquals(S1Angle.radians(PI))) {
      // If the permissible range of buffer angles includes Pi then we might as well take advantage
      // of that.
      this.absRadius = S1ChordAngle.STRAIGHT;
    } else {
      this.absRadius = S1ChordAngle.fromS1Angle(absRadius);
      this.vertexStep = S1ChordAngle.fromS1Angle(getMaxEdgeSpan(absRadius, requestedError));

      // We take extra care to ensure that points are buffered as regular polygons. The step angle
      // is adjusted up slightly to ensure that we don't wind up with a tiny extra edge.
      this.pointStep = S1ChordAngle.fromRadians(2 * PI / ceil(2 * PI / vertexStep.radians()) + 1e-15);

      // Edges are buffered only if the buffer radius (including permissible error) is less than 90
      // degrees.
      S1Angle edgeRadius = S1Angle.radians(M_PI_2).sub(absRadius);
      if (edgeRadius.greaterThan(maxError)) {
        this.edgeStep = S1ChordAngle.fromS1Angle(getMaxEdgeSpan(edgeRadius, requestedError));
      }
    }

    // The buffered output should include degeneracies (i.e., isolated points and/or sibling edge
    // pairs) only if (1) the user specified a non-negative buffer radius, and (2) the adjusted
    // buffer radius is zero. The only purpose of keeping degeneracies is to allow points/polylines
    // in the input geometry to be converted back to points/polylines in the output if the client so
    // desires.
    S2WindingOperation.Options windingOptions =
        new S2WindingOperation.Options(options.snapFunction());
    windingOptions.setIncludeDegeneracies(
        bufferSign == 0 && options.bufferRadius().greaterOrEquals(S1Angle.ZERO));
    op.init(resultLayer, windingOptions);
  }

  /** Returns the current Options of this S2BufferOperation. */
  public Options options() {
    return options;
  }

  // Each call below represents a different input layer. Note that if the buffer radius is negative,
  // then at most one input layer is allowed (ignoring any layers that contain only points and
  // polylines).

  /** Adds an input layer containing a single point. */
  public void addPoint(S2Point point) {
    // If bufferRadius < 0, points are discarded.
    if (bufferSign < 0) {
      return;
    }

    // Buffering by 180 degrees or more always yields the full polygon. (We don't need to worry
    // about buffering by 180 degrees yielding a degenerate hole because errorFraction is always
    // positive.
    if (absRadius.greaterOrEquals(S1ChordAngle.STRAIGHT)) {
      addFullPolygon();
      return;
    }

    // If bufferRadius == 0, points are converted into degenerate loops.
    if (bufferSign == 0) {
      path.add(point);
    } else {
      // Since S1ChordAngle can only represent angles between 0 and 180 degrees, we generate the
      // circle in four 90 degree increments.
      setInputVertex(point);
      S2Point start = S2.ortho(point);
      S1ChordAngle angle = S1ChordAngle.ZERO;
      for (int quadrant = 0; quadrant < 4; ++quadrant) {
        // Generate 90 degrees of the circular arc. Normalize "rotateDir" at each iteration to avoid
        // magnifying normalization errors in "point".
        S2Point rotateDir = point.crossProd(start).normalize();
        for (; angle.lessThan(S1ChordAngle.RIGHT); angle = S1ChordAngle.add(angle, pointStep)) {
          S2Point dir = getPointOnRay(start, rotateDir, angle);
          addOffsetVertex(getPointOnRay(point, dir, absRadius));
        }
        angle = S1ChordAngle.sub(angle, S1ChordAngle.RIGHT);
        start = rotateDir;
      }
      closeBufferRegion();
    }
    outputPath();
  }

  /**
   * Adds an input layer containing a polyline. Note the following:
   *
   * <ul>
   *   <li>Polylines with 0 or 1 vertices are considered to be empty.
   *   <li>A polyline with 2 identical vertices is equivalent to a point.
   *   <li>Polylines have end caps. See {@link Options#endCapStyle()}.
   *   <li>One-sided polyline buffering is supported. See {@link Options#polylineSide()}.
   * </ul>
   */
  public void addPolyline(List<S2Point> polyline) {
    // Left-sided buffering is supported by reversing the polyline and then buffering on the right.
    if (options.polylineSide() == PolylineSide.LEFT) {
      List<S2Point> reversed = Lists.reverse(polyline);
      polyline = reversed;
    }

    // If bufferRadius < 0, polylines are discarded.
    if (bufferSign < 0) {
      return;
    }

    // Polylines with 0 or 1 vertices are defined to have no edges.
    int n = polyline.size();
    if (n <= 1) {
      return;
    }
    // At least two vertices, one edge.

    // Polylines with one degenerate edge are treated as points.
    if (n == 2 && polyline.get(0).equalsPoint(polyline.get(1))) {
      addPoint(polyline.get(0));
      return;
    }

    // Buffering by 180 degrees or more always yields the full polygon.
    if (absRadius.greaterOrEquals(S1ChordAngle.STRAIGHT)) {
      addFullPolygon();
      return;
    }

    // If bufferRadius == 0, polylines are converted into degenerate loops.
    if (bufferSign == 0) {
      path.addAll(polyline.subList(0, polyline.size() - 2));
      path.addAll(Lists.reverse(polyline).subList(0, polyline.size() - 2));
    } else {
      // Otherwise we buffer each side of the polyline separately.
      setInputVertex(polyline.get(0));
      addStartCap(polyline.get(0), polyline.get(1));
      for (int i = 0; i < n - 2; ++i) {
        bufferEdgeAndVertex(polyline.get(i), polyline.get(i + 1), polyline.get(i + 2));
      }
      addEdgeArc(polyline.get(n - 2), polyline.get(n - 1));
      addEndCap(polyline.get(n - 2), polyline.get(n - 1));

      if (options.polylineSide() == PolylineSide.BOTH) {
        for (int i = n - 3; i >= 0; --i) {
          bufferEdgeAndVertex(polyline.get(i + 2), polyline.get(i + 1), polyline.get(i));
        }
        addEdgeArc(polyline.get(1), polyline.get(0));
        closeBufferRegion();
      } else {
        // The other side of the polyline is not buffered. Note that for PolylineSide.LEFT, the
        // polyline direction has been reversed.
        path.addAll(Lists.reverse(polyline));
        // Don't call closeBufferRegion() since the path has already been closed.
      }
    }
    outputPath();
  }

  /**
   * Adds an input layer containing a loop. Note the following:
   *
   * <ul>
   *   <li>A loop with no vertices is empty.
   *   <li>A loop with 1 vertex is equivalent to a point.
   *   <li>The interior of the loop is on its left.
   *   <li>Buffering a self-intersecting loop produces undefined results.
   * </ul>
   */
  public void addLoop(List<S2Point> loop) {
    if (loop.isEmpty()) {
      return;
    }
    // TODO(torrey): What about one and two vertex loops? S2Loop must have three.
    bufferLoop(loop);

    refWinding += S2ShapeUtil.containsBruteForce(new S2Loop(loop), refPoint) ? 1 : 0;
    numPolygonLayers += 1;
  }

  /**
   * Adds an input layer containing the given shape. Shapes are handled as points, polylines, or
   * polygons according to the rules above. In addition note the following:
   *
   * <ul>
   *   <li>Polygon holes may be degenerate (e.g., consisting of a single vertex or entirely of
   *       sibling pairs such as ABCBDB).
   *   <li>Full polygons are supported. Note that since full polygons do not have a boundary, they
   *       are not affected by buffering.
   */
  public void addShape(S2Shape shape) {
    bufferShape(shape);
    refWinding += S2ShapeUtil.containsBruteForce(shape, refPoint) ? 1 : 0;
    numPolygonLayers += (shape.dimension() == 2) ? 1 : 0;
  }

  /**
   * Adds an input layer containing all of the shapes in the given index.
   *
   * <p>REQUIRES: The interiors of polygons must be disjoint from all other indexed geometry,
   * including other polygon interiors. (S2BooleanOperation also requires this.)
   */
  public void addShapeIndex(S2ShapeIndex index) {
    int maxDimension = -1;
    for (S2Shape shape : index.getShapes()) {
      maxDimension = max(maxDimension, shape.dimension());
      bufferShape(shape);
    }
    refWinding += new S2ContainsPointQuery(index).contains(refPoint) ? 1 : 0;
    numPolygonLayers += maxDimension == 2 ? 1 : 0;
  }

  /**
   * Computes the union of the buffered input shapes and sends the output polygon to the S2Builder
   * layer specified in the constructor. Returns true on success and otherwise sets "error"
   * appropriately.
   *
   * <p>Note that if the buffer radius is negative, only a single input layer is allowed (ignoring
   * any layers that contain only points and polylines).
   *
   * <p>The algorithm below essentially computes the offset curve of the original boundary, and uses
   * this curve to divide the sphere into regions of constant winding number. Since winding numbers
   * on the sphere are relative rather than absolute (see s2winding_operation.h), we also need to
   * keep track of the desired winding number at a fixed reference point. The initial winding number
   * for this point is the number of input shapes that contain it. We then update it during the
   * buffering process by imagining a "sweep edge" that extends from the current point A on the
   * input boundary to the corresponding point B on the offset curve. As we process an input loop
   * and generate the corresponding offset curve, the sweep edge moves continuously and covers the
   * entire buffer region (i.e., the region added to or subtracted from the input geometry). We
   * increase the winding number of the reference point by one whenever it crosses the sweep edge
   * from left to right, and we decrease the winding number by one whenever it crosses the sweep
   * edge from right to left.
   *
   * <p>Concave vertices require special handling, because the corresponding offset curve can leave
   * behind regions whose winding number is zero or negative. We handle this by splicing the concave
   * vertex into the offset curve itself; this effectively terminates the current buffer region and
   * starts a new one, such that the region of overlap is counted twice (i.e., its winding number
   * increases by two). The end result is the same as though we had computed the union of a sequence
   * of buffered convex boundary segments. This trick is described in the following paper: "Polygon
   * Offsetting by Computing Winding Numbers" (Chen and McMains, Proceedings of IDETC/CIE 2005).
   *
   * <p>TODO(ericv): The algorithm below is much faster than, say, computing the union of many
   * buffered edges. However further improvements are possible. In particular, there is an
   * unimplemented optimization that would make it much faster to buffer concave boundaries when the
   * buffer radius is large.
   */
  public boolean build(S2Error error) {
    if (bufferSign < 0 && numPolygonLayers > 1) {
      error.init(
          S2Error.Code.FAILED_PRECONDITION,
          "Negative buffer radius requires at most one polygon layer");
      return false;
    }
    return op.build(refPoint, refWinding, S2WindingOperation.WindingRule.POSITIVE, error);
  }

  private S1Angle getMaxEdgeSpan(S1Angle radius, S1Angle requestedError) {
    // If the allowable radius range spans Pi/2 then we can use edges as long as we like, however we
    // always use at least 3 edges to approximate a circle.
    // TODO(torrey): Just use doubles for step and minRadius
    S1Angle step = S1Angle.radians(2 * PI / 3 + 1e-15);
    S1Angle minRadius = radius.sub(requestedError);
    assert minRadius.greaterOrEquals(S1Angle.ZERO);
    if (radius.radians() < M_PI_2) {
      step =
          S1Angle.min(
              step, S1Angle.radians(2 * acos(tan(minRadius.radians()) / tan(radius.radians()))));
    } else if (minRadius.radians() > M_PI_2) {
      step =
          S1Angle.min(
              step, S1Angle.radians(2 * acos(tan(radius.radians()) / tan(minRadius.radians()))));
    }
    return step;
  }

  private void setInputVertex(S2Point newA) {
    if (haveInputStart) {
      assert haveOffsetStart;
      updateRefWinding(sweepA, sweepB, newA);
    } else {
      inputStart = newA;
      haveInputStart = true;
    }
    sweepA = newA;
  }

  /**
   * Adds the point "newB" to the offset path. Also advances the sweep edge AB // by moving its
   * second vertex B to "newB" and updating the winding number of // the reference point if
   * necessary (see introduction).
   */
  private void addOffsetVertex(S2Point newB) {
    path.add(newB);
    if (haveOffsetStart) {
      assert haveInputStart;
      updateRefWinding(sweepA, sweepB, newB);
    } else {
      offsetStart = newB;
      haveOffsetStart = true;
    }
    sweepB = newB;
  }

  /**
   * Finishes buffering the current loop by advancing the sweep edge back to its starting location,
   * updating the winding number of the reference point if necessary.
   */
  private void closeBufferRegion() {
    if (haveOffsetStart && haveInputStart) {
      updateRefWinding(sweepA, sweepB, inputStart);
      updateRefWinding(inputStart, sweepB, offsetStart);
    }
  }

  /**
   * Outputs the current buffered path (which is assumed to be a loop), and resets the state to
   * prepare for buffering a new loop.
   */
  private void outputPath() {
    op.addLoop(path);
    path.clear(); // Does not change capacity.
    haveInputStart = false;
    haveOffsetStart = false;
  }

  /**
   * Given a triangle ABC that has just been covered by the sweep edge AB, updates the winding
   * number of the reference point if necessary.
   */
  private void updateRefWinding(S2Point a, S2Point b, S2Point c) {
    // TODO(ericv): This code could be made much faster by maintaining a bounding plane that
    // separates the current sweep edge from the reference point. Whenever the sweepA or sweepB is
    // updated we would just need to check that the new vertex is still on the opposite side of the
    // bounding plane (i.e., one dot product). If not, we test the current triangle using the code
    // below and then compute a new bounding plane.
    //
    // Another optimization would be to choose the reference point to be 90 degrees away from the
    // first input vertex, since then triangle tests would not be needed unless the input geometry
    // spans more than 90 degrees. This would involve adding a new flag have_refPoint rather than
    // always choosing the reference point to be S2.Origin().
    //
    // According to profiling these optimizations are not currently worthwhile, but this is worth
    // revisiting if and when other improvements are made.
    int sign = S2Predicates.sign(a, b, c);
    if (sign == 0) {
      return;
    }
    boolean inside = angleContainsVertex(a, b, c) == (sign > 0);
    EdgeCrosser crosser = new EdgeCrosser(b, refPoint);
    inside ^= crosser.edgeOrVertexCrossing(a, b);
    inside ^= crosser.edgeOrVertexCrossing(b, c);
    inside ^= crosser.edgeOrVertexCrossing(c, a);
    if (inside) {
      refWinding += sign;
    }
  }

  /** Ensures that the output will be the full polygon. */
  private void addFullPolygon() {
    refWinding += 1;
  }

  /**
   * Returns the edge normal for the given edge AB. The sign is chosen such that the normal is on
   * the right of AB if bufferSign > 0, and on the left of AB if bufferSign < 0.
   */
  private S2Point getEdgeAxis(S2Point a, S2Point b) {
    assert bufferSign != 0;
    return robustCrossProd(b, a).normalize().mul(bufferSign);
  }

  /**
   * Adds a semi-open offset arc around vertex V. The arc proceeds CCW from "start" to "end" (both
   * of which must be perpendicular to V).
   */
  private void addVertexArc(S2Point v, S2Point start, S2Point end) {
    // Make sure that we output at least one point even when span == 0.
    S2Point rotateDir = v.crossProd(start).normalize().mul(bufferSign);
    S1ChordAngle angle = S1ChordAngle.ZERO;
    S1ChordAngle span = new S1ChordAngle(start, end);
    do {
      S2Point dir = getPointOnRay(start, rotateDir, angle);
      addOffsetVertex(getPointOnRay(v, dir, absRadius));
      angle = S1ChordAngle.add(angle, vertexStep);
    } while (angle.lessThan(span));
  }

  /** Closes the semi-open arc generated by addVertexArc(). */
  private void closeVertexArc(S2Point v, S2Point end) {
    addOffsetVertex(getPointOnRay(v, end, absRadius));
  }

  /** Adds a semi-open offset arc for the given edge AB. */
  private void addEdgeArc(S2Point a, S2Point b) {
    S2Point abAxis = getEdgeAxis(a, b);
    if (edgeStep.isZero()) {
      // If the buffer radius is more than 90 degrees, edges do not contribute to the buffered
      // boundary. Instead we force the offset path to pass through a vertex located at the edge
      // normal. This is similar to the case of concave vertices (below) where it is necessary to
      // route the offset path through the concave vertex to ensure that the winding numbers in all
      // output regions have the correct sign.
      addOffsetVertex(abAxis);
    } else {
      // Make sure that we output at least one point even when span == 0.
      S2Point rotateDir = a.crossProd(abAxis).normalize().mul(bufferSign);
      S1ChordAngle angle = S1ChordAngle.ZERO;
      S1ChordAngle span = new S1ChordAngle(a, b);
      do {
        S2Point p = getPointOnRay(a, rotateDir, angle);
        addOffsetVertex(getPointOnRay(p, abAxis, absRadius));
        angle = S1ChordAngle.add(angle, edgeStep);
      } while (angle.lessThan(span));
    }
    setInputVertex(b);
  }

  /** Closes the semi-open arc generated by addEdgeArc(). */
  private void closeEdgeArc(S2Point a, S2Point b) {
    if (!edgeStep.isZero()) {
      addOffsetVertex(getPointOnRay(b, getEdgeAxis(a, b), absRadius));
    }
  }

  /**
   * Buffers the edge AB and the vertex B. (The vertex C is used to determine the range of angles
   * that should be buffered at B.)
   *
   * <p>TODO(ericv): Let A* denote the possible offset points of A with respect to the edge AB for
   * buffer radii in the range specified by "radius" and "errorFraction". Rather than requiring that
   * the path so far terminates at a point in A*, as you might expect, instead we only require that
   * the path terminates at a point X such that for any point Y in A*, the edge XY does not leave
   * the valid buffer zone of the previous edge and vertex.
   */
  private void bufferEdgeAndVertex(S2Point a, S2Point b, S2Point c) {
    assert !a.equalsPoint(b);
    assert !b.equalsPoint(c);
    assert bufferSign != 0;

    // For left (convex) turns we need to add an offset arc. For right (concave) turns we connect
    // the end of the current offset path to the vertex itself and then to the start of the offset
    // path for the next edge. Note that A == C is considered to represent a convex (left) turn.
    addEdgeArc(a, b);
    if (bufferSign * S2Predicates.sign(a, b, c) >= 0) {
      // The boundary makes a convex turn. If there is no following edge arc then we need to
      // generate a closed vertex arc.
      S2Point start = getEdgeAxis(a, b);
      S2Point end = getEdgeAxis(b, c);
      addVertexArc(b, start, end);
      if (edgeStep.isZero()) {
        closeVertexArc(b, end);
      }
    } else {
      // The boundary makes a concave turn. It is tempting to simply connect the end of the current
      // offset path to the start of the offset path for the next edge, however this can create
      // output regions where the winding number is incorrect. A solution that always works is to
      // terminate the current offset path and start a new one by connecting the two offset paths
      // through the input vertex whenever it is concave. We first need to close the previous
      // semi-open edge arc if necessary.
      closeEdgeArc(a, b);
      addOffsetVertex(b); // Connect through the input vertex.
    }
  }

  /**
   * Given a polyline that starts with the edge AB, adds an end cap (as specified by endCapStyle()
   * and polylineSide()) for the vertex A.
   */
  private void addStartCap(S2Point a, S2Point b) {
    S2Point axis = getEdgeAxis(a, b);
    if (options.endCapStyle() == EndCapStyle.FLAT) {
      // One-sided flat end caps require no additional vertices since the "offset curve" for the
      // opposite side is simply the reversed polyline.
      if (options.polylineSide() == PolylineSide.BOTH) {
        addOffsetVertex(getPointOnRay(a, axis.neg(), absRadius));
      }
    } else {
      assert options.endCapStyle() == EndCapStyle.ROUND;
      if (options.polylineSide() == PolylineSide.BOTH) {
        // The end cap consists of a semicircle.
        addVertexArc(a, axis.neg(), axis);
      } else {
        // The end cap consists of a quarter circle. Note that for PolylineSide.LEFT, the polyline
        // direction has been reversed.
        addVertexArc(a, axis.crossProd(a).normalize(), axis);
      }
    }
  }

  /**
   * Given a polyline that ends with the edge AB, adds an end cap (as specified by endCapStyle() and
   * polylineSide()) for the vertex B.
   */
  private void addEndCap(S2Point a, S2Point b) {
    S2Point axis = getEdgeAxis(a, b);
    if (options.endCapStyle() == EndCapStyle.FLAT) {
      closeEdgeArc(a, b); // Close the previous semi-open edge arc if necessary.
    } else {
      assert options.endCapStyle() == EndCapStyle.ROUND;
      if (options.polylineSide() == PolylineSide.BOTH) {
        // The end cap consists of a semicircle.
        addVertexArc(b, axis, axis.neg());
      } else {
        // The end cap consists of a quarter circle. We close the arc since it will be followed by
        // the reversed polyline vertices. Note that for PolylineSide.LEFT, the polyline direction
        // has been reversed.
        S2Point end = b.crossProd(axis).normalize();
        addVertexArc(b, axis, end);
        closeVertexArc(b, end);
      }
    }
  }

  /** Helper function that buffers the given loop. */
  private void bufferLoop(List<S2Point> loop) {
    // Empty loops always yield an empty path.
    if (loop.isEmpty()) {
      return;
    }

    // Loops with one degenerate edge are treated as points.
    // TODO(torrey): What about a loop with two identical vertices?
    if (loop.size() == 1) {
      addPoint(loop.get(0));
      return;
    }

    // Buffering by 180 degrees or more always yields the full polygon.
    // Buffering by -180 degrees or more always yields the empty polygon.
    if (absRadius.greaterOrEquals(S1ChordAngle.STRAIGHT)) {
      if (bufferSign > 0) {
        addFullPolygon();
      }
      return;
    }

    // If bufferRadius == 0, the loop is passed through unchanged.
    if (bufferSign == 0) {
      path.addAll(loop);
    } else {
      setInputVertex(loop.get(0));
      // Allow indexing past or around the end of the loop.
      List<S2Point> loopSpan = S2ShapeUtil.s2PointLoopList(loop);
      for (int i = 0; i < loop.size(); ++i) {
        bufferEdgeAndVertex(loopSpan.get(i), loopSpan.get(i + 1), loopSpan.get(i + 2));
      }
      closeBufferRegion();
    }
    outputPath();
  }

  private void bufferShape(S2Shape shape) {
    int dimension = shape.dimension();
    int numChains = shape.numChains();
    for (int c = 0; c < numChains; ++c) {
      if (shape.getChainLength(c) == 0) {
        continue;
      }
      if (dimension == 0) {
        addPoint(shape.getChainVertex(c, 0));
      } else {
        if (dimension == 1) {
          addPolyline(shape.chain(c));
        } else {
          bufferLoop(shape.chain(c));
        }
      }
    }
  }
}
