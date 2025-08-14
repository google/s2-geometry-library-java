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

import static com.google.common.collect.Iterables.concat;
import static com.google.common.collect.Iterables.transform;
import static com.google.common.geometry.S2.skipAssertions;
import static com.google.common.geometry.S2Projections.MIN_WIDTH;
import static java.lang.Math.PI;
import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;
import static java.lang.Math.sqrt;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.base.Preconditions;
import com.google.common.base.Predicate;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.google.common.collect.TreeMultimap;
import com.google.common.geometry.PrimitiveArrays.Bytes;
import com.google.common.geometry.PrimitiveArrays.Cursor;
import com.google.common.geometry.S2Builder.SnapFunction;
import com.google.common.geometry.S2BuilderGraph.PolylineType;
import com.google.common.geometry.S2BuilderSnapFunctions.IdentitySnapFunction;
import com.google.common.geometry.S2BuilderSnapFunctions.S2CellIdSnapFunction;
import com.google.common.geometry.S2CrossingEdgeQuery.Edges;
import com.google.common.geometry.S2Projections.FaceSiTi;
import com.google.common.geometry.S2Shape.MutableEdge;
import com.google.common.geometry.S2ValidationQueries.S2LegacyValidQuery;
import com.google.errorprone.annotations.CanIgnoreReturnValue;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Serializable;
import java.util.AbstractList;
import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.concurrent.atomic.AtomicInteger;
import jsinterop.annotations.JsIgnore;
import jsinterop.annotations.JsMethod;
import jsinterop.annotations.JsType;
import org.jspecify.annotations.Nullable;

/**
 * An S2Polygon is an S2Region object that represents a polygon. A polygon is defined by zero or
 * more loops; recall that the interior of a loop is defined to be its left-hand side (see {@link
 * S2Loop}.) Unlike S2Point and S2Polyline, S2Polygons are mutable.
 *
 * <p>There are two different conventions for creating an S2Polygon:
 *
 * <p>First, {@link #initNested(List)} expects the input loops to be nested hierarchically. The
 * polygon interior then consists of the set of points contained by an odd number of loops. So for
 * example, a circular region with a hole in it would be defined as two CCW loops, with one loop
 * containing the other. The loops can be provided in any order.
 *
 * <p>When the orientation of the input loops is unknown, the nesting requirement is typically met
 * by calling {@link S2Loop#normalize()} on each loop (which inverts the loop if necessary so that
 * it encloses at most half the sphere). But in fact any set of loops can be used as long as (1)
 * there is no pair of loops that cross, and (2) there is no pair of loops whose union is the entire
 * sphere.
 *
 * <p>Second, {@link #initOriented(List)} expects the input loops to be oriented such that the
 * polygon interior is on the left-hand side of every loop. So for example, a circular region with a
 * hole in it would be defined using a CCW outer loop and a CW inner loop. The loop orientations
 * must all be consistent; for example, it is not valid to have one CCW loop nested inside another
 * CCW loop, because the region between the two loops is on the left-hand side of one loop and the
 * right-hand side of the other.
 *
 * <p>These requirements are not enforced in production code for performance reasons, but operations
 * on invalid polygons may return unexpected results. Clients are recommended to check any polygons
 * that might be invalid with {@link #isValid()} or {@link #findValidationError(S2Error)} before
 * operating on them.
 *
 * <p>When Java assertions are enabled, initNested(), initOriented(), and other public constructors
 * and factory methods will assert that the resulting polygon is valid. To deliberately construct
 * invalid polygons for testing, use {@link GeometryTestCase#uncheckedCreate(Supplier)} or {@link
 * GeometryTestCase#uncheckedInitialize(Runnable)}.
 *
 * <p>Most clients will not call these init methods directly; instead they should use {@link
 * S2Builder}, which has better support for dealing with imperfect data.
 *
 * <p>When the polygon is initialized, the given loops are automatically converted into a canonical
 * form consisting of "shells" and "holes". Shells and holes are both oriented CCW, and are nested
 * hierarchically. The loops are reordered to correspond to a preorder traversal of the nesting
 * hierarchy; initOriented may also invert some loops. The set of input S2Loop objects is always
 * preserved; the caller can use this to determine how the loops were reordered, if desired.
 *
 * <p>Polygons may represent any region of the sphere with a polygonal boundary, including the
 * entire sphere (known as the "full" polygon). The full polygon consists of a single canonical full
 * loop, as returned by {@link S2Loop#full()}, whereas the empty polygon has no loops at all.
 *
 * <p>Valid polygons have the following restrictions:
 *
 * <ul>
 *   <li>Vertices must be unit length.
 *   <li>Loops may not cross, i.e. the boundary of a loop may not intersect both the interior and
 *       exterior of any other loop.
 *   <li>Loops may not share edges, i.e. if a loop contains an edge AB, then no other loop may
 *       contain AB or BA.
 *   <li>Loops may share vertices, however no vertex may appear twice in a single loop (see {@link
 *       S2Loop}).
 *   <li>No loop may be empty. The full loop may appear only in the full polygon. Otherwise, loops
 *       must have at least three vertices.
 * </ul>
 *
 * <p>Note that S2Polygon has stricter validity rules for full and empty loops than {@link S2Loop}:
 * S2Polygon only accepts the *canonical* full and empty loops as valid, as returned by {@link
 * S2Loop#full()} and {@link S2Loop#empty()}, and considers any other single-vertex loop invalid. In
 * contrast, S2Loop considers every single-vertex loop to be valid, and either empty or full. The
 * S2Polygon convention is safer, as "accidental" single-vertex loops may arise from snapping or
 * other operations.
 *
 * <p>These restrictions are enforced by some constructors and methods when assertions are enabled.
 * Clients are recommended to check any polygons that might be invalid with {@link #isValid()}
 * before operating on them, otherwise unexpected results may be returned.
 *
 * @author shakusa@google.com (Steven Hakusa) ported from util/geometry
 * @author ericv@google.com (Eric Veach) original author
 */
@JsType
@SuppressWarnings({"Assertion", "IdentifierName"})
public final class S2Polygon implements S2Region, Comparable<S2Polygon>, Serializable {

  /** Version number of the lossless encoding format for S2Polygon. */
  private static final byte LOSSLESS_ENCODING_VERSION = 1;

  /** Version number of the compressed encoding format for S2Polygon. */
  private static final byte COMPRESSED_ENCODING_VERSION = 4;

  /**
   * Special S2Loop instance used as the "parent" of the polygon's outermost shell loops while
   * constructing the loop nesting hierarchy. The content is unimportant: comparison is by identity.
   */
  public static final S2Loop ROOT = new S2Loop(new S2Point[] {});

  /** Returns false for all shapes. */
  private static boolean reverseNone(S2Shape input) {
    return false;
  }

  /** Returns true for S2Loops for which {@link S2Loop#isHole()} is true. */
  private static boolean reverseHoles(S2Shape input) {
    if (input instanceof S2Loop) {
      S2Loop loop = (S2Loop) input;
      return loop.isHole();
    } else {
      return false;
    }
  }

  /**
   * The loops of this polygon. There is no total ordering of the loops, but a nested loop always
   * follows its containing loop, and all loops between parent and child are nested somewhere under
   * the parent.
   */
  private final List<S2Loop> loops = Lists.newArrayList();

  /**
   * {@code bound} is a conservative bound on all points contained by this polygon: If
   * A.contains(P), then A.bound.contains(new S2LatLng(P)).
   */
  private S2LatLngRect bound;

  /**
   * Since "bound" is not exact, it is possible that a polygon A contains another polygon B whose
   * bounds are slightly larger. "subregionBound" has been expanded sufficiently to account for this
   * error, i.e. if A.Contains(B), then A.subregionBound.contains(B.bound).
   */
  private S2LatLngRect subregionBound;

  /**
   * The spatial index for this S2Polygon. Contains either no shapes, or a single S2Polygon.Shape.
   */
  @VisibleForTesting transient S2ShapeIndex index;

  /**
   * In general we build the index the first time it is needed, but we make an exception for
   * contains(S2Point) because this method has a simple brute force implementation that is
   * relatively cheap. For this one method we keep track of the number of calls made and only build
   * the index once enough calls have been made that we think an index would be worthwhile.
   *
   * <p>NOTE: This counts downward: it is initialized to the maximum allowed based on the number of
   * vertices, decremented for each call, and an index is built when it reaches zero.
   */
  private final AtomicInteger unindexedContainsCalls = new AtomicInteger();

  /**
   * True if initOriented() was called and the given loops had inconsistent orientations (i.e., it
   * is not possible to construct a polygon such that the interior is on the left-hand side of all
   * loops). We need to remember this error so that it can be returned later by
   * findValidationError(), since it is not possible to detect this error once the polygon has been
   * initialized. This field is not preserved by encode/decode.
   */
  private boolean errorInconsistentLoopOrientations = false;

  /** Total number of vertices in all loops. */
  private int numVertices = 0;

  /**
   * Creates an empty polygon and then calls {@link #initNested(List)} with the given loops. The
   * resulting polygon is asserted to be valid, if assertions are enabled. Takes ownership of the
   * given loops, which must not be subsequently modified by the caller.
   *
   * <p>Among other requirements, valid loops must have at least three vertices, or be a *canonical*
   * full or empty loop. This is stricter than the S2Loop conventions. See the class Javadoc for
   * details.
   */
  public static S2Polygon fromLoops(List<S2Loop> loops) {
    S2Polygon polygon = new S2Polygon();
    polygon.uncheckedInitNested(loops);
    polygon.assertValid();
    return polygon;
  }

  /**
   * Creates an empty polygon and then calls {@link #initNested(List)} with the given loops. The
   * resulting polygon is asserted to be valid, if assertions are enabled. Takes ownership of the
   * given loops, which must not be subsequently modified by the caller.
   *
   * <p>Among other requirements, valid loops must have at least three vertices, or be a *canonical*
   * full or empty loop. This is stricter than the S2Loop conventions. See the class Javadoc for
   * details.
   */
  @JsMethod(name = "fromLoopArray")
  public static S2Polygon fromLoops(S2Loop... loops) {
    return fromLoops(Arrays.asList(loops));
  }

  /**
   * Creates an empty polygon. It can be made non-empty by calling {@link #init(List)} or one of the
   * other initialization methods.
   */
  @JsIgnore
  public S2Polygon() {
    bound = S2LatLngRect.empty();
    subregionBound = S2LatLngRect.empty();
    initIndex(); // No need to check validity on an empty polygon.
  }

  /**
   * Convenience constructor that creates an S2Polygon with a single loop corresponding to the given
   * cell.
   */
  @JsIgnore
  public S2Polygon(S2Cell cell) {
    loops.add(new S2Loop(cell));
    initOneLoop(); // No need to check validity for a loop created from a cell.
  }

  /**
   * Creates an empty polygon and then calls {@link #initNested(List)} with the given loops. Clears
   * the given list. The resulting polygon is asserted to be valid, if assertions are enabled. See
   * also {@link #fromLoops(List)} which is similar, but does not clear the given list.
   *
   * <p>Among other requirements, valid loops must have at least three vertices, or be a *canonical*
   * full or empty loop. This is stricter than the S2Loop conventions. See the class Javadoc for
   * details.
   */
  @JsIgnore
  public S2Polygon(List<S2Loop> loops) {
    uncheckedInitNested(loops);
    loops.clear();

    // TODO(user): Enable assertions that polygons are valid here.
    // assertValid();
  }

  /** Convenience constructor that calls {@link #initOneLoop(S2Loop)}. */
  @JsIgnore
  public S2Polygon(S2Loop loop) {
    initOneLoop(loop);
  }

  /**
   * Initializes this polygon from a single loop. Converts a canonical empty loop into a valid empty
   * polygon, which has no loops. A canonical full loop results in a valid full polygon. Otherwise,
   * the loop must have at least three vertices to be considered valid. This is stricter than the
   * S2Loop conventions. See the class Javadoc for details.
   *
   * <p>The resulting polygon is asserted to be valid, if assertions are enabled.
   */
  public void initOneLoop(S2Loop loop) {
    clearLoops();
    if (isCanonicalEmptyLoop(loop)) {
      // A canonical empty loop produces a valid empty polygon, with no loops.
      initLoopProperties();
      return;
    }

    loops.add(loop);
    initOneLoop();
    // TODO(user): Enable assertions that polygons are valid here.
    // assertValid();
  }

  /** Given that this.loops contains a single non-empty loop, initializes all other fields. */
  private void initOneLoop() {
    assert loops.size() == 1;
    S2Loop loop = loops.get(0);
    loop.setDepth(0);
    errorInconsistentLoopOrientations = false;
    numVertices = loop.numVertices();
    bound = loop.getRectBound();
    subregionBound = loop.getSubregionBound();
    initIndex();
  }

  /** Copy constructor. */
  @JsIgnore
  public S2Polygon(S2Polygon src) {
    copy(src);
  }

  /**
   * Returns true if the loop is a canonical empty loop, with a single vertex S2Loop.EMPTY_VERTEX,
   * as opposed to a single vertex loop that (accidentally) happens to be {@link S2Loop#isEmpty()}.
   * See also {@link #isCanonicalFullLoop(S2Loop)}.
   */
  @VisibleForTesting
  static boolean isCanonicalEmptyLoop(S2Loop loop) {
    return loop.numVertices() == 1 && loop.vertex(0).equalsPoint(S2Loop.EMPTY_VERTEX);
  }

  /**
   * Returns true if the loop is a canonical full loop, with a single vertex S2Loop.FULL_VERTEX, as
   * opposed to an single vertex loop that (accidentally) happens to be {@link S2Loop#isFull()}.
   *
   * <p>Other than canonical full and empty loops, loops must have at least three vertices to be
   * valid, according to S2Polygon conventions. S2Loop has looser conventions, and considers *any*
   * single-vertex loop to be valid, and either empty or full.
   */
  @VisibleForTesting
  static boolean isCanonicalFullLoop(S2Loop loop) {
    return loop.numVertices() == 1 && loop.vertex(0).equalsPoint(S2Loop.FULL_VERTEX);
  }

  /**
   * Returns true if the loop has fewer than three vertices, except if it is a canonical empty or
   * full loop.
   */
  static boolean isInvalidLoopTooFewVertices(S2Loop loop) {
    int numVertices = loop.numVertices();
    return numVertices == 0
        || (numVertices == 1 && !isCanonicalEmptyLoop(loop) && !isCanonicalFullLoop(loop))
        || numVertices == 2;
  }

  /**
   * Initializes this polygon to a copy of the given polygon, including making copies of the
   * (mutable) loops.
   */
  void copy(S2Polygon src) {
    clearLoops();
    for (int i = 0; i < src.numLoops(); ++i) {
      loops.add(new S2Loop(src.loop(i)));
    }

    // Don't copy errorInconsistentLoopOrientations, since this is not a property of the polygon,
    // but only of the way the polygon was constructed.
    this.numVertices = src.numVertices;
    this.unindexedContainsCalls.set(src.unindexedContainsCalls.get());
    this.bound = src.bound;
    this.subregionBound = src.subregionBound;
    initIndex();
  }

  /** Initializes the polygon index, which should contain no shapes. */
  private void initIndex() {
    // The index should be empty: either null, freshly created, or cleared by the caller.
    if (index != null) {
      assert index.getShapes().isEmpty();
    } else {
      index = new S2ShapeIndex();
    }
    index.add(shape());

    // See S2Loop for the details behind 'maxUnindexedContainsCalls'.
    int maxUnindexedContainsCalls;
    if (numVertices <= 8) {
      maxUnindexedContainsCalls = 10;
    } else if (numVertices <= 8192) {
      maxUnindexedContainsCalls = 50;
    } else if (numVertices <= 50000) {
      maxUnindexedContainsCalls = 10;
    } else {
      maxUnindexedContainsCalls = 2;
    }
    this.unindexedContainsCalls.set(maxUnindexedContainsCalls);
  }

  /**
   * Throws an AssertionError with a descriptive message, if this polygon is invalid, and assertions
   * are enabled. Only actually checks validity if assertions are enabled, so the cost of
   * findValidationError() is not incurred in production code.
   */
  void assertValid() {
    if (skipAssertions) {
      return;
    }

    assert assertValidHelper();
  }

  private boolean assertValidHelper() {
    S2Error error = new S2Error();
    if (findValidationError(error)) {
      throw new AssertionError(error.toString());
    }
    return true;
  }

  /** Returns the same instance after initializing transient fields. */
  @CanIgnoreReturnValue
  private Object readResolve() {
    initIndex();
    // We assume the serialized polygon was valid.
    return this;
  }

  /**
   * Two polygons are equal if and only if they have exactly the same loops. The loops must appear
   * in the same order, and corresponding loops must have the same linear vertex ordering (i.e.,
   * cyclic rotations are not allowed).
   */
  @Override
  public boolean equals(Object o) {
    if (o instanceof S2Polygon) {
      return equalsPolygon((S2Polygon) o);
    }
    return false;
  }

  /**
   * Return true if two polygons have exactly the same loops. The loops must appear in the same
   * order, and corresponding loops must have the same linear vertex ordering (i.e., cyclic
   * rotations are not allowed).
   */
  public boolean equalsPolygon(S2Polygon b) {
    if (numLoops() != b.numLoops()) {
      return false;
    }
    if (numVertices != b.numVertices) {
      return false;
    }
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop aLoop = loop(i);
      S2Loop bLoop = b.loop(i);
      if ((bLoop.depth() != aLoop.depth()) || !bLoop.equals(aLoop)) {
        return false;
      }
    }
    return true;
  }

  @Override
  public int hashCode() {
    return bound.hashCode();
  }

  /**
   * Comparator (needed by Comparable interface). For two polygons to be compared as equal:
   *
   * <ul>
   *   <li>They must have the same number of loops
   *   <li>The loops must be ordered in the same way (this is guaranteed by the total ordering
   *       imposed by {@link #sortValueLoops})
   *   <li>Loops must be logically equivalent (even if ordered with a different starting point, e.g.
   *       ABCD and BCDA).
   * </ul>
   */
  @Override
  public int compareTo(S2Polygon other) {
    // If number of loops differ, use that.
    if (this.numLoops() != other.numLoops()) {
      return this.numLoops() - other.numLoops();
    }
    for (int i = 0; i < this.numLoops(); ++i) {
      int compare = this.loops.get(i).compareTo(other.loops.get(i));
      if (compare != 0) {
        return compare;
      }
    }
    return 0;
  }

  /** An alias for {@link #initNested(List)}. */
  public void init(List<S2Loop> loops) {
    initNested(loops);
  }

  /**
   * Initializes this polygon from a set of hierarchically nested loops. The polygon interior
   * consists of the points contained by an odd number of loops. (Recall that a loop contains the
   * set of points on its left-hand side.) Asserts that the resulting polygon is valid, if
   * assertions are enabled.
   *
   * <p>This method takes ownership of the given loops and clears the given list. (See also {@link
   * #fromLoops(List)}, which does not clear the given list). It then figures out the loop nesting
   * hierarchy and assigns every loop a depth. Shells have even depths, and holes have odd depths.
   * Note that the loops are reordered so the hierarchy can be traversed more easily (see {@link
   * #getParent(int)}, {@link #getLastDescendant(int)}, and {@link S2Loop#depth()}). Loops for which
   * {@link #isCanonicalEmptyLoop(S2Loop)} are ignored, but otherwise invalid loops will result in
   * an invalid S2Polygon, and cause assertion failures, if assertions are enabled.
   *
   * <p>This method may be called more than once, in which case any existing loops are deleted
   * before being replaced by the input loops.
   *
   * <p>Among other requirements, valid loops must have at least three vertices, or be a *canonical*
   * full or empty loop. This is stricter than the S2Loop conventions. See the class Javadoc for
   * details.
   */
  public void initNested(List<S2Loop> nestedLoops) {
    uncheckedInitNested(nestedLoops);
    nestedLoops.clear();

    // TODO(user): Enable assertions that polygons are valid here.
    // assertValid();
  }

  /**
   * As above, but package private, does not assert that the polygon is valid, and does not clear or
   * modify the given list of loops. However, loops will be modified: they may be inverted and their
   * depth will be set.
   */
  void uncheckedInitNested(List<S2Loop> nestedLoops) {
    if (nestedLoops.size() == 1) {
      initOneLoop(nestedLoops.get(0));
      return;
    }

    clearLoops();

    // There are multiple loops, so nesting must be computed.
    IdentityHashMap<S2Loop, List<S2Loop>> loopMap = Maps.newIdentityHashMap();
    // Add an entry for the root of the map.
    loopMap.put(ROOT, new ArrayList<>());
    // Insert each loop into the loop map, reorganizing as needed at each step so that every loop is
    // mapped to its immediate nested children, and the root key has all the outermost
    // shells as children.
    for (S2Loop loop : nestedLoops) {
      if (!isCanonicalEmptyLoop(loop)) {
        insertLoop(loop, loopMap);
      }
    }

    // Sort all of the lists of loops; in this way we guarantee a total ordering on loops in the
    // polygon. Loops will be sorted by their natural ordering, while also preserving the
    // requirement that each loop is immediately followed by its descendants in the nesting
    // hierarchy.
    //
    // TODO(user): as per CL 18750833 code review comments: This should work for now, but I
    // think it's possible to guarantee the correct order inside insertLoop by searching for the
    // correct position in the children list before inserting.
    sortValueLoops(loopMap);

    // Build "this.loops" in depth-first traversal order from the loop map.
    initLoops(loopMap);

    // Computes numVertices, bound, subregionBound, and initializes the index.
    initLoopProperties();
  }

  /** Builds this.loops in depth-first traversal order from the loop map. */
  @SuppressWarnings("ReferenceEquality") // IdentityHashMap, ROOT is special.
  private void initLoops(IdentityHashMap<S2Loop, List<S2Loop>> loopMap) {
    ArrayDeque<S2Loop> loopStack = new ArrayDeque<>();
    loopStack.addFirst(ROOT);
    int depth = -1;
    while (!loopStack.isEmpty()) {
      S2Loop loop = loopStack.removeFirst();

      if (loop != ROOT) {
        depth = loop.depth();
        loops.add(loop);
      }
      List<S2Loop> children = loopMap.get(loop);
      for (int i = children.size() - 1; i >= 0; --i) {
        S2Loop child = children.get(i);
        assert child != null;
        child.setDepth(depth + 1);
        loopStack.addFirst(child);
      }
    }
  }

  /**
   * Like {@link #initNested(List)}, but expects loops to be oriented such that the polygon interior
   * is on the left-hand side of all loops. This implies that shells and holes should have opposite
   * orientations in the input to this method. (During initialization, loops representing holes will
   * automatically be inverted.) Asserts that the resulting polygon is valid, if assertions are
   * enabled.
   */
  public void initOriented(List<S2Loop> orientedLoops) {
    // Here is the algorithm:
    //
    // 1. Remember which of the given loops contain S2.origin().
    //
    // 2. Invert loops as necessary to ensure that they are nestable (i.e., no loop contains the
    //    complement of any other loop). This may result in a set of loops corresponding to the
    //    complement of the given polygon, but we will fix that problem later.
    //
    //    We make the loops nestable by first normalizing all the loops (i.e., inverting any loops
    //    whose turning angle is negative). This handles all loops except those whose turning angle
    //    is very close to zero (within the maximum error tolerance). Any such loops are inverted
    //    if and only if they contain S2.origin(). (In theory this step is only necessary if there
    //    are at least two such loops.)  The resulting set of loops is guaranteed to be nestable.
    //
    // 3. Build the polygon. This yields either the desired polygon or its complement.
    //
    // 4. If there is at least one loop, we find a loop L that is adjacent to S2.origin() (where
    //    "adjacent" means that there exists a path connecting S2.origin() to some vertex of L such
    //    that the path does not cross any loop). There may be a single such adjacent loop, or
    //    there may be several (in which case they should all have the same contains_origin()
    //    value). We choose L to be the loop containing the origin whose depth is greatest, or
    //    loop(0) (a top-level shell) if no such loop exists.
    //
    // 5. If (L originally contained origin) != (polygon contains origin), we invert the polygon.
    //    This is done by inverting a top-level shell whose turning angle is minimal and then fixing
    //    the nesting hierarchy. Note that because we normalized all the loops initially, this step
    //    is only necessary if the polygon requires at least one non-normalized loop to represent
    //    it.
    clearLoops();
    Set<S2Loop> containedOrigin = Sets.newIdentityHashSet();
    for (S2Loop loop : orientedLoops) {
      if (isCanonicalEmptyLoop(loop)) {
        continue;
      }
      if (loop.containsOrigin()) {
        containedOrigin.add(loop);
      }
      double angle = loop.getTurningAngle();
      if (abs(angle) > S2.getTurningAngleMaxError(loop.numVertices())) {
        // Normalize the loop.
        if (angle < 0) {
          loop.invert();
        }
      } else {
        // Ensure that the loop does not contain the origin.
        if (loop.containsOrigin()) {
          loop.invert();
        }
      }
    }

    initNested(orientedLoops); // initNested asserts validity, if assertions are enabled.

    if (numLoops() > 0) {
      S2Loop originLoop = loop(0);
      boolean polygonContainsOrigin = false;
      for (int i = 0; i < numLoops(); ++i) {
        if (loop(i).containsOrigin()) {
          polygonContainsOrigin ^= true;
          originLoop = loop(i);
        }
      }
      if (containedOrigin.contains(originLoop) != polygonContainsOrigin) {
        invert();
      }
    }

    // Verify that the original loops had consistent shell/hole orientations. Each original loop L
    // should have been inverted if and only if it now represents a hole.
    for (S2Loop loop : this.loops) {
      if ((containedOrigin.contains(loop) != loop.containsOrigin()) != loop.isHole()) {
        // There is no point in saving the loop index, because the error is a property of the entire
        // set of loops. In general there is no way to determine which ones are incorrect.
        errorInconsistentLoopOrientations = true;
        // Assertion of validity normally happens in initNested(), called above, but this error is
        // detected too late for that.
        // TODO(user): Enable assertions that polygons are valid here.
        // assert S2.skipAssertions || !errorInconsistentLoopOrientations;
      }
    }
  }

  /**
   * Computes numVertices, bound, subregionBound, and the index. Note this is not called by
   * decodeUncompressed.
   */
  private void initLoopProperties() {
    numVertices = 0;
    S2LatLngRect.Builder builder = S2LatLngRect.Builder.empty();
    for (S2Loop loop : loops) {
      if (loop.depth() == 0) {
        builder.union(loop.getRectBound());
      }
      numVertices += loop.numVertices();
    }
    bound = builder.build();
    subregionBound = S2EdgeUtil.RectBounder.expandForSubregions(bound);
    initIndex();
  }

  /**
   * Initializes a polygon from a set of {@link S2Loop}s. The loops must form a valid polygon, which
   * is checked if assertions are enabled. Among other requirements, valid loops must have at least
   * three vertices, or be a *canonical* full or empty loop. This is stricter than the S2Loop
   * conventions. See the class Javadoc for details.
   *
   * <p>Unlike {@link #init} this method assumes the caller already knows the nesting of loops
   * within other loops. The passed-in map maps from parents to their immediate child loops, with
   * {@code SD2Polygon.ROOT} mapping to the list of top-most shell loops. Immediate child loops must
   * be completely spatially contained within their parent loop, but not contained in any other
   * loop, except for ancestors of the parent. This method avoids the cost of determining nesting
   * internally, but if the passed in nesting is wrong, future operations on the S2Polygon may be
   * arbitrarily incorrect.
   *
   * <p>The passed-in loops become owned by the S2Polygon and should not be modified by the caller
   * after calling this method.
   *
   * @param nestedLoops loops with nesting.
   */
  public void initWithNestedLoops(IdentityHashMap<S2Loop, List<S2Loop>> nestedLoops) {
    // Now that we use ROOT instead of NULL to identify the root of the loop hierarchy, we have to
    // catch callers using the old convention.
    if (nestedLoops.containsKey(null)) {
      throw new IllegalArgumentException(
          "Use S2Polygon.ROOT, not null, to identify the root of the nesting hierarchy.");
    }

    clearLoops();
    initLoop(ROOT, -1, nestedLoops);

    // Empty the map as an indication we have taken ownership of the loops.
    nestedLoops.clear();

    initLoopProperties();
    // TODO(user): Enable assertions that polygons are valid here.
    // assertValid();
  }

  /** Appends the loops of this polygon to the given list and resets this polygon to be empty. */
  public void release(List<S2Loop> loops) {
    loops.addAll(this.loops);
    clearLoops();
  }

  /** Clears the loops list and other state, and resets the index. */
  private void clearLoops() {
    loops.clear();
    if (index != null) {
      index.reset();
    }
    numVertices = 0;
    bound = S2LatLngRect.empty();
    subregionBound = S2LatLngRect.empty();
    errorInconsistentLoopOrientations = false;
  }

  /**
   * Returns a new polygon containing explicitly the given loops, for testing. Allows construction
   * of invalid polygons, but if that is all you need, consider {@code
   * GeometryTestCase#unsafeCreate(() -> new S2Polygon(loops));} instead.
   *
   * <p>This factory method differs from the normal constructor in that it does not assert polygon
   * validity, does not skip empty loops, does not set loop depths, and does not reorder the given
   * loops according to nesting. It does set the polygon numVertices, bound, subregionBound, and
   * index as usual.
   */
  static S2Polygon fromExplicitLoops(List<S2Loop> loops) {
    S2Polygon polygon = new S2Polygon();
    polygon.clearLoops();
    for (S2Loop loop : loops) {
      polygon.loops.add(loop);
    }
    polygon.initLoopProperties();
    return polygon;
  }

  /**
   * Returns true if the given loops form a valid polygon, including checking whether the loops
   * themselves are valid.
   */
  @JsIgnore
  public static boolean isValid(List<S2Loop> loops) {
    S2Polygon maybeValid = new S2Polygon();
    maybeValid.uncheckedInitNested(Lists.newArrayList(loops));
    return maybeValid.isValid();
  }

  /**
   * Returns true if each loop on this polygon is valid, and if the relationships between all loops
   * are valid. Note that validity is checked automatically during initialization when Java
   * assertions are enabled, as they normally are for unit tests.
   *
   * <p>Specifically, this verifies that {@link S2Loop#isValid} is true for each {@link S2Loop}, and
   * additionally that loops do not cross each other, no loop is empty, the full loop only appears
   * in the full polygon, and that the loop nesting hierarchy is valid.
   *
   * <p>Note that an empty polygon is valid.
   */
  public boolean isValid() {
    S2Error error = new S2Error();
    return !findValidationError(error);
  }

  /**
   * Returns true if this is *not* a valid polygon and sets {@code error} appropriately. Otherwise,
   * returns false and leaves {@code error} unchanged.
   */
  public boolean findValidationError(S2Error error) {
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop loop = loop(i);
      if (isInvalidLoopTooFewVertices(loop)) {
        error.init(
            S2Error.Code.LOOP_NOT_ENOUGH_VERTICES,
            "Loop %d: has only %d vertices, at least three are required.",
            i,
            loop.numVertices());
        return true;
      }
    }

    // S2LegacyValidQuery doesn't have access to the loop depth information via the S2Shape API,
    // so we'll validate it manually. We just have to check that the depth values are non-negative,
    // and that we don't skip depths.
    int lastDepth = -1;
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop loop = loop(i);
      int depth = loop.depth();
      if (depth < 0 || depth > lastDepth + 1) {
        error.init(
            S2Error.Code.POLYGON_INVALID_LOOP_DEPTH,
            Platform.formatString("Loop %d: invalid loop depth (%d)", i, depth));
        return true;
      }
      lastDepth = depth;

      // S2LegacyValidQuery will go into an infinite loop if there are NaN vertices. Check they
      // are unit-length like S2Loop.findValidationErrorNoIndex does for the no-validation-query
      // case.
      for (int j = 0; j < loop.numVertices(); ++j) {
        if (!S2.isUnitLength(loop.vertex(j))) {
          error.init(
              S2Error.Code.NOT_UNIT_LENGTH,
              Platform.formatString("Loop %d: Vertex %d is not unit length", i, j));
          return true;
        }
      }
    }

    // Check whether initOriented detected inconsistent loop orientations.
    if (errorInconsistentLoopOrientations) {
      error.init(
          S2Error.Code.POLYGON_INCONSISTENT_LOOP_ORIENTATIONS,
          "Inconsistent loop orientations detected");
      return true;
    }

    // Finally, check that the geometry is topologically valid.
    S2LegacyValidQuery query = new S2LegacyValidQuery();
    return !query.validate(index, error);
  }

  /** Returns true if this polygon contains no loops, and thus no area and no edges. */
  public boolean isEmpty() {
    return loops.isEmpty();
  }

  /**
   * Returns true if this is the full polygon, which has a single canonical full loop. That is, a
   * loop with a single vertex that is equal to the single vertex in {@link S2Loop#full()}. Note
   * that this is stricter than {@link S2Loop#isFull()}, which returns true for any single-vertex
   * loop where the vertex has a negative z-value.
   */
  public boolean isFull() {
    return loops.size() == 1 && isCanonicalFullLoop(loops.get(0));
  }

  /**
   * Returns the number of loops in this polygon. A full polygon has one loop. An empty polygon has
   * no loops.
   */
  public int numLoops() {
    return loops.size();
  }

  /**
   * Returns the total number of vertices in all loops. A full polygon has one vertex, {@code
   * S2Loop.FULL_VERTEX}. An empty polygon has no vertices.
   */
  public int numVertices() {
    return numVertices;
  }

  /**
   * Returns the loop at the given index. Note that during initialization, the given loops are
   * reordered according to a preorder traversal of the loop nesting hierarchy. This implies that
   * every loop is immediately followed by its descendants. This hierarchy can be traversed using
   * the methods {@link #getParent(int)}, {@link #getLastDescendant(int)}, and {@link
   * S2Loop#depth()}.
   */
  public S2Loop loop(int k) {
    return loops.get(k);
  }

  /** Returns a view of the list of {@link S2Loop}s that make up this S2Polygon. */
  public List<S2Loop> getLoops() {
    return new AbstractList<S2Loop>() {
      @Override
      public int size() {
        return loops.size();
      }

      @Override
      public S2Loop get(int index) {
        return loop(index);
      }
    };
  }

  /** Returns the index of this polygon. */
  public S2ShapeIndex index() {
    return index;
  }

  /** Returns the index of the parent of loop {@code k}, or -1 if it has no parent. */
  public int getParent(int k) {
    int depth = loop(k).depth();
    if (depth == 0) {
      return -1; // Optimization.
    }
    while (--k >= 0 && loop(k).depth() >= depth) {
      // spin
    }
    return k;
  }

  /**
   * Returns the index of the last loop that is contained within loop {@code k}. Returns {@code
   * numLoops() - 1} if {@code k < 0}. Note that loops are indexed according to a preorder traversal
   * of the nesting hierarchy, so the immediate children of loop {@code k} can be found by iterating
   * over loops {@code (k+1)..getLastDescendant(k)} and selecting those whose depth is equal to
   * {@code (loop(k).depth() + 1)}.
   */
  public int getLastDescendant(int k) {
    if (k < 0) {
      return numLoops() - 1;
    }
    int depth = loop(k).depth();
    while (++k < numLoops() && loop(k).depth() > depth) {
      // spin
    }
    return k - 1;
  }

  /**
   * Note that the S2Point centroid, if computed, is not normalized. It is simply the sum of all the
   * loop centroids, which are scaled by loop area.
   */
  private S2AreaCentroid getAreaCentroid(boolean doCentroid) {
    double areaSum = 0;
    S2Point centroidSum = S2Point.ZERO;
    for (int i = 0; i < numLoops(); ++i) {
      S2AreaCentroid areaCentroid = doCentroid ? loop(i).getAreaAndCentroid() : null;
      double loopArea = doCentroid ? areaCentroid.getArea() : loop(i).getArea();

      int loopSign = loop(i).sign();
      areaSum += loopSign * loopArea;
      if (doCentroid && !loop(i).isEmptyOrFull()) {
        S2Point currentCentroid = areaCentroid.getCentroid();
        if (currentCentroid != null) {
          centroidSum =
              new S2Point(
                  centroidSum.x + loopSign * currentCentroid.x,
                  centroidSum.y + loopSign * currentCentroid.y,
                  centroidSum.z + loopSign * currentCentroid.z);
        }
      }
    }

    return new S2AreaCentroid(areaSum, doCentroid ? centroidSum : null);
  }

  /**
   * Returns the area of the polygon interior, i.e. the region on the left side of an odd number of
   * loops (the area is between 0 and 4*Pi) and the true centroid of the polygon, weighted by the
   * area of the polygon (see s2.h for details on centroids). Note that the centroid might not be
   * contained by the polygon, and is not normalized.
   */
  public S2AreaCentroid getAreaAndCentroid() {
    return getAreaCentroid(true);
  }

  /**
   * Returns the area of the polygon interior, i.e. the region on the left side of an odd number of
   * loops. The return value is between 0 and 4*Pi.
   */
  public double getArea() {
    return getAreaCentroid(false).getArea();
  }

  /**
   * Returns the true centroid of the polygon, multiplied by the area of the polygon (see S2.java
   * for details on centroids). Note that the centroid might not be contained by the polygon, and is
   * not normalized.
   */
  public S2Point getCentroid() {
    return getAreaCentroid(true).getCentroid();
  }

  /**
   * If all of the polygon's vertices happen to be the centers of S2Cells at some level, then
   * returns that level, otherwise returns -1. See also {@link #initToSnapped(S2Polygon, int)} and
   * {@link S2BuilderSnapFunctions.S2CellIdSnapFunction}. Returns -1 if the polygon has no vertices.
   */
  public int getSnapLevel() {
    int snapLevel = -1;
    for (S2Loop loop : loops) {
      for (int j = 0; j < loop.numVertices(); j++) {
        S2Point p = loop.vertex(j);
        FaceSiTi faceSiTi = S2Projections.xyzToFaceSiTi(p);
        int level = S2Projections.levelIfCenter(faceSiTi, p);
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
    }
    return snapLevel;
  }

  /**
   * Computes the level at which most of the vertices are snapped. If multiple levels have the same
   * maximum number of vertices snapped to it, the first one (lowest level number / largest area /
   * smallest encoding length) will be chosen, so this is desired. Returns -1 for unsnapped
   * polygons.
   *
   * <p>See also {@link #initToSnapped(S2Polygon, int)} and {@link
   * S2BuilderSnapFunctions.S2CellIdSnapFunction}.
   */
  int getBestSnapLevel() {
    int[] histogram = new int[S2CellId.MAX_LEVEL + 1];
    for (S2Loop loop : loops) {
      for (S2Point p : loop.vertices()) {
        FaceSiTi faceSiTi = S2Projections.xyzToFaceSiTi(p);
        int level = S2Projections.levelIfCenter(faceSiTi, p);
        // Level is -1 for unsnapped points.
        if (level >= 0) {
          histogram[level]++;
        }
      }
    }
    int snapLevel = 0;
    for (int i = 1; i < histogram.length; i++) {
      if (histogram[i] > histogram[snapLevel]) {
        snapLevel = i;
      }
    }
    if (histogram[snapLevel] == 0 && !isEmpty()) {
      // This is an unsnapped polygon.
      return -1;
    }
    return snapLevel;
  }

  /**
   * Returns the shortest distance from a point P to this polygon, given as the angle formed between
   * P, the origin, and the nearest point on the polygon to P. This angle in radians is equivalent
   * to the arclength along the unit sphere. The given point must be normalized (unit length).
   *
   * <p>If the point is contained by the polygon, the distance returned is 0. If the polygon is
   * empty, returns {@link S1Angle.INFINITY}.
   */
  public S1Angle getDistance(S2Point p) {
    if (contains(p)) {
      return S1Angle.radians(0);
    }
    return getDistanceToBoundary(p);
  }

  /**
   * Returns the distance from the given point to the polygon boundary. If the polygon is empty or
   * full, returns {@link S1Angle.INFINITY}, since the polygon has no boundary. The given point must
   * be unit length.
   */
  public S1Angle getDistanceToBoundary(S2Point p) {
    S2ClosestEdgeQuery.Builder builder = S2ClosestEdgeQuery.builder().setIncludeInteriors(false);
    S2ClosestEdgeQuery.Query query = builder.build(index);
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target = new S2ClosestEdgeQuery.PointTarget<>(p);
    return query.getDistance(target).toAngle();
  }

  /**
   * Returns the overlap fraction of polygon b on polygon a, i.e. the ratio of area of intersection
   * to the area of polygon a. As a special case, if polygon 'a' is empty, the result is 1.
   */
  public static double getOverlapFraction(S2Polygon a, S2Polygon b) {
    S2Polygon intersection = new S2Polygon();
    intersection.initToIntersection(a, b);
    double intersectionArea = intersection.getArea();
    double aArea = a.getArea();
    if (aArea > 0) {
      return intersectionArea >= aArea ? 1 : intersectionArea / aArea;
    } else {
      // Special case: consider 0/0 to be 1.
      return (intersectionArea == 0) ? 1 : 0;
    }
  }

  /**
   * Old implementation of getOverlapFraction, which will be removed when all callers have migrated.
   *
   * @deprecated Use {@link #getOverlapFraction(S2Polygon, S2Polygon)} instead.
   */
  @Deprecated
  public static double getOverlapFractionOld(S2Polygon a, S2Polygon b) {
    S2Polygon intersection = new S2Polygon();
    intersection.initToIntersectionOld(a, b);
    double intersectionArea = intersection.getArea();
    double aArea = a.getArea();
    if (aArea > 0) {
      return intersectionArea >= aArea ? 1 : intersectionArea / aArea;
    } else {
      // Special case: consider 0/0 to be 1.
      return (intersectionArea == 0) ? 1 : 0;
    }
  }

  /**
   * Returns the IoU of polygon b and polygon a, i.e. the ratio of area of intersection to the area
   * of the union of polygon a and polygon b.
   */
  public static double getIntersectionOverUnion(S2Polygon a, S2Polygon b) {
    S2Polygon intersection = new S2Polygon();
    intersection.initToIntersection(a, b);
    double intersectionArea = intersection.getArea();
    double unionArea = union(ImmutableList.of(a, b)).getArea();
    if (unionArea > 0) {
      return intersectionArea / unionArea;
    } else {
      return 0;
    }
  }

  /**
   * If the given point is contained by the polygon, it is returned. Otherwise, returns the closest
   * point on the polygon boundary to the given point. If the polygon is empty, returns the given
   * point.
   *
   * <p>Note that the result may or may not be contained by the polygon. The given point must be
   * normalized (unit length).
   */
  public S2Point project(S2Point p) {
    if (contains(p)) {
      return p;
    }
    return projectToBoundary(p);
  }

  /**
   * Return the closest point on the polygon boundary to the given point. If the polygon is empty or
   * full, return the input argument since the polygon has no boundary. The given point "x" should
   * be unit length.
   */
  public S2Point projectToBoundary(S2Point x) {
    if (isEmpty() || isFull()) {
      return x;
    }
    S2ClosestEdgeQuery.Builder builder = S2ClosestEdgeQuery.builder().setIncludeInteriors(false);
    S2ClosestEdgeQuery.Query query = builder.build(index);
    S2ClosestEdgeQuery.PointTarget<S1ChordAngle> target = new S2ClosestEdgeQuery.PointTarget<>(x);
    Optional<S2BestEdgesQueryBase.Result<S1ChordAngle>> edge = query.findClosestEdge(target);
    // A closest edge should always be found, since this polygon is not empty or full, and therefore
    // has a boundary.
    Preconditions.checkState(edge.isPresent());
    return query.project(x, edge.get());
  }

  /**
   * Returns true if this polygon contains the given other polygon, i.e., if polygon A contains all
   * points contained by polygon B. This is the old implementation, which will be removed as soon as
   * all callers have migrated.
   *
   * @deprecated Use {@link #contains(S2Polygon)} instead.
   */
  @Deprecated
  public boolean containsOld(S2Polygon b) {
    // It's worth checking bounding rectangles, since they are precomputed. Note that the first
    // bound has been expanded to account for possible numerical errors in the second bound.

    if (!subregionBound.contains(b.getRectBound())) {
      // It is possible that A contains B even though bound(A) does not contain bound(B). This can
      // only happen when polygon B has at least two outer shells and the union of the two bounds
      // spans all longitudes. For example, suppose that B consists of two shells with a longitude
      // gap between them, while A consists of one shell that surrounds both shells of B but goes
      // the other way around the sphere (so that it does not intersect the longitude gap).
      //
      // For convenience we just check whether B has at least two loops rather than two outer
      // shells.
      if (b.numLoops() == 1 || !bound.lng().union(b.bound.lng()).isFull()) {
        return false;
      }
    }

    // The following case is not handled by S2BooleanOperation because it only determines whether
    // the boundary of the result is empty (which does not distinguish between the full and empty
    // polygons).
    if (isEmpty() && b.isFull()) {
      return false;
    }

    // Polygon A contains B iff B does not intersect the complement of A. From the intersection
    // algorithm below, this means that the complement of A must exclude the entire boundary of B,
    // and B must exclude all shell boundaries of the complement of A. (It can be shown that B must
    // then exclude the entire boundary of the complement of A.)  The first call below returns false
    // if the boundaries cross, therefore the second call does not need to check for any crossing
    // edges (which makes it cheaper).
    return containsBoundary(b) && b.excludesNonCrossingComplementShells(this);
  }

  /**
   * Returns true if this polygon contains the given other polygon, i.e., if polygon A contains all
   * points contained by polygon B.
   */
  public boolean contains(S2Polygon b) {
    // It's worth checking bounding rectangles, since they are precomputed. Note that the first
    // bound has been expanded to account for possible numerical errors in the second bound.
    if (!subregionBound.contains(b.getRectBound())) {
      // It is possible that A contains B even though bound(A) does not contain bound(B). This can
      // only happen when polygon B has at least two outer shells and the union of the two bounds
      // spans all longitudes. For example, suppose that B consists of two shells with a longitude
      // gap between them, while A consists of one shell that surrounds both shells of B but goes
      // the other way around the sphere (so that it does not intersect the longitude gap).
      //
      // For convenience we just check whether B has at least two loops rather than two outer
      // shells.
      if (b.numLoops() == 1 || !bound.lng().union(b.bound.lng()).isFull()) {
        return false;
      }
    }

    // The following case is not handled by S2BooleanOperation because it only determines whether
    // the boundary of the result is empty (which does not distinguish between the full and empty
    // polygons).
    if (isEmpty() && b.isFull()) {
      return false;
    }

    return S2BooleanOperation.contains(index, b.index);
  }

  /**
   * Returns true if this polygon (A) approximately contains the given other polygon (B). This is
   * true if it is possible to move the vertices of B no further than "vertexMergeRadius" such that
   * A contains the modified B.
   *
   * <p>For example, the empty polygon will contain any polygon whose maximum width is no more than
   * vertexMergeRadius.
   */
  public boolean approxContains(S2Polygon b, S1Angle vertexMergeRadius) {
    S2Polygon difference = new S2Polygon();
    difference.initToDifference(b, this, new IdentitySnapFunction(vertexMergeRadius));
    return difference.numLoops() == 0;
  }

  /**
   * Returns true if this polygon (A) approximately contains the given other polygon (B). This is
   * true if it is possible to move the vertices of B no further than "vertexMergeRadius" such that
   * A contains the modified B.
   *
   * <p>For example, the empty polygon will contain any polygon whose maximum width is no more than
   * vertexMergeRadius.
   *
   * @deprecated Use {@link #approxContains(S2Polygon, S1Angle)} instead.
   */
  @Deprecated
  public boolean approxContainsOld(S2Polygon b, S1Angle vertexMergeRadius) {
    S2Polygon difference = new S2Polygon();
    difference.initToDifferenceSloppy(b, this, vertexMergeRadius);
    return difference.numLoops() == 0;
  }

  /**
   * Returns true if this polygon intersects the given other polygon, i.e., if there is a point that
   * is contained by both polygons. Note that polygons merely "touching" at a common point or shared
   * (but reversed) edge are not considered intersecting. This is the old implementation, which will
   * be removed as soon as all callers have migrated.
   *
   * @deprecated Use {@link #intersects(S2Polygon)} instead.
   */
  @Deprecated
  public boolean intersectsOld(S2Polygon b) {
    // If both polygons have one loop, use the more efficient S2Loop method.
    // Note that S2Loop.intersects does its own bounding rectangle check.
    if (numLoops() == 1 && b.numLoops() == 1) {
      return loop(0).intersects(b.loop(0));
    }

    // Otherwise if neither polygon has holes, we can still use the more efficient S2Loop.intersects
    // method. The polygons intersect if and only if some pair of loop regions intersect.
    if (!bound.intersects(b.getRectBound())) {
      return false;
    }
    if (!hasHoles() && !b.hasHoles()) {
      for (S2Loop loop : b.loops) {
        if (anyLoopIntersects(loop)) {
          return true;
        }
      }
      return false;
    }

    // Polygon A is disjoint from B if A excludes the entire boundary of B and B excludes all shell
    // boundaries of A. (It can be shown that B must then exclude the entire boundary of A.)  The
    // first call below returns false if the boundaries cross, therefore the second call does not
    // need to check for crossing edges.
    return !excludesBoundary(b) || !b.excludesNonCrossingShells(this);
  }

  /**
   * Returns true if this polygon intersects the given other polygon, i.e., if there is a point that
   * is contained by both polygons. Note that polygons merely "touching" at a common point or shared
   * (but reversed) edge are not considered intersecting.
   */
  public boolean intersects(S2Polygon b) {
    // It's worth checking bounding rectangles, since they are precomputed.
    if (!this.bound.intersects(b.bound)) {
      return false;
    }

    // The following case is not handled by S2BooleanOperation because it only determines whether
    // the boundary of the result is empty (which does not distinguish between the full and empty
    // polygons).
    if (isFull() && b.isFull()) {
      return true;
    }
    return S2BooleanOperation.intersects(index, b.index);
  }

  /**
   * Clips the boundary of A to the interior of B, and adds the resulting edges to {@code builder}.
   * Shells are directed CCW and holes are directed clockwise. If {@code reverseA} is true, these
   * directions are reversed in polygon A. If {@code invertB} is true, the boundary of A is clipped
   * to the exterior rather than the interior of B. If {@code addSharedEdges} is true, then the
   * output will include any edges that are shared between A and B (both edges must be in the same
   * direction after any edge reversals are taken into account).
   *
   * @deprecated This method is only used by other deprecated methods.
   */
  @Deprecated
  private static void clipBoundary(
      final S2Polygon a,
      boolean reverseA,
      final S2Polygon b,
      boolean invertB,
      boolean addSharedEdges,
      S2PolygonBuilder builder) {
    EdgeClipper clipper = new EdgeClipper(b.index, addSharedEdges, S2Polygon::reverseHoles);
    List<ParametrizedS2Point> intersections = Lists.newArrayList();
    for (S2Loop aLoop : a.loops) {
      int n = aLoop.numVertices();
      int dir = (aLoop.isHole() ^ reverseA) ? -1 : 1;
      boolean inside = b.contains(aLoop.vertex(0)) ^ invertB;
      for (int j = (dir > 0) ? 0 : n; n > 0; --n, j += dir) {
        S2Point a0 = aLoop.vertex(j);
        S2Point a1 = aLoop.vertex(j + dir);
        clipper.clipEdge(a0, a1, intersections);

        if (inside) {
          intersections.add(new ParametrizedS2Point(0, a0));
        }
        inside = ((intersections.size() & 0x1) == 0x1);
        assert (b.contains(a1) ^ invertB == inside);
        if (inside) {
          intersections.add(new ParametrizedS2Point(1, a1));
        }

        Collections.sort(intersections);
        for (int k = 0; k < intersections.size(); k += 2) {
          S2Point x = intersections.get(k).getPoint();
          S2Point y = intersections.get(k + 1).getPoint();
          if (x.equalsPoint(y)) {
            continue;
          }
          builder.addEdge(x, y);
        }
        intersections.clear();
      }
    }
  }

  /** Returns the total number of vertices in all loops. */
  public int getNumVertices() {
    return this.numVertices;
  }

  /** Initializes this polygon to the complement of the given polygon. */
  public void initToComplement(S2Polygon a) {
    clearLoops();
    copy(a);
    invert();
  }

  /**
   * Returns a new polygon that describes the outline of the given cell union. In principle this
   * polygon should exactly contain the cell union and this polygon's inverse should not intersect
   * the cell union, but rounding issues may cause this not to be the case.
   */
  public static S2Polygon fromCellUnionBorder(S2CellUnion cells) {
    // We use S2Builder to compute the union. Due to rounding errors, we can't compute an exact
    // union - when a small cell is adjacent to a larger cell, the shared edges can fail to line up
    // exactly. Two cell edges cannot come closer then MIN_WIDTH, so if we have S2Builder snap edges
    // within half that distance, then we should always merge shared edges without merging different
    // edges.
    S1Angle snapRadius = S1Angle.radians(0.5 * MIN_WIDTH.getValue(S2CellId.MAX_LEVEL));

    S2Builder builder = new S2Builder.Builder(new IdentitySnapFunction(snapRadius)).build();
    S2PolygonLayer layer = new S2PolygonLayer();
    builder.startLayer(layer);
    // If there are no loops, either the cell union is empty, or it consists of all six faces. In
    // latter case, the result should be the full polygon rather than the empty one.
    builder.addIsFullPolygonPredicate(g -> !cells.isEmpty());

    for (S2CellId id : cells) {
      builder.addLoop(new S2Loop(new S2Cell(id)));
    }

    S2Error error = new S2Error();
    boolean unused = builder.build(error);
    assert error.ok() : error.text();
    return layer.getPolygon();
  }

  /**
   * Initializes the polygon from input polygon "a" using the given S2Builder. If the result has an
   * empty boundary (no loops), also decides whether the result should be the full polygon rather
   * than the empty one based on the area of the input polygon.
   */
  private void initFromBuilder(final S2Polygon a, S2Builder builder) {
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(this);
    builder.startLayer(layer);
    // If there are no loops, the result should be the full polygon rather than the empty one if the
    // input area and bounds are more than a hemisphere.
    builder.addIsFullPolygonPredicate(g -> a.bound.area() > 2 * PI && a.getArea() > 2 * PI);
    builder.addPolygon(a);
    S2Error error = new S2Error();
    boolean unused = builder.build(error);
    // TODO(torrey): Perhaps this should throw S2Exception, rather than only asserting? C++ only
    // does LOG(DFATAL), FWIW.
    assert error.ok() : error.text();
  }

  /**
   * Snaps the vertices of the given polygon using the given SnapFunction. For example, {@code
   * IntLatLngSnapFunction(6)} snaps to E6 coordinates. This can change the polygon topology
   * (merging loops, for example), but the resulting polygon is guaranteed to be valid, and no
   * vertex will move by more than snapFunction.snapRadius(). See {@link S2Builder} for other
   * guarantees (e.g., minimum edge-vertex separation).
   *
   * <p>Note that this method is a thin wrapper over S2Builder, so if you are starting with data
   * that is not in S2Polygon format (e.g., integer E7 coordinates) then it is faster to just use
   * S2Builder directly.
   */
  @JsMethod(name = "initToSnapFunction")
  public void initToSnapped(final S2Polygon input, SnapFunction snapFunction) {
    S2Builder builder = new S2Builder.Builder(snapFunction).build();
    initFromBuilder(input, builder);
  }

  /**
   * Convenience function that snaps the vertices to S2CellId centers at the given level. If
   * snapping to S2 leaves, at MAX_LEVEL (30), this will simplify the polygon with a tolerance of
   * {@code S2Projections.maxDiag.getValue(S2CellId.MAX_LEVEL)}, or approximately 0.13 microdegrees,
   * or 1.5cm on the surface of the Earth. Polygons can be efficiently encoded after they have been
   * snapped.
   */
  public void initToSnapped(final S2Polygon polygon, int snapLevel) {
    initToSnapped(polygon, new S2CellIdSnapFunction(snapLevel));
  }

  /** Inverts this polygon (replacing it by its complement.) */
  private void invert() {
    // Inverting any one loop will invert the polygon. The best loop to invert is the one whose
    // area is largest, since this yields the smallest area after inversion. The loop with the
    // largest area is always at depth 0. The descendants of this loop all have their depth reduced
    // by 1, while the former siblings of this loop all have their depth increased by 1.

    // The empty and full polygons are handled specially.
    if (isEmpty()) {
      loops.add(S2Loop.full());
    } else if (isFull()) {
      clearLoops();
    } else {
      // Find the loop whose area is largest (i.e., whose turning angle is smallest), minimizing
      // calls to getTurningAngle(). In particular, for polygons with a single shell at level 0
      // there is not need to call getTurningAngle() at all. (This method is relatively expensive.)
      int best = -1;
      double bestAngle = 0;
      for (int i = 1; i < numLoops(); ++i) {
        S2Loop loop = loop(i);
        if (loop.depth() == 0) {
          // We defer computing the turning angle of loop 0 until we discover that the polygon has
          // another top-level shell.
          if (best == -1) {
            best = 0;
            bestAngle = loop(best).getTurningAngle();
          }
          double angle = loop.getTurningAngle();
          if (angle < bestAngle) {
            best = i;
            bestAngle = angle;
          }
        }
      }

      if (best < 0) {
        best = 0;
      }

      // Build the new loops list, starting with the inverted loop.
      loop(best).invert();
      List<S2Loop> newLoops = Lists.newArrayListWithCapacity(numLoops());
      newLoops.add(loop(best));

      // Add the former siblings of this loop as descendants.
      int lastBest = getLastDescendant(best);
      for (int i = 0; i < numLoops(); ++i) {
        if (i < best || i > lastBest) {
          S2Loop loop = loop(i);
          loop.setDepth(loop.depth() + 1);
          newLoops.add(loop);
        }
      }

      // Add the former children of this loop as siblings.
      for (int i = 0; i < numLoops(); ++i) {
        if (i > best && i <= lastBest) {
          S2Loop loop = loop(i);
          loop.setDepth(loop.depth() - 1);
          newLoops.add(loop);
        }
      }
      Preconditions.checkState(loops.size() == newLoops.size());
      clearLoops();
      loops.addAll(newLoops);
    }

    index.reset();
    initLoopProperties();
    // Asserting validity of the resulting polygon is not necessary, assuming the input was valid.
  }

  /**
   * Initializes this polygon to the result of the given boolean operation on the two provided
   * polygons, returning true on success, or returning false and setting 'error' on failure.
   *
   * <p>The polygons must be valid. This is checked, if assertions are enabled.
   */
  private boolean initToOperation(
      S2BooleanOperation.OpType opType,
      SnapFunction snapFunction,
      final S2Polygon a,
      final S2Polygon b,
      S2Error error) {
    // S2BooleanOperation on invalid inputs is likely to fail.
    a.assertValid();
    b.assertValid();

    S2BooleanOperation.Builder builder = new S2BooleanOperation.Builder(snapFunction);
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(this);
    S2BooleanOperation op = builder.build(opType, layer);
    return op.build(a.index(), b.index(), error);
  }

  /**
   * Initializes this polygon to the result of the given boolean operation on the two provided
   * polygons. Returns true on success, false on failure. The polygons must be valid.
   *
   * @throws AssertionError if input polygons are invalid or the operation fails, and Java
   *     assertions are enabled.
   */
  private boolean initToOperation(
      S2BooleanOperation.OpType opType,
      SnapFunction snapFunction,
      final S2Polygon a,
      final S2Polygon b) {
    S2Error error = new S2Error();
    boolean success = initToOperation(opType, snapFunction, a, b, error);
    assert success : "initToOperation failed: " + error;
    return success;
  }

  /**
   * Initializes this polygon to the intersection of the given two polygons. Uses a default snap
   * function, which is the IdentitySnapFunction with a snap radius of {@link
   * S2EdgeUtil#INTERSECTION_MERGE_RADIUS} (equal to about 1.8e-15 radians or 11 nanometers on the
   * Earth's surface). This means that vertices may be positioned arbitrarily, but vertices that are
   * extremely close together can be merged together. The reason for a non-zero default snap radius
   * is that it helps to eliminate narrow cracks and slivers when T-vertices are present. For
   * example, adjacent S2Cells at different levels do not share exactly the same boundary, so there
   * can be a narrow crack between them. If a polygon is intersected with those cells and the pieces
   * are unioned together, the result would have a narrow crack unless the snap radius is set to a
   * non-zero value.
   *
   * <p>Note that if you want to encode the vertices in a lower-precision representation (such as
   * S2CellIds or E7), it is much better to use a suitable SnapFunction rather than rounding the
   * vertices yourself, because this will create self-intersections unless you ensure that the
   * vertices and edges are sufficiently well-separated first. In particular you need to use a snap
   * function whose {@link SnapFunction#minEdgeVertexSeparation()} is at least twice the maximum
   * distance that a vertex can move when rounded.
   *
   * <p>Returns true on success, false otherwise. Should not fail if the polygons are valid.
   *
   * @throws AssertionError if the operation fails and Java assertions are enabled.
   */
  @CanIgnoreReturnValue
  public boolean initToIntersection(final S2Polygon a, final S2Polygon b) {
    return initToIntersection(a, b, new IdentitySnapFunction(S2EdgeUtil.INTERSECTION_MERGE_RADIUS));
  }

  /**
   * Initializes this polygon to the intersection of the given two polygons, using the given
   * "snapFunction". Returns true if successful, false otherwise. Should not fail if the polygons
   * are valid.
   *
   * @throws AssertionError if the operation fails and Java assertions are enabled.
   */
  @JsMethod(name = "initToSnappedIntersection")
  @CanIgnoreReturnValue
  public boolean initToIntersection(
      final S2Polygon a, final S2Polygon b, SnapFunction snapFunction) {
    if (!a.bound.intersects(b.bound)) {
      initNested(new ArrayList<>());
      return true;
    }
    return initToOperation(S2BooleanOperation.OpType.INTERSECTION, snapFunction, a, b);
  }

  /**
   * Initializes this polygon to the intersection of the given two polygons, using the given
   * "snapFunction".
   *
   * <p>The snapFunction allows you to specify a minimum spacing between output vertices, and/or
   * that the vertices should be snapped to a discrete set of points (e.g. S2CellId centers or E7
   * lat/lng coordinates). Any snap function can be used, including the IdentitySnapFunction with a
   * snapRadius of zero (which preserves the input vertices exactly).
   *
   * <p>The boundary of the output polygon before snapping is guaranteed to be accurate to within
   * {@link S2EdgeUtil#INTERSECTION_ERROR} of the exact result. Snapping can move the boundary by an
   * additional distance that depends on the snap function. Finally, any degenerate portions of the
   * output polygon are automatically removed (i.e., regions that do not contain any points) since
   * S2Polygon does not allow such regions.
   *
   * <p>See {@link S2Builder} and {@link S2SnapFunctions} for more details on snap functions. For
   * example, you can snap to E7 coordinates by setting "snapFunction" to {@code new
   * S2BuilderSnapFunctions.IntLatLngSnapFunction(7)}.
   *
   * <p>Returns true on success, or returns false and sets "error" appropriately otherwise. However,
   * this function should never return an error provided that both input polygons are valid (i.e.,
   * {@link #isValid()} returns true).
   */
  @JsMethod(name = "initToSnappedIntersectionOrError")
  @CanIgnoreReturnValue
  public boolean initToIntersection(
      final S2Polygon a, final S2Polygon b, SnapFunction snapFunction, S2Error error) {
    if (!a.bound.intersects(b.bound)) {
      initNested(new ArrayList<>());
      return true;
    }
    return initToOperation(S2BooleanOperation.OpType.INTERSECTION, snapFunction, a, b, error);
  }

  /**
   * Initializes this polygon to the intersection of the given two polygons. Uses a default merge
   * radius, which is just large enough to compensate for errors that occur when computing
   * intersection points between edges ({@link S2EdgeUtil#DEFAULT_INTERSECTION_TOLERANCE}).
   *
   * <p>If you are going to convert the resulting polygon to a lower-precision format, it is
   * necessary to increase the merge radius. See {@link #initToIntersectionSloppy(S2Polygon,
   * S2Polygon, S1Angle)} below.
   *
   * @deprecated Use {@link #initToIntersection(S2Polygon, S2Polygon) or one of the variants above.
   */
  @Deprecated
  public void initToIntersectionOld(final S2Polygon a, final S2Polygon b) {
    initToIntersectionSloppy(a, b, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  /**
   * Initializes this polygon to the intersection of the given two polygons. The {@code
   * vertexMergeRadius} determines how close two vertices must be to be merged together and how
   * close a vertex must be to an edge in order to be spliced into it (see {@link S2PolygonBuilder}
   * for details).
   *
   * <p>If you are going to convert the resulting polygon to a lower-precision format, it is
   * necessary to set the merge radius appropriately in order to get a valid result after rounding
   * (i.e., no duplicate vertices, etc.) For example, if you are going to convert them to {@code
   * geostore.PolygonProto} format, then {@code S1Angle.e7(1)} is a good value for {@code
   * vertexMergeRadius}.
   *
   * @deprecated Use {@link #initToIntersection(S2Polygon, S2Polygon, SnapFunction) or one of the
   * variants above. An IdentitySnapFunction(vertexMergeRadius) will produce similar results.
   */
  @Deprecated
  public void initToIntersectionSloppy(
      final S2Polygon a, final S2Polygon b, S1Angle vertexMergeRadius) {
    clearLoops();
    if (a.isEmpty() || b.isEmpty()) {
      // From the check above we are already empty.
    } else if (a.isFull()) {
      copy(b);
    } else if (b.isFull()) {
      copy(a);
    } else {
      if (!a.bound.intersects(b.bound)) {
        return;
      }

      // We want the boundary of A clipped to the interior of B, plus the boundary of B clipped to
      // the interior of A, plus one copy of any directed edges that are in both boundaries.

      S2PolygonBuilder.Options options =
          S2PolygonBuilder.Options.DIRECTED_XOR.toBuilder()
              .setMergeDistance(vertexMergeRadius)
              .build();
      S2PolygonBuilder builder = new S2PolygonBuilder(options);
      clipBoundary(a, false, b, false, true, builder);
      clipBoundary(b, false, a, false, false, builder);
      builder.assemblePolygon(this, null);

      // If the result had a non-empty boundary then we are done. Unfortunately, if the boundary is
      // empty then there are two possible results: the empty polygon or the full polygon. This
      // choice would be trivial to resolve except for the existence of "vertex_merge_radius" and
      // also numerical errors when computing edge intersection points. In particular:
      //
      //  - The intersection of two non-full polygons may be full. For example, one or both
      //    polygons may have tiny cracks that are eliminated due to vertex merging/edge splicing.
      //
      //  - The intersection of two polygons that both contain S2.origin() (or any other point) may
      //    be empty. For example, both polygons may have tiny shells that surround the common
      //    point but that are eliminated.
      //
      //  - Even before any vertex merging/edge splicing, the computed boundary edges are not useful
      //    in distinguishing almost-full polygons from almost-empty due to numerical errors in
      //    computing edge intersections. Such errors can reverse the orientation of narrow cracks
      //    or slivers.
      //
      // So instead we fall back to heuristics. Essentially we compute the minimum and maximum
      // intersection area based on the areas of the two input polygons. If only one of {0, 4*Pi}
      // is possible then we return that result. If neither is possible (before vertex merging,
      // etc) then we return the one that is closest to being possible. (It never true that
      // both are possible.)
      if (numLoops() == 0) {
        // We know that both polygons are non-empty due to the initial bounds check. By far the
        // most common case is that the intersection is empty, so we want to make that case fast.
        // The intersection area satisfies:
        //
        //   max(0, A + B - 4*Pi) <= Intersection(A, B) <= min(A, B)
        //
        // where A, B refer to a polygon and/or its area. Note that if either A or B is at most
        // 2*Pi, the result must be "empty". We can use the bounding rectangle areas as upper
        // bounds on the polygon areas.
        if (a.bound.area() <= 2 * PI || b.bound.area() <= 2 * PI) {
          return;
        }
        double aArea = a.getArea();
        double bArea = b.getArea();
        double minArea = max(0.0, aArea + bArea - 4 * PI);
        double maxArea = min(aArea, bArea);
        if (minArea > 4 * PI - maxArea) {
          invert();
        }
      }
    }
  }

  /**
   * Initializes this polygon to the union of the given two polygons. Uses a default snap function,
   * which is the IdentitySnapFunction with a snap radius of {@link
   * S2EdgeUtil#INTERSECTION_MERGE_RADIUS} (equal to about 1.8e-15 radians or 11 nanometers on the
   * Earth's surface). This means that vertices may be positioned arbitrarily, but vertices that are
   * extremely close together can be merged together. The reason for a non-zero default snap radius
   * is that it helps to eliminate narrow cracks and slivers when T-vertices are present. For
   * example, adjacent S2Cells at different levels do not share exactly the same boundary, so there
   * can be a narrow crack between them. If a polygon is intersected with those cells and the pieces
   * are unioned together, the result would have a narrow crack unless the snap radius is set to a
   * non-zero value.
   *
   * <p>Note that if you want to encode the vertices in a lower-precision representation (such as
   * S2CellIds or E7), it is much better to use a suitable SnapFunction rather than rounding the
   * vertices yourself, because this will create self-intersections unless you ensure that the
   * vertices and edges are sufficiently well-separated first. In particular you need to use a snap
   * function whose {@link SnapFunction#minEdgeVertexSeparation()} is at least twice the maximum
   * distance that a vertex can move when rounded.
   *
   * <p>Returns true on success, false otherwise. Should not fail if the polygons are valid.
   *
   * @throws AssertionError if the operation fails or an input polygon is invalid, and Java
   *     assertions are enabled.
   */
  @CanIgnoreReturnValue
  public boolean initToUnion(final S2Polygon a, final S2Polygon b) {
    return initToUnion(a, b, new IdentitySnapFunction(S2EdgeUtil.INTERSECTION_MERGE_RADIUS));
  }

  /**
   * Like {@link #initToUnion(S2Polygon, S2Polygon, SnapFunction, S2Error)} below, initializes this
   * polygon to the union of the given two polygons. Returns true on success, false otherwise.
   * Should not fail if the polygons are valid.
   *
   * @throws AssertionError if the operation fails or an input polygon is invalid, and Java
   *     assertions are enabled.
   */
  @JsMethod(name = "initToSnappedUnion")
  @CanIgnoreReturnValue
  public boolean initToUnion(final S2Polygon a, final S2Polygon b, SnapFunction snapFunction) {
    return initToOperation(S2BooleanOperation.OpType.UNION, snapFunction, a, b);
  }

  /**
   * Initializes this polygon to the union of the given two polygons, and snaps the result using the
   * given "snapFunction". The snapFunction allows you to specify a minimum spacing between output
   * vertices, and/or that the vertices should be snapped to a discrete set of points (e.g. S2CellId
   * centers or E7 lat/lng coordinates). Any snap function can be used, including the
   * IdentitySnapFunction with a snapRadius of zero (which preserves the input vertices exactly).
   *
   * <p>The boundary of the output polygon before snapping is guaranteed to be accurate to within
   * {@link S2EdgeUtil#INTERSECTION_ERROR} of the exact result. Snapping can move the boundary by an
   * additional distance that depends on the snap function. Finally, any degenerate portions of the
   * output polygon are automatically removed (i.e., regions that do not contain any points) since
   * S2Polygon does not allow such regions.
   *
   * <p>See {@link S2Builder} and {@link S2SnapFunctions} for more details on snap functions. For
   * example, you can snap to E7 coordinates by setting "snapFunction" to {@code new
   * S2BuilderSnapFunctions.IntLatLngSnapFunction(7)}.
   *
   * <p>Returns true on success, or returns false and sets "error" appropriately otherwise. However,
   * this function should never return an error provided that both input polygons are valid (i.e.,
   * {@link #isValid()} returns true). The provided polygons are checked for validity when
   * assertions are enabled.
   */
  @JsMethod(name = "initToSnappedUnionOrError")
  public boolean initToUnion(
      final S2Polygon a, final S2Polygon b, SnapFunction snapFunction, S2Error error) {
    return initToOperation(S2BooleanOperation.OpType.UNION, snapFunction, a, b, error);
  }

  /**
   * Initializes this polygon to the union of the given two polygons, and snaps the result so that
   * vertices within "vertexMergeRadius" of each other are merged together. The provided polygons
   * must be valid, which is is checked when assertions are enabled.
   *
   * @throws IllegalArgumentException if initToUnion fails, which should only happen if an input
   *     polygon is invalid.
   */
  @JsMethod(name = "initToUnionMergeRadius")
  public void initToUnion(final S2Polygon a, final S2Polygon b, S1Angle vertexMergeRadius) {
    S2Error error = new S2Error();
    boolean success = initToUnion(a, b, new IdentitySnapFunction(vertexMergeRadius), error);
    if (!success) {
      throw new IllegalArgumentException(error.text());
    }
  }

  /**
   * Initializes this polygon to the union of the given two polygons. The {@code vertexMergeRadius}
   * determines how close two vertices must be to be merged together and how close a vertex must be
   * to an edge in order to be spliced into it (see {@link S2PolygonBuilder} for details).
   *
   * <p>If you are going to convert the resulting polygon to a lower-precision format, it is
   * necessary to set the merge radius appropriately in order to get a valid result after rounding
   * (i.e., no duplicate vertices, etc.) For example, if you are going to convert them to {@code
   * geostore.PolygonProto} format, then {@code S1Angle.e7(1)} is a good value for {@code
   * vertexMergeRadius}.
   *
   * @deprecated Use {@link #initToUnion(S2Polygon, S2Polygon, SnapFunction)}. A SnapFunction of
   *     IdentitySnapFunction(vertexMergeRadius) should provide similar results.
   */
  @Deprecated
  public void initToUnionSloppy(final S2Polygon a, final S2Polygon b, S1Angle vertexMergeRadius) {
    clearLoops();
    if (a.isEmpty() || b.isFull()) {
      copy(b);
    } else if (b.isEmpty() || a.isFull()) {
      copy(a);
    } else {
      // We want the boundary of A clipped to the exterior of B, plus the boundary of B clipped to
      // the exterior of A, plus one copy of any directed edges that are in both boundaries.
      S2PolygonBuilder.Options options =
          S2PolygonBuilder.Options.DIRECTED_XOR.toBuilder()
              .setMergeDistance(vertexMergeRadius)
              .build();
      S2PolygonBuilder builder = new S2PolygonBuilder(options);
      clipBoundary(a, false, b, true, true, builder);
      clipBoundary(b, false, a, true, false, builder);
      builder.assemblePolygon(this, null);

      if (numLoops() == 0) {
        // See comments in InitToIntersectionSloppy(). In this case, the union area satisfies:
        //
        //   max(A, B) <= Union(A, B) <= min(4*Pi, A + B)
        //
        // The most common case is that neither input polygon is empty, but the union is empty due
        // to vertex merging/simplification.
        if (a.bound.area() + b.bound.area() <= 2 * PI) {
          return;
        }
        double aArea = a.getArea();
        double bArea = b.getArea();
        double minArea = max(aArea, bArea);
        double maxArea = min(4 * PI, aArea + bArea);
        if (minArea > 4 * PI - maxArea) {
          invert();
        }
      }
    }
  }

  /**
   * Initializes this polygon to the difference (A - B) of the given two polygons. Uses a default
   * snap function, which is the IdentitySnapFunction with a snap radius of {@link
   * S2EdgeUtil#INTERSECTION_MERGE_RADIUS} (equal to about 1.8e-15 radians or 11 nanometers on the
   * Earth's surface). This means that vertices may be positioned arbitrarily, but vertices that are
   * extremely close together can be merged together. The reason for a non-zero default snap radius
   * is that it helps to eliminate narrow cracks and slivers when T-vertices are present. For
   * example, adjacent S2Cells at different levels do not share exactly the same boundary, so there
   * can be a narrow crack between them. If a polygon is intersected with those cells and the pieces
   * are unioned together, the result would have a narrow crack unless the snap radius is set to a
   * non-zero value.
   *
   * <p>Note that if you want to encode the vertices in a lower-precision representation (such as
   * S2CellIds or E7), it is much better to use a suitable SnapFunction rather than rounding the
   * vertices yourself, because this will create self-intersections unless you ensure that the
   * vertices and edges are sufficiently well-separated first. In particular you need to use a snap
   * function whose {@link SnapFunction#minEdgeVertexSeparation()} is at least twice the maximum
   * distance that a vertex can move when rounded.
   *
   * <p>Returns true on success, false otherwise. Should not fail if the polygons are valid.
   *
   * @throws AssertionError if the operation fails and Java assertions are enabled.
   */
  @CanIgnoreReturnValue
  public boolean initToDifference(final S2Polygon a, final S2Polygon b) {
    return initToDifference(a, b, new IdentitySnapFunction(S2EdgeUtil.INTERSECTION_MERGE_RADIUS));
  }

  /**
   * Like {@link #initToDifference(S2Polygon, S2Polygon, SnapFunction)} below, initializes this
   * polygon to the difference (A - B) of the given two polygons, using the given "snapFunction".
   *
   * <p>Returns true on success, false otherwise. Should not fail if the polygons are valid.
   *
   * @throws AssertionError if the operation fails and Java assertions are enabled.
   */
  @JsMethod(name = "initToSnappedDifference")
  @CanIgnoreReturnValue
  public boolean initToDifference(final S2Polygon a, final S2Polygon b, SnapFunction snapFunction) {
    return initToOperation(S2BooleanOperation.OpType.DIFFERENCE, snapFunction, a, b);
  }

  /**
   * Initializes this polygon to the difference (A - B) of the given two polygons, using the given
   * "snapFunction". The snapFunction allows you to specify a minimum spacing between output
   * vertices, and/or that the vertices should be snapped to a discrete set of points (e.g. S2CellId
   * centers or E7 lat/lng coordinates). Any snap function can be used, including the
   * IdentitySnapFunction with a snapRadius of zero (which preserves the input vertices exactly).
   *
   * <p>The boundary of the output polygon before snapping is guaranteed to be accurate to within
   * {@link S2EdgeUtil#INTERSECTION_ERROR} of the exact result. Snapping can move the boundary by an
   * additional distance that depends on the snap function. Finally, any degenerate portions of the
   * output polygon are automatically removed (i.e., regions that do not contain any points) since
   * S2Polygon does not allow such regions.
   *
   * <p>See {@link S2Builder} and {@link S2SnapFunctions} for more details on snap functions. For
   * example, you can snap to E7 coordinates by setting "snapFunction" to {@code new
   * S2BuilderSnapFunctions.IntLatLngSnapFunction(7)}.
   *
   * <p>Returns true on success, or returns false and sets "error" appropriately otherwise. However,
   * note that this function should never return an error provided that both input polygons are
   * valid (i.e., {@link #isValid()} returns true).
   */
  @JsMethod(name = "initToSnappedDifferenceOrError")
  @CanIgnoreReturnValue
  public boolean initToDifference(
      final S2Polygon a, final S2Polygon b, SnapFunction snapFunction, S2Error error) {
    return initToOperation(S2BooleanOperation.OpType.DIFFERENCE, snapFunction, a, b, error);
  }

  /**
   * Initializes this polygon to the difference (A - B) of the given two polygons. The {@code
   * vertexMergeRadius} determines how close two vertices must be to be merged together and how
   * close a vertex must be to an edge in order to be spliced into it (see {@link S2PolygonBuilder}
   * for details).
   *
   * <p>If you are going to convert the resulting polygon to a lower-precision format, it is
   * necessary to set the merge radius appropriately in order to get a valid result after rounding
   * (i.e., no duplicate vertices, etc.) For example, if you are going to convert them to {@code
   * geostore.PolygonProto} format, then {@code S1Angle.e7(1)} is a good value for {@code
   * vertexMergeRadius}.
   *
   * @deprecated Use {@link #initToDifference(S2Polygon, S2Polygon, SnapFunction)}. A SnapFunction
   *     of IdentitySnapFunction(vertexMergeRadius) should provide similar results.
   */
  @Deprecated
  public void initToDifferenceSloppy(
      final S2Polygon a, final S2Polygon b, S1Angle vertexMergeRadius) {
    clearLoops();
    if (a.isEmpty() || b.isFull()) {
      // From the check above, this polygon is already empty.
    } else if (a.isFull()) {
      initToComplement(b);
    } else {
      // Note that we cannot short circuit the b.isEmpty() case because even something and nothing
      // might have no difference if the difference falls within the merge distance.

      // We want the boundary of A clipped to the exterior of B, plus the reversed boundary of B
      // clipped to the interior of A, plus one copy of any edge in A that is also a reverse edge in
      // B.
      S2PolygonBuilder.Options options =
          S2PolygonBuilder.Options.DIRECTED_XOR.toBuilder()
              .setMergeDistance(vertexMergeRadius)
              .build();
      S2PolygonBuilder builder = new S2PolygonBuilder(options);
      clipBoundary(a, false, b, true, true, builder);
      clipBoundary(b, true, a, false, false, builder);
      builder.assemblePolygon(this, null);

      if (numLoops() == 0) {
        // See comments in InitToIntersectionSloppy(). In this case, the difference area satisfies:
        //
        //   max(0, A - B) <= Difference(A, B) <= min(A, 4*Pi - B)
        //
        // By far the most common case is that result is empty.
        if (a.bound.area() <= 2 * PI || b.bound.area() >= 2 * PI) {
          return;
        }
        double aArea = a.getArea();
        double bArea = b.getArea();
        double minArea = max(0.0, aArea - bArea);
        double maxArea = min(aArea, 4 * PI - bArea);
        if (minArea > 4 * PI - maxArea) {
          invert();
        }
      }
    }
  }

  /**
   * Initializes this polygon to the symmetric difference of the given two polygons. Uses a default
   * snap function, which is the IdentitySnapFunction with a snap radius of {@link
   * S2EdgeUtil#INTERSECTION_MERGE_RADIUS} (equal to about 1.8e-15 radians or 11 nanometers on the
   * Earth's surface). This means that vertices may be positioned arbitrarily, but vertices that are
   * extremely close together can be merged together. The reason for a non-zero default snap radius
   * is that it helps to eliminate narrow cracks and slivers when T-vertices are present. For
   * example, adjacent S2Cells at different levels do not share exactly the same boundary, so there
   * can be a narrow crack between them. If a polygon is intersected with those cells and the pieces
   * are unioned together, the result would have a narrow crack unless the snap radius is set to a
   * non-zero value.
   *
   * <p>Note that if you want to encode the vertices in a lower-precision representation (such as
   * S2CellIds or E7), it is much better to use a suitable SnapFunction rather than rounding the
   * vertices yourself, because this will create self-intersections unless you ensure that the
   * vertices and edges are sufficiently well-separated first. In particular you need to use a snap
   * function whose {@link SnapFunction#minEdgeVertexSeparation()} is at least twice the maximum
   * distance that a vertex can move when rounded.
   *
   * <p>Returns true on success, false otherwise. Should not fail if the polygons are valid.
   *
   * @throws AssertionError if the operation fails and Java assertions are enabled.
   */
  @CanIgnoreReturnValue
  public boolean initToSymmetricDifference(final S2Polygon a, final S2Polygon b) {
    return initToSymmetricDifference(
        a, b, new IdentitySnapFunction(S2EdgeUtil.INTERSECTION_MERGE_RADIUS));
  }

  /**
   * As {@link #initToSymmetricDifference(S2Polygon, S2Polygon, SnapFunction)} below, initializes
   * this polygon to the symmetric difference of the given two polygons. Returns true on success,
   * false otherwise. Should not fail if the polygons are valid.
   *
   * @throws AssertionError if the operation fails and Java assertions are enabled.
   */
  @JsMethod(name = "initToSnappedSymmetricDifference")
  @CanIgnoreReturnValue
  public boolean initToSymmetricDifference(
      final S2Polygon a, final S2Polygon b, SnapFunction snapFunction) {
    return initToOperation(S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE, snapFunction, a, b);
  }

  /**
   * initializes this polygon to the symmetric difference of the given two polygons, using the given
   * "snapFunction". The snapFunction allows you to specify a minimum spacing between output
   * vertices, and/or that the vertices should be snapped to a discrete set of points (e.g. S2CellId
   * centers or E7 lat/lng coordinates). Any snap function can be used, including the
   * IdentitySnapFunction with a snapRadius of zero (which preserves the input vertices exactly).
   *
   * <p>The boundary of the output polygon before snapping is guaranteed to be accurate to within
   * {@link S2EdgeUtil#INTERSECTION_ERROR} of the exact result. Snapping can move the boundary by an
   * additional distance that depends on the snap function. Finally, any degenerate portions of the
   * output polygon are automatically removed (i.e., regions that do not contain any points) since
   * S2Polygon does not allow such regions.
   *
   * <p>See {@link S2Builder} and {@link S2SnapFunctions} for more details on snap functions. For
   * example, you can snap to E7 coordinates by setting "snapFunction" to {@code new
   * S2BuilderSnapFunctions.IntLatLngSnapFunction(7)}.
   *
   * <p>Returns true on success, or returns false and sets "error" appropriately otherwise. However,
   * note that this function should never return an error provided that both input polygons are
   * valid (i.e., {@link #isValid()} returns true).
   */
  @JsMethod(name = "initToSnappedSymmetricDifferenceOrError")
  @CanIgnoreReturnValue
  public boolean initToSymmetricDifference(
      final S2Polygon a, final S2Polygon b, SnapFunction snapFunction, S2Error error) {
    return initToOperation(
        S2BooleanOperation.OpType.SYMMETRIC_DIFFERENCE, snapFunction, a, b, error);
  }

  /**
   * Initializes this polygon to a polygon that contains fewer vertices and is within tolerance of
   * the polygon a, with some caveats.
   *
   * <p>If {@code snapToCellCenters} is true, the vertices of this polygon will be snapped to the
   * centers of cells at the smallest level that is guaranteed to result in a valid polygon given
   * the specified tolerance.
   *
   * <ul>
   *   <li>If there is a very small island in the original polygon, it may disappear completely.
   *       Thus some parts of the original polygon may not be close to the simplified polygon. Those
   *       parts are small, though, and arguably don't need to be kept.
   *   <li>However, if there are dense islands, they may all disappear, instead of replacing them by
   *       a big simplified island.
   *   <li>For small tolerances (compared to the polygon size), it may happen that the simplified
   *       polygon has more vertices than the original, if the first step of the simplification
   *       creates too many self-intersections. One can construct unrealistic cases where that
   *       happens to an extreme degree.
   * </ul>
   */
  public void initToSimplified(S2Polygon a, S1Angle tolerance, boolean snapToCellCenters) {
    SnapFunction snapFunction;
    if (snapToCellCenters) {
      // Get the cell level such that vertices will not move by more than 'tolerance'.
      int level = S2CellIdSnapFunction.levelForMaxSnapRadius(tolerance);
      snapFunction = new S2CellIdSnapFunction(level);
    } else {
      // Will snap to points whose distance from the input point is no greater than 'tolerance'.
      snapFunction = new IdentitySnapFunction(tolerance);
    }

    initToSimplified(a, snapFunction);
  }

  /**
   * Returns this polygon if {@link #getNumVertices()} is under 'maxVertices', or returns a copy of
   * this polygon that has been simplified to at most 'maxVertices' distinct vertices. Note that if
   * the given 'maxVertices' is less than 24, then 24 is used instead.
   *
   * <p>The result is simplified using {@link S2CellIdSnapFunction} at the highest level that snaps
   * vertices enough to get under the given limit. The result may be smaller than requested, and may
   * be quantized in appearance due to the rectilinear nature of the S2CellId coordinates at a given
   * level.
   */
  public S2Polygon simplify(int maxVertices) {
    maxVertices = max(24, maxVertices);
    if (getNumVertices() <= maxVertices) {
      return this;
    }
    Iterable<S2Point> vertices = concat(transform(getLoops(), S2Loop::vertices));
    for (int level = S2CellId.MAX_LEVEL; level >= 0; level--) {
      S2CellIdSnapFunction snapFunction = new S2CellIdSnapFunction(level);
      int n = 0;
      S2Point last = S2Point.ZERO;
      for (S2Point v : vertices) {
        v = snapFunction.snapPoint(v);
        if (!v.equalsPoint(last)) {
          last = v;
          n++;
          if (n > maxVertices) {
            break; // try the next level
          }
        }
      }
      if (n <= maxVertices) {
        S2Polygon simplified = new S2Polygon();
        simplified.initToSimplified(this, snapFunction);
        return simplified;
      }
    }
    throw new IllegalStateException("Unable to simplify");
  }

  /**
   * Snaps the input polygon according to the given "snapFunction" and reduces the number of
   * vertices if possible, while ensuring that no vertex moves further than
   * snapFunction.snapRadius().
   *
   * <p>A zero snap radius will leave the input geometry unmodified.
   *
   * <p>Simplification works by replacing nearly straight chains of short edges with longer edges,
   * in a way that preserves the topology of the input polygon up to the creation of degeneracies.
   * This means that loops or portions of loops may become degenerate, in which case they are
   * removed.
   *
   * <p>For example, if there is a very small island in the original polygon, it may disappear
   * completely. (Even if there are dense islands, they could all be removed rather than being
   * replaced by a larger simplified island if more area is covered by water than land.)
   *
   * <p>What's more, since we snap at the same time that we simplify, edges that come within the
   * snap radius of a vertex may have a vertex inserted resulting in a "pinch" that forces S2Builder
   * to produce multiple output loops.
   */
  @JsMethod(name = "initToSimplifiedSnapFunction")
  public void initToSimplified(S2Polygon a, SnapFunction snapFunction) {
    S2Builder builder = new S2Builder.Builder(snapFunction).setSimplifyEdgeChains(true).build();
    initFromBuilder(a, builder);
  }

  /**
   * Like {@link #initToSimplified}, except that any vertices or edges on the boundary of the given
   * S2Cell are preserved if possible. This method requires that the polygon has already been
   * clipped so that it does not extend outside the cell by more than "boundaryTolerance". In other
   * words, it operates on polygons that have already been intersected with a cell.
   *
   * <p>Typically this method is used in geometry-processing pipelines that intersect polygons with
   * a collection of S2Cells and then process those cells in parallel, where each cell generates
   * some geometry that needs to be simplified. In contrast, if you just need to simplify the
   * _input_ geometry then it is easier and faster to do the simplification before computing the
   * intersection with any S2Cells.
   *
   * <p>"boundaryTolerance" specifies how close a vertex must be to the cell boundary to be kept.
   * The default tolerance of 1e-15 is large enough to handle any reasonable way of interpolating
   * points along the cell boundary, such as S2.getIntersection(), S2.interpolate(), or direct (u,v)
   * interpolation using S2.faceUVtoXYZ(). However, if the vertices have been snapped to a
   * lower-precision representation (e.g., S2CellId centers or E7 coordinates) then you will need to
   * set this tolerance explicitly. For example, if the vertices were snapped to E7 coordinates then
   * "boundaryTolerance" should be set to {@code IntLatLngSnapFunction.minSnapRadiusForExponent(7);}
   *
   * <p>Degenerate portions of loops are always removed, so if a vertex on the cell boundary belongs
   * only to degenerate regions then it will not be kept. For example, if the input polygon is a
   * narrow strip of width less than "snapRadius" along one side of the cell, then the entire loop
   * may become degenerate and be removed.
   *
   * <p>REQUIRES: all vertices of "a" are within "boundaryTolerance" of "cell".
   */
  public void initToSimplifiedInCell(
      S2Polygon a, final S2Cell cell, S1Angle snapRadius, S1Angle boundaryTolerance) {
    // The polygon to be simplified consists of "boundary edges" that follow the cell boundary and
    // "interior edges" that do not. We want to simplify the interior edges while leaving the
    // boundary edges unchanged. It's not sufficient to call S2Builder.forceVertex() on all boundary
    // vertices. For example, suppose the polygon includes a triangle ABC where all three vertices
    // are on the cell boundary and B is a cell corner. Then if interior edge AC snaps to vertex B,
    // this loop would become degenerate and be removed. Similarly, we don't want boundary edges to
    // snap to interior vertices, since this also would cause portions of the polygon along the
    // boundary to be removed.
    //
    // Instead we use a two-pass algorithm. In the first pass, we simplify *only* the interior
    // edges, using forceVertex() to ensure that any edge endpoints on the cell boundary do not
    // move. In the second pass, we add the boundary edges (which are guaranteed to still form loops
    // with the interior edges) and build the output polygon.
    //
    // Note that in theory, simplifying the interior edges could create an intersection with one of
    // the boundary edges, since if two interior edges intersect very near the boundary then the
    // intersection point could be slightly outside the cell (by at most INTERSECTION_ERROR). This
    // is the *only* way that a self-intersection can be created, and it is expected to be extremely
    // rare. Nevertheless we use a small snap radius in the second pass in order to eliminate any
    // such self-intersections.
    //
    // We also want to preserve the cyclic vertex order of loops, so that the original polygon can
    // be reconstructed when no simplification is possible (i.e., idempotency). In order to do this,
    // we group the input edges into a sequence of polylines. Each polyline contains only one type
    // of edge (interior or boundary). We use S2Builder to simplify the interior polylines, while
    // the boundary polylines are passed through unchanged. Each interior polyline is in its own
    // S2Builder layer in order to keep the edges in sequence. This lets us ensure that in the
    // second pass, the edges are added in their original order so that S2PolygonLayer can
    // reconstruct the original loops.

    // We want an upper bound on how much "u" or "v" can change when a point on the boundary of the
    // S2Cell is moved away by up to "boundaryTolerance". Inverting this, instead we could compute a
    // lower bound on how far a point can move away from an S2Cell edge when "u" or "v" is changed
    // by a given amount. The latter quantity is simply (MIN_WIDTH.deriv() / 2) under the
    // S2_LINEAR_PROJECTION model, where we divide by 2 because we want the bound in terms of
    // {@code (u = 2 * s - 1)} rather than "s" itself. Consulting s2metrics.cc, this value is
    // sqrt(2/3)/2 = sqrt(1/6). Going back to the original problem, this gives:
    double boundaryToleranceUV = sqrt(6) * boundaryTolerance.radians();

    // The first pass yields a collection of simplified polylines that preserve the original cyclic
    // vertex order.
    List<S2Polyline> polylines = simplifyEdgesInCell(a, cell, boundaryToleranceUV, snapRadius);

    // The second pass eliminates any intersections between interior edges and boundary edges, and
    // then assembles the edges into a polygon.
    S2Builder.Builder options =
        new S2Builder.Builder(
            new IdentitySnapFunction(S1Angle.radians(S2EdgeUtil.INTERSECTION_ERROR)));
    options.setIdempotent(false); // Force snapping up to the given radius
    S2Builder builder = options.build();
    S2PolygonLayer layer = new S2PolygonLayer();
    layer.setPolygon(this);
    builder.startLayer(layer);
    // If there are no loops, the result should be the full polygon rather than the empty one if the
    // input area and bounds are more than a hemisphere.
    builder.addIsFullPolygonPredicate(g -> a.bound.area() > 2 * PI && a.getArea() > 2 * PI);

    for (S2Polyline polyline : polylines) {
      builder.addPolyline(polyline);
    }

    S2Error error = new S2Error();
    if (!builder.build(error)) {
      throw new IllegalArgumentException("Could not build polygon: " + error);
    }
  }

  /**
   * Just like {@link #initToSimplifiedInCell(S2Polygon, S2Cell, S1Angle, S1Angle)} above, but uses
   * a default boundaryTolerance of 1e-15 radians.
   */
  @JsIgnore
  public void initToSimplifiedInCell(S2Polygon a, final S2Cell cell, S1Angle snapRadius) {
    initToSimplifiedInCell(a, cell, snapRadius, S1Angle.radians(1e-15));
  }

  /** See comments in {@link #initToSimplifiedInCell}. */
  private List<S2Polyline> simplifyEdgesInCell(
      S2Polygon a, S2Cell cell, double toleranceUV, S1Angle snapRadius) {
    S2Builder.Builder options = new S2Builder.Builder(new IdentitySnapFunction(snapRadius));
    options.setSimplifyEdgeChains(true);
    S2Builder builder = options.build();

    // The output consists of a sequence of polylines. Polylines consisting of interior edges are
    // simplified using S2Builder, while polylines consisting of boundary edges are returned
    // unchanged.
    List<S2Polyline> outputPolylines = new ArrayList<>();
    List<S2PolylineLayer> interiorPolylines = new ArrayList<>();
    for (int i = 0; i < a.numLoops(); ++i) {
      S2Loop aLoop = a.loop(i);
      S2Point v0 = aLoop.orientedVertex(0);
      byte mask0 = getCellEdgeIncidenceMask(cell, v0, toleranceUV);
      boolean inInterior = false; // Was the last edge an interior edge?
      for (int j = 1; j <= aLoop.numVertices(); ++j) {
        S2Point v1 = aLoop.orientedVertex(j);
        byte mask1 = getCellEdgeIncidenceMask(cell, v1, toleranceUV);
        if ((mask0 & mask1) != 0) {
          // This is an edge along the cell boundary. Such edges do not get simplified; we add them
          // directly to the output. (We create a separate polyline for each edge to keep things
          // simple.)  We call forceVertex on all boundary vertices to ensure that they don't move,
          // and so that nearby interior edges are snapped to them.
          assert !inInterior;
          builder.forceVertex(v1);
          outputPolylines.add(new S2Polyline(new S2Point[] {v0, v1}));
        } else {
          // This is an interior edge. If this is the first edge of an interior chain, then start a
          // new S2Builder layer. Also ensure that any polyline vertices on the boundary do not
          // move, so that they will still connect with any boundary edge(s) there.
          if (!inInterior) {
            S2PolylineLayer layer = new S2PolylineLayer();
            builder.startLayer(layer);
            interiorPolylines.add(layer);
            inInterior = true;
          }
          builder.addEdge(v0, v1);
          if (mask1 != 0) {
            builder.forceVertex(v1);
            inInterior = false; // Terminate this polyline.
          }
        }
        v0 = v1;
        mask0 = mask1;
      }
    }

    // Process the interior edges and gather the resulting polylines.
    S2Error error = new S2Error();
    boolean unused = builder.build(error);
    assert error.ok() : "InitToSimplifiedInCell failed: " + error.text();

    for (S2PolylineLayer layer : interiorPolylines) {
      outputPolylines.add(layer.getPolyline());
    }

    return outputPolylines;
  }

  /**
   * Given a point "p" inside an S2Cell or on its boundary, return a mask indicating which of the
   * S2Cell edges the point lies on. All boundary comparisons are to within a maximum "u" or "v"
   * error of "toleranceUV". Bit "i" in the result is set if and only "p" is incident to the edge
   * corresponding to S2Cell.edge(i).
   */
  byte getCellEdgeIncidenceMask(S2Cell cell, S2Point p, double toleranceUV) {
    byte mask = 0;
    @Nullable R2Vector uv = S2Projections.faceXyzToUv(cell.face(), p);
    if (uv != null) {
      R2Rect bound = cell.getBoundUV();
      assert bound.expanded(toleranceUV).contains(uv);
      if (abs(uv.get(1) - bound.y().lo()) <= toleranceUV) {
        mask |= 1;
      }
      if (abs(uv.get(0) - bound.x().hi()) <= toleranceUV) {
        mask |= 2;
      }
      if (abs(uv.get(1) - bound.y().hi()) <= toleranceUV) {
        mask |= 4;
      }
      if (abs(uv.get(0) - bound.x().lo()) <= toleranceUV) {
        mask |= 8;
      }
    }
    return mask;
  }

  /**
   * Takes a set of possibly intersecting edges, stored in the S2ShapeIndex, and breaks the edges
   * into small pieces so that there is no intersection anymore, and adds all these edges to the
   * builder.
   *
   * @deprecated Instead of using the deprecated S2PolygonBuilder, use S2Builder, and set the
   *     splitCrossingEdges option.
   */
  @Deprecated
  public static void breakEdgesAndAddToBuilder(S2ShapeIndex index, S2PolygonBuilder builder) {
    // If there are self intersections, we add the pieces separately.
    // addSharedEdges ("true" below) can be false or true: it makes no difference due to the way we
    // call ClipEdge.
    EdgeClipper clipper = new EdgeClipper(index, true, S2Polygon::reverseNone);
    List<ParametrizedS2Point> intersections = Lists.newArrayList();
    MutableEdge edge = new MutableEdge();
    for (S2Shape shape : index.getShapes()) {
      int numEdges = shape.numEdges();
      for (int e = 0; e < numEdges; e++) {
        shape.getEdge(e, edge);
        if (!edge.a.equalsPoint(edge.b)) {
          intersections.add(new ParametrizedS2Point(0, edge.a));
          clipper.clipEdge(edge.a, edge.b, intersections);
          intersections.add(new ParametrizedS2Point(1, edge.b));
          Collections.sort(intersections);
          for (int k = 0; k + 1 < intersections.size(); ++k) {
            S2Point p1 = intersections.get(k).getPoint();
            S2Point p2 = intersections.get(k + 1).getPoint();
            if (!p1.equalsPoint(p2)) {
              builder.addEdge(p1, p2);
            }
          }
          intersections.clear();
        }
      }
    }
  }

  /**
   * Returns a polygon that is the union of the given polygons.
   *
   * <p>This is a convenience wrapper for the union method below which provides a default
   * snapFunction, with a snap radius of {@link S2EdgeUtil#INTERSECTION_MERGE_RADIUS}, (equal to
   * about 1.8e-15 radians or 11 nanometers on the Earth's surface). This means that vertices may be
   * positioned arbitrarily, but vertices that are extremely close together can be
   * merged together. The reason for a non-zero default snap radius is that it helps to eliminate
   * narrow cracks and slivers when T-vertices are present. For example, adjacent S2Cells at
   * different levels do not share exactly the same boundary, so there can be a narrow crack between
   * them. If a polygon is intersected with those cells and the pieces are unioned together, the
   * result would have a narrow crack unless the snap radius is set to a non-zero value.
   */
  public static S2Polygon union(Iterable<S2Polygon> polygons) {
    return union(polygons, new IdentitySnapFunction(S2EdgeUtil.INTERSECTION_MERGE_RADIUS));
  }

  /**
   * Returns a polygon that is the union of the given polygons. The current implementation
   * repeatedly unions pairs of polygons and snaps each intermediate result with the given
   * snapFunction, using {@link #initToUnion(S2Polygon, S2Polygon, SnapFunction)}.
   *
   * @throws AssertionError if the operation fails or an input polygon is invalid, and Java
   * assertions are enabled.
   */
  @JsMethod(name = "unionSnapped")
  public static S2Polygon union(Iterable<S2Polygon> polygons, SnapFunction snapFunction) {
    // Create a priority queue of polygons in order of number of vertices.
    TreeMultimap<Integer, S2Polygon> queue = TreeMultimap.create();
    for (S2Polygon polygon : polygons) {
      queue.put(polygon.getNumVertices(), polygon);
    }

    // Repeatedly union the two smallest polygons and add the result to the queue until we have a
    // single polygon to return.
    Set<Map.Entry<Integer, S2Polygon>> queueSet = queue.entries();
    while (queueSet.size() > 1) {
      // Pop the two simplest polygons from queue.
      queueSet = queue.entries();
      Iterator<Map.Entry<Integer, S2Polygon>> smallestIter = queueSet.iterator();

      Map.Entry<Integer, S2Polygon> smallest = smallestIter.next();
      S2Polygon aPolygon = smallest.getValue();
      smallestIter.remove();

      smallest = smallestIter.next();
      S2Polygon bPolygon = smallest.getValue();
      smallestIter.remove();

      // Union and add result back to queue.
      S2Polygon unionPolygon = new S2Polygon();
      boolean unused = unionPolygon.initToUnion(aPolygon, bPolygon, snapFunction);

      queue.put(unionPolygon.getNumVertices(), unionPolygon);
    }

    if (queue.isEmpty()) {
      return new S2Polygon();
    } else {
      return queue.get(queue.asMap().firstKey()).first();
    }
  }

  /**
   * Returns a polygon that is the union of the given polygons; combines vertices that form edges
   * that are within {@link S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE} of each other.
   *
   * @deprecated Use {@link #union(Iterable<S2Polygon>)} instead. Note that the new implementation
   * is built on S2Builder and S2BooleanOperation, instead of S2PolygonBuilder, and the results are
   * not identical.
   */
  @JsIgnore
  @Deprecated
  public static S2Polygon unionSloppy(Iterable<S2Polygon> polygons) {
    return unionSloppy(polygons, S2EdgeUtil.DEFAULT_INTERSECTION_TOLERANCE);
  }

  /**
   * Returns a polygon that is the union of the given polygons; combines vertices that form edges
   * that are almost identical, as defined by {@code vertexMergeRadius}.
   *
   * @deprecated Use {@link #union(Iterable<S2Polygon>, SnapFunction)} instead. Note that the new
   *     implementation is built on S2Builder and S2BooleanOperation, instead of S2PolygonBuilder,
   *     and the results are not identical.
   */
  @Deprecated
  public static S2Polygon unionSloppy(Iterable<S2Polygon> polygons, S1Angle vertexMergeRadius) {
    // Effectively create a priority queue of polygons in order of number of vertices. Repeatedly
    // union the two smallest polygons and add the result to the queue until we have a single
    // polygon to return.

    // map: # of vertices -> polygon
    TreeMultimap<Integer, S2Polygon> queue = TreeMultimap.create();

    for (S2Polygon polygon : polygons) {
      queue.put(polygon.getNumVertices(), polygon);
    }

    Set<Map.Entry<Integer, S2Polygon>> queueSet = queue.entries();
    while (queueSet.size() > 1) {
      // Pop two simplest polygons from queue.
      queueSet = queue.entries();
      Iterator<Map.Entry<Integer, S2Polygon>> smallestIter = queueSet.iterator();

      Map.Entry<Integer, S2Polygon> smallest = smallestIter.next();
      S2Polygon aPolygon = smallest.getValue();
      smallestIter.remove();

      smallest = smallestIter.next();
      S2Polygon bPolygon = smallest.getValue();
      smallestIter.remove();

      // Union and add result back to queue.
      S2Polygon unionPolygon = new S2Polygon();
      unionPolygon.initToUnionSloppy(aPolygon, bPolygon, vertexMergeRadius);
      queue.put(unionPolygon.getNumVertices(), unionPolygon);
    }

    if (queue.isEmpty()) {
      return new S2Polygon();
    } else {
      return queue.get(queue.asMap().firstKey()).first();
    }
  }

  /**
   * Returns true if this polygon contains the given polyline. Note that a polygon contains its
   * boundary edges, but not reversed boundary edges.
   */
  @JsMethod(name = "containsPolyline")
  public boolean contains(S2Polyline b) {
    return subtractFromPolyline(b).isEmpty();
    // TODO(user): This would be much faster if implemented as S2BooleanOperation.isEmpty(),
    // similar to the "disjoint" method below, but with OpType.DIFFERENCE. However, the
    // S2BooleanOperation options need to be set correctly to maintain the existing, desired
    // boundary containment behavior. With default options, boundary edges are not contained.
  }

  /**
   * Returns true if this polygon and the given polyline are disjoint, i.e. have no points in
   * common.
   */
  public boolean disjoint(S2Polyline b) {
    // TODO(user): Add this method to the C++ implementation.
    S2ShapeIndex polylineIndex = new S2ShapeIndex();
    polylineIndex.add(b.shape());
    return S2BooleanOperation.isEmpty(
        S2BooleanOperation.OpType.INTERSECTION, index(), polylineIndex);
  }

  /**
   * Returns true if this polygon intersects the given polyline, i.e. have at least one point in
   * common.
   */
  public boolean intersectsPolyline(S2Polyline b) {
    // TODO(user): Add this method to the C++ implementation.
    return !disjoint(b);
  }

  /**
   * Performs the given S2BooleanOperation.OpType on this polygon and the given S2Polyline. Snaps
   * the resulting edges with the given SnapFunction and sends them to an S2PolylineVectorLayer,
   * which assembles the edges into polylines. The polylines are returned without being validated.
   * Note that degenerate polylines are discarded.
   */
  public List<S2Polyline> operationWithPolyline(
      S2BooleanOperation.OpType opType, SnapFunction snapFunction, S2Polyline a) {
    S2PolylineVectorLayer.Options layerOptions = new S2PolylineVectorLayer.Options();
    layerOptions.setPolylineType(PolylineType.WALK);
    S2PolylineVectorLayer layer = new S2PolylineVectorLayer(layerOptions);

    S2ShapeIndex lineIndex = new S2ShapeIndex();
    lineIndex.add(a.shape());

    S2BooleanOperation op = new S2BooleanOperation.Builder(snapFunction).build(opType, layer);
    S2Error error = new S2Error();
    boolean built = op.build(lineIndex, index(), error);
    if (built && error.ok()) {
      return layer.getPolylines();
    }

    // Failures should not happen, as S2BooleanOperation.build() does not set any errors, nor does
    // S2PolylineVectorLayer when not validating results.
    throw new IllegalStateException(
        "OperationWithPolyline failed for " + opType + ": " + error.text());
  }

  /**
   * Intersects this polygon with the polyline "in", and then assembles the results into zero or
   * more polylines. The returned polylines are ordered in the order they would be encountered by
   * traversing "in" from beginning to end.
   *
   * <p>Note that degenerate polylines are discarded, so it is possible to have an S2Polyline and
   * S2Polygon where polygon.intersects(polyline) is true, but the result of
   * polygon.intersectWithPolyline(polyline) is empty.
   */
  public List<S2Polyline> intersectWithPolyline(S2Polyline in) {
    return intersectWithPolyline(
        in, new IdentitySnapFunction(S2EdgeUtil.INTERSECTION_MERGE_RADIUS));
  }

  /**
   * Intersects this polygon with the polyline "in", snaps the resulting edges with the given
   * SnapFunction, and then assembles the results into zero or more polylines. The returned
   * polylines are ordered in the order they would be encountered by traversing "in" from beginning
   * to end.
   *
   * <p>Note that degenerate polylines are discarded, so it is possible to have an S2Polyline and
   * S2Polygon where polygon.intersects(polyline) is true, but the result of
   * polygon.intersectWithPolyline(polyline) is empty.
   */
  @JsMethod(name = "intersectWithPolylineSnapped")
  public List<S2Polyline> intersectWithPolyline(S2Polyline in, SnapFunction snapFunction) {
    return operationWithPolyline(S2BooleanOperation.OpType.INTERSECTION, snapFunction, in);
  }

  /**
   * Similar to {@link #intersectWithPolyline}, except that vertices will be dropped as necessary to
   * ensure that all adjacent vertices in the sequence obtained by concatenating the output
   * polylines will be farther than {@code vertexMergeRadius} apart. Note that this can change the
   * number of output polylines and/or yield single-vertex polylines.
   *
   * @deprecated Use intersectWithPolyline(in, new IdentitySnapFunction(vertexMergeRadius)) but note
   *     it has slightly different behavior.
   */
  @Deprecated
  public List<S2Polyline> intersectWithPolylineSloppy(S2Polyline in, S1Angle vertexMergeRadius) {
    return internalClipPolyline(false, in, vertexMergeRadius);
  }

  /**
   * Subtracts this polygon from the given {@link S2Polyline}, and snaps the resulting edges
   * using an IdentitySnapFunction with the minimum merge radius. Assembles the results into zero or
   * more polylines.
   *
   * <p>The polylines are ordered in the order they would be encountered by traversing the
   * given polyline {@code in} from beginning to end. Note that degenerate polylines are discarded.
   */
  public List<S2Polyline> subtractFromPolyline(S2Polyline in) {
    return subtractFromPolyline(in, S2EdgeUtil.INTERSECTION_MERGE_RADIUS);
  }

  /**
   * Subtracts this polygon from the given {@link S2Polyline}, and snaps the resulting edges using
   * an IdentitySnapFunction with the given snapRadius. Assembles the results into zero or more
   * polylines.
   *
   * <p>The returned polylines are ordered in the order they would be encountered by traversing the
   * given polyline {@code in} from beginning to end. Note that degenerate polylines are discarded.
   * (This method is named ApproxSubtractFromPolyline in the C++ implementation.)
   */
  @JsMethod(name = "approxSubtractFromPolyline")
  public List<S2Polyline> subtractFromPolyline(S2Polyline in, S1Angle snapRadius) {
    return subtractFromPolyline(in, new IdentitySnapFunction(snapRadius));
  }

  /**
   * Subtracts this polygon from the given {@link S2Polyline}, and snaps the resulting edges using
   * the given SnapFunction. Assembles the results into zero or more polylines.
   *
   * <p>The returned polylines are ordered in the order they would be encountered by traversing the
   * given polyline {@code in} from beginning to end. Note that degenerate polylines are discarded.
   */
  @JsMethod(name = "subtractFromPolylineSnapped")
  public List<S2Polyline> subtractFromPolyline(S2Polyline in, S2Builder.SnapFunction snapFunction) {
    return operationWithPolyline(S2BooleanOperation.OpType.DIFFERENCE, snapFunction, in);
  }

  /**
   * Similar to {@link #intersectWithPolylineSloppy}, but subtracts this polygon from the given
   * polyline.
   *
   * @deprecated Use {@link #subtractFromPolyline(S2Polyline)} instead.
   */
  @Deprecated
  public List<S2Polyline> subtractFromPolylineSloppy(S2Polyline in, S1Angle vertexMergeRadius) {
    return internalClipPolyline(true, in, vertexMergeRadius);
  }

  /**
   * Clips the {@link S2Polyline} {@code a} to the interior of this polygon. The resulting
   * polyline(s) will be returned. If {@code invert} is {@code true}, we clip {@code a} to the
   * exterior of this polygon instead. Vertices will be dropped such that adjacent vertices will not
   * be closer than {@code mergeRadius}.
   *
   * <p>We do the intersection/subtraction by walking the polyline edges. For each edge, we compute
   * all intersections with the polygon boundary and sort them in increasing order of distance along
   * that edge. We then divide the intersection points into pairs, and output a clipped polyline
   * segment for each one. We keep track of whether we're inside or outside of the polygon at all
   * times to decide which segments to output.
   *
   * @deprecated Used only by the deprecated methods subtractFromPolylineSloppy and
   * intersectWithPolylineSloppy.
   */
  @Deprecated
  private List<S2Polyline> internalClipPolyline(boolean invert, S2Polyline a, S1Angle mergeRadius) {
    List<S2Polyline> out = Lists.newArrayList();
    EdgeClipper clipper = new EdgeClipper(index(), true, S2Polygon::reverseNone);
    List<ParametrizedS2Point> intersections = Lists.newArrayList();
    List<S2Point> vertices = Lists.newArrayList();
    int n = a.numVertices();
    boolean inside = contains(a.vertex(0)) ^ invert;
    for (int j = 0; j < n - 1; j++) {
      S2Point a0 = a.vertex(j);
      S2Point a1 = a.vertex(j + 1);
      clipper.clipEdge(a0, a1, intersections);
      if (inside) {
        intersections.add(new ParametrizedS2Point(0, a0));
      }
      inside = (intersections.size() & 1) != 0;
      // assert ((contains(a1) ^ invert) == inside);
      if (inside) {
        intersections.add(new ParametrizedS2Point(1, a1));
      }
      Collections.sort(intersections);
      // At this point we have a sorted array of vertex pairs representing the edge(s) obtained
      // after clipping (a0, a1) against the polygon.
      for (int k = 0; k < intersections.size(); k += 2) {
        S2Point v0 = intersections.get(k).getPoint();
        S2Point v1 = intersections.get(k + 1).getPoint();
        if (!v0.equalsPoint(v1)) {
          // If the gap from the previous vertex to this one is large enough, start a new polyline.
          if (!vertices.isEmpty()
              && vertices.get(vertices.size() - 1).angle(v0) > mergeRadius.radians()) {
            out.add(new S2Polyline(vertices));
            vertices.clear();
          }
          // Append this segment to the current polyline, ignoring any vertices that are too close
          // to the previous vertex.
          if (vertices.isEmpty()) {
            vertices.add(v0);
          }
          if (vertices.get(vertices.size() - 1).angle(v1) > mergeRadius.radians()) {
            vertices.add(v1);
          }
        }
      }
      intersections.clear();
    }
    if (!vertices.isEmpty()) {
      out.add(new S2Polyline(vertices));
    }
    return out;
  }

  /**
   * Return true if every loop of this polygon shares at most one vertex with its parent loop. Every
   * polygon has a unique normalized form. A polygon can be normalized by passing it through
   * S2Builder (with no snapping) in order to reconstruct the polygon from its edges.
   *
   * <p>Generally there is no reason to convert polygons to normalized form. It is mainly useful for
   * testing in order to compare whether two polygons have exactly the same interior, even when they
   * have a different loop structure. For example, a diamond nested within a square (touching at
   * four points) could be represented as a square with a diamond-shaped hole, or as four triangles.
   * Methods such as {@link #boundaryApproxEquals(S2Polygon, double)} will report these polygons as
   * being different (because they have different boundaries) even though they contain the same
   * points. However if they are both converted to normalized form (the "four triangles" version)
   * then they can be compared more easily.
   */
  @SuppressWarnings("ReferenceEquality") // See comment above "parent != lastParent".
  public boolean isNormalized() {
    // TODO(user) per ericv@: The condition tested here is insufficient. The correct
    // condition is that each *connected component* of child loops can share at most one vertex
    // with their parent loop. Example: suppose loop A has children B, C, D, and the following pairs
    // are connected: AB, BC, CD, DA. Then the polygon is not normalized.
    Set<S2Point> vertices = Sets.newHashSet();
    S2Loop lastParent = null;
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop child = loop(i);
      if (child.depth() == 0) {
        continue;
      }
      S2Loop parent = loop(getParent(i));
      // Test if the loops are different; we can use identity since any two loops in a valid polygon
      // can only be equal by value if they're equal by reference.
      if (parent != lastParent) {
        vertices.clear();
        for (int j = 0; j < parent.numVertices(); ++j) {
          vertices.add(parent.vertex(j));
        }
        lastParent = parent;
      }
      int count = 0;
      for (int j = 0; j < child.numVertices(); ++j) {
        if (vertices.contains(child.vertex(j))) {
          ++count;
        }
      }
      if (count > 1) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns true if two polygons have the same boundary. More precisely, this method requires that
   * both polygons have loops with the same cyclic vertex order and the same nesting hierarchy.
   * (This implies that vertices may be cyclically rotated between corresponding loops, and the loop
   * ordering may be different between the two polygons as long as the nesting hierarchy is the
   * same.)
   */
  @VisibleForTesting
  boolean boundaryEquals(S2Polygon b) {
    if (numLoops() != b.numLoops()) {
      return false;
    }

    for (int i = 0; i < numLoops(); ++i) {
      S2Loop aLoop = loop(i);
      boolean success = false;
      for (int j = 0; j < numLoops(); ++j) {
        S2Loop bLoop = b.loop(j);
        if ((bLoop.depth() == aLoop.depth()) && bLoop.boundaryEquals(aLoop)) {
          success = true;
          break;
        }
      }
      if (!success) {
        return false;
      }
    }
    return true;
  }

  /**
   * Return true if two polygons are approximately equal to within the given tolerance. This is true
   * if it is possible to move the vertices of the two polygons so that they contain the same set of
   * points.
   *
   * <p>Note that according to this model, small regions less than "tolerance" in width do not need
   * to be considered, since these regions can be collapsed into degenerate loops (which contain no
   * points) by moving their vertices.
   *
   * <p>This model is not as strict as using the Hausdorff distance would be, and it is also not as
   * strict as {@link #boundaryNear(S2Polygon, double)} (defined below). However, it is a good
   * choice for comparing polygons that have been snapped, simplified, unioned, etc, since these
   * operations use a model similar to this one (i.e., degenerate loops or portions of loops are
   * automatically removed).
   */
  public boolean approxEquals(S2Polygon b, S1Angle tolerance) {
    // TODO(user): This can be implemented more cheaply with S2Builder, by simply adding all
    // the edges from one polygon, adding the reversed edge from the other polygon, and turning on
    // the options to split edges and discard sibling pairs. Then the polygons are approximately
    // equal if the output graph has no edges.
    S2Polygon symmetricDifference = new S2Polygon();
    symmetricDifference.initToSymmetricDifference(b, this, new IdentitySnapFunction(tolerance));
    return symmetricDifference.isEmpty();
  }

  /**
   * Returns true if two polygons have the same boundary, except for vertex perturbations. Both
   * polygons must have loops with the same cyclic vertex order and the same nesting hierarchy, but
   * the vertex locations are allowed to differ by up to {@code maxErrorRadians} radians. Note: This
   * method is intended for testing purposes.
   */
  @VisibleForTesting
  public boolean boundaryApproxEquals(S2Polygon b, double maxErrorRadians) {
    if (numLoops() != b.numLoops()) {
      return false;
    }

    // For now, we assume that there is at most one candidate match for each loop. (So far this
    // method is just used for testing.)
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop aLoop = loop(i);
      boolean success = false;
      for (int j = 0; j < numLoops(); ++j) {
        S2Loop bLoop = b.loop(j);
        if (bLoop.depth() == aLoop.depth() && bLoop.boundaryApproxEquals(aLoop, maxErrorRadians)) {
          success = true;
          break;
        }
      }
      if (!success) {
        return false;
      }
    }
    return true;
  }

  /**
   * Returns true if two polygons have boundaries that are within {@code maxErrorRadians} of each
   * other along their entire lengths. More precisely, there must be a bijection between the two
   * sets of loops such that {@code aLoop.boundaryNear(bLoop)} is true for each aLoop in {@code
   * this} and each bLoop in {@code b}.
   */
  boolean boundaryNear(S2Polygon b, double maxErrorRadians) {
    if (numLoops() != b.numLoops()) {
      return false;
    }

    // For now, we assume that there is at most one candidate match for each loop.
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop aLoop = loop(i);
      boolean success = false;
      for (int j = 0; j < numLoops(); ++j) {
        S2Loop bLoop = b.loop(j);
        if (bLoop.depth() == aLoop.depth() && bLoop.boundaryNear(aLoop, maxErrorRadians)) {
          success = true;
          break;
        }
      }
      if (!success) {
        return false;
      }
    }
    return true;
  }

  // S2Region interface (see S2Region.java for details):

  /**
   * Returns a spherical cap that bounds this loop. It may be expanded slightly such that if the
   * loop contains a point P, then the bound contains P also.
   */
  @Override
  public S2Cap getCapBound() {
    return bound.getCapBound();
  }

  /**
   * Returns a fairly tight bounding latitude-longitude rectangle. It is not guaranteed to be as
   * tight as possible, to ensure that if the loop contains a point P, then the bound contains P
   * also.
   */
  @Override
  public S2LatLngRect getRectBound() {
    return bound;
  }

  @Override
  public void getCellUnionBound(Collection<S2CellId> results) {
    new S2ShapeIndexRegion(index).getCellUnionBound(results);
  }

  /**
   * If this method returns true, the region completely contains the given cell. Otherwise, either
   * the region does not contain the cell or the containment relationship could not be determined.
   */
  @Override
  @JsIgnore
  public boolean contains(S2Cell cell) {
    return new S2ShapeIndexRegion(index()).contains(cell);
  }

  /**
   * If this method returns false, the region does not intersect the given cell. Otherwise, either
   * region intersects the cell, or the intersection relationship could not be determined.
   */
  @Override
  public boolean mayIntersect(S2Cell target) {
    return new S2ShapeIndexRegion(index()).mayIntersect(target);
  }

  /**
   * Returns true if this polygon contains the given point. Note that not all vertices of a polygon
   * boundary are necessarily contained by that polygon. A shared point, or shared (but reversed)
   * edge between two "touching" polygons is contained by just one of those polygons.
   *
   * <p>The point {@code p} must be normalized (unit length). Note that this changed in late 2024.
   */
  @Override
  @JsMethod(name = "containsPoint")
  public boolean contains(S2Point p) {
    assert S2.isUnitLength(p);
    // Note(ericv): A bounds check slows down the C++ version of this method by about 50%. It is
    // worthwhile only when it might allow us to delay building the index.
    if (!index.isFresh() && !bound.contains(p)) {
      return false;
    }

    // For small polygons it is faster to just check all the crossings. Otherwise we keep track of
    // the number of calls to contains() and only build the index once enough calls have been made
    // so that we think it is worth the effort. See S2Loop.contains(S2Point) for detailed comments.
    int maxBruteForceVertices = 32;
    if (getNumVertices() <= maxBruteForceVertices
        || (!index.isFresh() && unindexedContainsCalls.decrementAndGet() > 0)) {
      boolean inside = false;
      for (int i = 0; i < numLoops(); ++i) {
        // Use brute force to avoid building the loops' S2ShapeIndex.
        boolean loopContainsP = loop(i).bruteForceContains(p);
        // Set 'inside' to the bitwise XOR of 'inside' and 'loopContainsP'.
        inside ^= loopContainsP;
      }
      return inside;
    }

    // Otherwise we use a ContainsPointQuery to determine point containment.
    S2ContainsPointQuery query = new S2ContainsPointQuery(index);
    return query.contains(p);
  }

  /** For each map entry, sorts the value list. */
  private static void sortValueLoops(Map<S2Loop, List<S2Loop>> loopMap) {
    for (List<S2Loop> value : loopMap.values()) {
      Collections.sort(value);
    }
  }

  /** Add the given "newLoop" to the loopMap to the deepest loop containing it. */
  private static void insertLoop(S2Loop newLoop, Map<S2Loop, List<S2Loop>> loopMap) {
    // Insert a new loop list to hold the children of newLoop (if any).
    List<S2Loop> newChildren = loopMap.computeIfAbsent(newLoop, k -> new ArrayList<S2Loop>());
    S2Loop parent = ROOT;

    // Find the deepest-nested loop containing "newLoop". That parent's "children" list is where
    // newLoop will be added.
    List<S2Loop> children = null;
    for (boolean done = false; !done; ) {
      children = loopMap.get(parent);
      done = true;
      for (S2Loop child : children) {
        if (child.containsNested(newLoop)) {
          parent = child;
          done = false;
          break;
        }
      }
    }

    // Some of the other children of the parent loop may now be children of the new loop.
    for (int i = 0; i < children.size(); ) {
      S2Loop child = children.get(i);
      if (newLoop.containsNested(child)) {
        newChildren.add(child);
        children.remove(i);
      } else {
        ++i;
      }
    }
    children.add(newLoop);
  }

  /**
   * If the given loop is not {@code ROOT}, it must be in the given loopMap, and is assigned the
   * given depth. Then the children of the loop are recursively initialized with depth+1. If the
   * given loop is {@code ROOT}, then the top-level loops in the loopMap (children of the ROOT loop)
   * are initialized recursively.
   */
  @SuppressWarnings("ReferenceEquality") // ROOT is a specific instance.
  private void initLoop(S2Loop loop, int depth, IdentityHashMap<S2Loop, List<S2Loop>> loopMap) {
    Preconditions.checkNotNull(loop); // Can't fail. Existing initLoop() callers check this.
    if (loop != ROOT) {
      loop.setDepth(depth);
      loops.add(loop);
    }
    List<S2Loop> children = loopMap.get(loop);
    if (children != null) {
      for (S2Loop child : children) {
        initLoop(child, depth + 1, loopMap);
      }
    }
  }

  /**
   * Returns +1 if this polygon (A) contains the boundary of B, -1 if A excludes the boundary of B,
   * and 0 if the boundaries of A and B cross.
   */
  int compareBoundary(S2Loop b) {
    int result = -1;
    for (S2Loop loop : loops) {
      // If B crosses any loop of A, the result is 0. Otherwise the result changes sign each time B
      // is contained by a loop of A.
      result *= -loop.compareBoundary(b);
      if (result == 0) {
        break;
      }
    }
    return result;
  }

  /** Returns true if this polygon (A) contains the entire boundary of B. */
  private boolean containsBoundary(S2Polygon b) {
    for (S2Loop loop : b.loops) {
      if (compareBoundary(loop) <= 0) {
        return false;
      }
    }
    return true;
  }

  /** Returns true if this polygon (A) excludes the entire boundary of B. */
  private boolean excludesBoundary(S2Polygon b) {
    for (S2Loop loop : b.loops) {
      if (compareBoundary(loop) >= 0) {
        return false;
      }
    }
    return true;
  }

  /**
   * Given a loop B whose boundaries do not cross any loop of this polygon, returns true if this
   * polygon contains the boundary of B. Shared edges are handled according to the rule described in
   * {@link S2Loop#containsNonCrossingBoundary(int, boolean)}.
   */
  private boolean containsNonCrossingBoundary(S2Loop b, boolean bReverse) {
    boolean inside = false;
    for (int i = 0; i < loops.size(); i++) {
      inside ^= loop(i).containsNonCrossingBoundary(b, bReverse);
    }
    return inside;
  }

  /**
   * Given a polygon B such that the boundary of this polygon does not cross any loop of B, returns
   * true if this polygon excludes all shell boundaries of B.
   */
  private boolean excludesNonCrossingShells(S2Polygon b) {
    for (S2Loop bLoop : b.loops) {
      if (!bLoop.isHole()) {
        if (containsNonCrossingBoundary(bLoop, false)) {
          return false;
        }
      }
    }
    return true;
  }

  /**
   * Given a polygon B such that the boundary of this polygon does not cross any loop of B, returns
   * true if this polygon excludes all shell boundaries of the complement of B.
   */
  private boolean excludesNonCrossingComplementShells(S2Polygon b) {
    // Special case to handle the complement of the empty or full polygons.
    if (b.isEmpty()) {
      return !isFull();
    }

    if (b.isFull()) {
      return true;
    }

    // Otherwise the complement of B may be obtained by inverting loop(0) and then swapping the
    // shell/hole status of all other loops. This implies that the shells of the complement consist
    // of loop 0 plus all the holes of the original polygon.
    for (int j = 0; j < b.numLoops(); ++j) {
      if (j > 0 && !b.loop(j).isHole()) {
        continue;
      }

      // The interior of the complement is to the right of loop 0, and to the left of the loops that
      // were originally holes.
      if (containsNonCrossingBoundary(b.loop(j), j == 0)) {
        return false;
      }
    }

    return true;
  }

  /** Returns true if any loop intersects the given loop. */
  private boolean anyLoopIntersects(S2Loop b) {
    for (S2Loop loop : loops) {
      if (loop.intersects(b)) {
        return true;
      }
    }
    return false;
  }

  /** Returns a human readable representation of the polygon. */
  @Override
  public String toString() {
    StringBuilder sb = new StringBuilder();
    sb.append("Polygon: (").append(numLoops()).append(") loops:\n");
    for (int i = 0; i < numLoops(); ++i) {
      S2Loop s2Loop = loop(i);
      sb.append("loop with depth ").append(s2Loop.depth()).append(" <\n");
      for (int v = 0; v < s2Loop.numVertices(); ++v) {
        S2Point s2Point = s2Loop.vertex(v);
        sb.append(s2Point.toDegreesString());
        // end of vertex
        sb.append("\n");
      }
      // end of loop
      sb.append(">\n");
    }
    return sb.toString();
  }

  /**
   * Encodes the polygon into an efficient, lossless binary representation, which can be decoded by
   * calling {@link S2Polygon#decode}. The encoding is byte-compatible with the C++ version of the
   * S2 library.
   *
   * <p>Benchmarks show that the compressed format takes about 2x the time to decode as the
   * uncompressed format.
   *
   * @param output The output stream into which the encoding should be written.
   * @throws IOException if there was a problem writing into the output stream.
   */
  @JsIgnore // OutputStream is not available to J2CL.
  public void encode(OutputStream output) throws IOException {
    // TODO(torrey): The C++ implementation takes a CodingHint (FAST or COMPACT) and for COMPACT,
    // just uses S2CellId.MAX_LEVEL rather than attempting to find the best snap level.
    int level = getNumVertices() == 0 ? S2CellId.MAX_LEVEL : getBestSnapLevel();
    LittleEndianOutput encoder = new LittleEndianOutput(output);
    if (level == -1) {
      encodeUncompressed(encoder);
    } else {
      encodeCompressed(level, encoder);
    }
  }

  /**
   * Encodes the polygon into an uncompressed binary representation, which can be decoded by calling
   * {@link S2Polygon#decode(InputStream)}. The encoding is byte-compatible with the C++ version of
   * the S2 library.
   *
   * @param encoder The output stream into which the encoding should be written.
   * @throws IOException if there was a problem writing into the output stream.
   */
  @VisibleForTesting // Visible for benchmarking.
  @JsIgnore // LittleEndianOutput is not available to J2CL.
  public void encodeUncompressed(LittleEndianOutput encoder) throws IOException {
    encoder.writeByte(LOSSLESS_ENCODING_VERSION);
    // Placeholder value for backward compatibility; previously stored the "owns_loops_" value.
    encoder.writeByte((byte) 1);
    // Encode obsolete "hasHoles" field for backwards compatibility.
    encoder.writeByte((byte) (hasHoles() ? 1 : 0));
    encoder.writeInt(loops.size());
    for (S2Loop loop : loops) {
      loop.encode(encoder);
    }
    bound.encode(encoder);
  }

  /** Returns true if any of the loops of this polygon are holes. */
  private boolean hasHoles() {
    boolean hasHoles = false;
    for (S2Loop loop : loops) {
      if (loop.isHole()) {
        hasHoles = true;
      }
    }
    return hasHoles;
  }

  private void encodeCompressed(int level, LittleEndianOutput encoder) throws IOException {
    encoder.writeByte(COMPRESSED_ENCODING_VERSION);
    encoder.writeByte((byte) level);
    encoder.writeVarint32(numLoops());
    for (int i = 0; i < numLoops(); i++) {
      loop(i).encodeCompressed(level, encoder);
    }
    // Do not write the bound or numVertices as they can be cheaply recomputed by decodeCompressed.
    // Microbenchmarks show the speed difference is inconsequential.
  }

  /**
   * Decodes a polygon that was encoded using {@link S2Polygon#encode}.
   *
   * <p>This method will never return null. It will either throw an exception or return a valid
   * {@link S2Polygon}.
   *
   * @param input The input stream containing the encoded polygon data.
   * @return the decoded {@link S2Polygon}.
   * @throws IOException if there was a problem reading from the input stream, or the decoded
   *     polygon is invalid.
   */
  @JsIgnore // InputStream is not available to J2CL.
  public static S2Polygon decode(InputStream input) throws IOException {
    LittleEndianInput decoder = new LittleEndianInput(input);
    byte version = decoder.readByte();
    switch (version) {
      case LOSSLESS_ENCODING_VERSION:
        return decodeUncompressed(decoder);
      case COMPRESSED_ENCODING_VERSION:
        return decodeCompressed(decoder);
      default:
        throw new IOException("Unsupported S2Polygon encoding version " + version);
    }
  }

  /**
   * Decodes from the given input stream, stopping at the first byte after the polygon.
   *
   * <p>Note that the decoded polygon is not checked for validity unless Java assertions are
   * enabled, which is normally only the case for tests. Clients may want to call {@link
   * S2Polygon#isValid()} on the result if they are decoding from an untrusted source.
   */
  @VisibleForTesting // Visible for benchmarking.
  @JsIgnore // LittleEndianInput wraps InputStream and is not available to J2CL.
  public static S2Polygon decodeUncompressed(LittleEndianInput decoder) throws IOException {
    S2Polygon result = new S2Polygon();
    result.clearLoops();

    // Ignore obsolete serialized owns_loops_ value.
    byte unused = decoder.readByte();
    // Ignore obsolete serialized has_holes_ value.
    unused = decoder.readByte();

    // Polygons with no loops are explicitly allowed here: a newly created polygon has zero loops
    // and such polygons encode and decode properly.
    int numLoops = decoder.readInt();
    if (numLoops < 0) {
      throw new IOException("Can only decode polygons with up to 2^31 - 1 loops. Got " + numLoops);
    }
    result.numVertices = 0;

    for (int i = 0; i < numLoops; i++) {
      S2Loop loop = S2Loop.decode(decoder);
      if (isCanonicalEmptyLoop(loop)) {
        continue;
      }
      result.numVertices += loop.numVertices();
      result.loops.add(loop);
    }

    result.bound = S2LatLngRect.decode(decoder);
    result.subregionBound = S2EdgeUtil.RectBounder.expandForSubregions(result.bound);
    result.initIndex();

    // TODO(user): Strengthen this assertion-only check to throw an IOException upon
    // decoding an invalid polygon.
    // result.assertValid();
    return result;
  }

  /**
   * Decodes from the given input stream, stopping at the first byte after the polygon.
   *
   * <p>Note that the decoded polygon is not checked for validity unless Java assertions are
   * enabled, which is normally only the case for tests. Clients may want to call {@link
   * S2Polygon#isValid()} on the result if they are decoding from an untrusted source.
   */
  private static S2Polygon decodeCompressed(LittleEndianInput decoder) throws IOException {
    int level = decoder.readByte();
    if (level > S2CellId.MAX_LEVEL || level < 0) {
      throw new IOException("Invalid level " + level);
    }
    // Polygons with no loops are explicitly allowed here: a newly created polygon has zero loops
    // and such polygons encode and decode properly.
    int numLoops = decoder.readVarint32();
    if (numLoops < 0) {
      throw new IOException("Can only decode polygons with up to 2^31 - 1 loops. Got " + numLoops);
    }
    S2Polygon result = new S2Polygon();
    result.clearLoops();

    // Instantiating with the S2Polygon(List<S2Loops>) constructor might cause loops to get
    // reordered, which would break idempotence across platforms. The loops are in exactly the order
    // that we want.
    for (int i = 0; i < numLoops; i++) {
      S2Loop loop = S2Loop.decodeCompressed(level, decoder);
      if (isCanonicalEmptyLoop(loop)) {
        continue;
      }
      result.loops.add(loop);
    }
    // Recompute other properties, like bound and numVertices, and initialize the index.
    result.initLoopProperties();
    // TODO(user): Strengthen this assertion-only check to actually throw an IOException upon
    // decoding an invalid polygon.
    // result.assertValid();
    return result;
  }

  /**
   * EdgeClipper finds all the intersections of a given edge with the edges contained in an
   * S2ShapeIndex. It is used to implement polygon operations such as intersection and union.
   *
   * <p>TODO(user): Remove this class when the last (deprecated) uses of it are removed.
   */
  private static class EdgeClipper {
    private final S2CrossingEdgeQuery query;
    private final boolean addSharedEdges;
    private final Predicate<S2Shape> reverseEdges;
    private final S2ShapeIndex index;

    /**
     * Initialize an EdgeClipper for the given S2ShapeIndex. If the query edge is the same as an
     * index edge (a "shared edge"), then the edge will be included in the output if and only if
     * {@code addSharedEdges} is true. The {@code reverseEdges} function allows the edges of any
     * index shape to be reversed before this test is performed (this is used to reverse the loop
     * orientation of "holes" in certain algorithms).
     */
    public EdgeClipper(
        S2ShapeIndex index, boolean addSharedEdges, Predicate<S2Shape> reverseEdges) {
      this.index = index;
      this.query = new S2CrossingEdgeQuery(index);
      this.addSharedEdges = addSharedEdges;
      this.reverseEdges = reverseEdges;
    }

    /**
     * Finds all points where the polygon B intersects the edge (a0, a1), and add the corresponding
     * parameter values (in the range [0,1]) to {@code intersections}. The result is unsorted.
     */
    public void clipEdge(S2Point a0, S2Point a1, List<ParametrizedS2Point> intersections) {
      Iterator<Edges> candidates = query.getCandidates(a0, a1).iterator();
      if (!candidates.hasNext()) {
        return;
      }

      // Iterate through the candidate loops, and then the candidate edges within each loop.
      S2EdgeUtil.EdgeCrosser crosser = new S2EdgeUtil.EdgeCrosser(a0, a1);
      MutableEdge edge = new MutableEdge();
      do {
        Edges edges = candidates.next();
        final S2Shape bShape = index.getShapes().get(edges.shapeId());
        while (!edges.isEmpty()) {
          int edgeId = edges.nextEdge();
          bShape.getEdge(edgeId, edge);
          int crossing = crosser.robustCrossing(edge.getStart(), edge.getEnd());
          if (crossing >= 0) {
            addIntersection(
                a0, a1, edge.getStart(), edge.getEnd(), bShape, crossing, intersections);
          }
        }
      } while (candidates.hasNext());
    }

    /**
     * Given two edges A and B such that robustCrossing(A, B) >= 0, determines if they intersect and
     * adds any intersection point to {@code intersections}. {@code bShape} is the S2Shape
     * containing edge B, and {@code crossing} is the result of {@code robustCrossing(A, B)}.
     */
    private void addIntersection(
        S2Point a0,
        S2Point a1,
        S2Point b0,
        S2Point b1,
        S2Shape bShape,
        int crossing,
        List<ParametrizedS2Point> intersections) {
      assert crossing >= 0;
      if (crossing > 0) {
        // There is a proper edge crossing.
        S2Point x = S2EdgeUtil.getIntersection(a0, a1, b0, b1);
        double t = S2EdgeUtil.getDistanceFraction(x, a0, a1);
        intersections.add(new ParametrizedS2Point(t, x));
      } else if (S2EdgeUtil.vertexCrossing(a0, a1, b0, b1)) {
        // There is a crossing at one of the vertices. The basic rule is simple: if a0 equals one of
        // the 'b' vertices, the crossing occurs at t = 0; otherwise, it occurs at t = 1.
        //
        // This has the effect that when two reversed edges exist (i.e., a0 == b1 and a1 == b0) and
        // each edge is clipped against the other, neither one is included in the output (which is
        // correct). However, when two shared edges exist (i.e., a0 == b0 and a1 == b1), both are
        // included in the output (which is incorrect). The 'addSharedEdges' flag gives explicit
        // over whether these shared edges are included in the output; if it is false, shared edges
        // are excluded by changing their intersection parameter from 0 to 1. This allows exactly
        // one copy of shared edges to be preserved, by calling clipBoundary() twice with opposite
        // values of the 'addSharedEdges' flag.
        double t = (a0.equalsPoint(b0) || a0.equalsPoint(b1)) ? 0 : 1;
        if (!addSharedEdges && a1.equalsPoint(reverseEdges.apply(bShape) ? b0 : b1)) {
          // Excludes (a0, a1) from the output.
          t = 1;
        }
        intersections.add(new ParametrizedS2Point(t, t == 0 ? a0 : a1));
      }
    }
  }

  /** Returns a shape wrapping this polygon. */
  public Shape shape() {
    if (numLoops() > Shape.MAX_LINEAR_SEARCH_LOOPS) {
      return binarySearchShape();
    } else {
      return linearSearchShape();
    }
  }

  /**
   * Returns an implementation that does a binary search to map an edge to one of a large number of
   * loops.
   */
  Shape binarySearchShape() {
    final int[] cumulativeEdges = new int[numLoops()];
    int numEdges = 0;
    for (int i = 0; i < loops.size(); i++) {
      cumulativeEdges[i] = numEdges;
      numEdges += loops.get(i).numVertices();
    }
    return new Shape(this) {
      private static final long serialVersionUID = 1L;

      @Override
      public void getEdge(int edgeId, MutableEdge result) {
        int start = Arrays.binarySearch(cumulativeEdges, edgeId);
        if (start < 0) {
          // Found a position just beyond the loop we want.
          start = -start - 2;
        }
        S2Loop loop = loop(start);
        edgeId -= cumulativeEdges[start];
        result.set(loop.orientedVertex(edgeId), loop.orientedVertex(edgeId + 1));
      }

      @Override
      public int getChainStart(int chainId) {
        Preconditions.checkElementIndex(chainId, numChains());
        return cumulativeEdges[chainId];
      }
    };
  }

  /**
   * Returns an implementation that does a linear search to map an edge to one of a number of loops.
   *
   * <p>When the number of loops is small, linear search is faster. Most often there is exactly one
   * loop and the getEdge and getChainStart loops execute zero times.
   */
  Shape linearSearchShape() {
    return new Shape(this) {
      private static final long serialVersionUID = 1L;

      @Override
      public void getEdge(int edgeId, MutableEdge result) {
        S2Loop loop = loop(0);
        for (int i = 1; edgeId >= loop.numVertices(); loop = loop(i++)) {
          edgeId -= loop.numVertices();
        }
        result.set(loop.orientedVertex(edgeId), loop.orientedVertex(edgeId + 1));
      }

      @Override
      public int getChainStart(int chainId) {
        Preconditions.checkElementIndex(chainId, numChains());
        int start = 0;
        for (int i = 0; i < chainId; i++) {
          start += loops.get(i).numVertices();
        }
        return start;
      }
    };
  }

  /** A fast {@link S2Coder} for polygons. */
  public static final S2Coder<S2Polygon> FAST_CODER =
      new S2Coder<S2Polygon>() {
        @Override
        public void encode(S2Polygon value, OutputStream output) throws IOException {
          value.encodeUncompressed(new LittleEndianOutput(output));
        }

        @Override
        public S2Polygon decode(Bytes data, Cursor cursor) throws IOException {
          return S2Polygon.decode(data.toInputStream(cursor));
        }

        @Override
        public boolean isLazy() {
          return false;
        }
      };

  /** A slower {@link S2Coder} for polygons that takes up less space. */
  public static final S2Coder<S2Polygon> COMPACT_CODER =
      new S2Coder<S2Polygon>() {
        @Override
        public void encode(S2Polygon value, OutputStream output) throws IOException {
          value.encode(output);
        }

        @Override
        public S2Polygon decode(Bytes data, Cursor cursor) throws IOException {
          return S2Polygon.decode(data.toInputStream(cursor));
        }

        @Override
        public boolean isLazy() {
          return false;
        }
      };

  /**
   * Wrapper class for indexing a polygon via {@link S2ShapeIndex}. This class has several subtypes
   * to store additional data based on the kind of polygon being indexed.
   *
   * <p>Note that unlike S2Polygon, the edges of this shape are directed such that the polygon
   * interior is always on the left.
   */
  @JsType
  public abstract static class Shape implements S2Shape, Serializable {
    private final S2Polygon polygon;

    /**
     * @param polygon the polygon to present as an {@link S2Shape}.
     */
    public Shape(S2Polygon polygon) {
      this.polygon = polygon;
    }

    /** A compact shape coder. */
    public static final S2Coder<S2Polygon.Shape> FAST_CODER =
        S2Polygon.FAST_CODER.delegating(Shape::polygon, S2Polygon::shape);

    /** A compact shape coder. */
    public static final S2Coder<S2Polygon.Shape> COMPACT_CODER =
        S2Polygon.COMPACT_CODER.delegating(Shape::polygon, S2Polygon::shape);

    // TODO(user): measure this with a benchmark.
    private static final int MAX_LINEAR_SEARCH_LOOPS = 5; // From benchmarks.
    private static final long serialVersionUID = 1L;

    public S2Polygon polygon() {
      return polygon;
    }

    @Override
    public int numEdges() {
      return polygon.isFull() ? 0 : polygon.numVertices;
    }

    @Override
    public boolean hasInterior() {
      return true;
    }

    @Override
    public boolean containsOrigin() {
      boolean containsOrigin = false;
      for (S2Loop loop : polygon.loops) {
        containsOrigin ^= loop.containsOrigin();
      }
      return containsOrigin;
    }

    @Override
    public int numChains() {
      return polygon.numLoops();
    }

    @Override
    public int getChainLength(int chainId) {
      Preconditions.checkElementIndex(chainId, numChains());
      // Full or empty S2Loops have one vertex, but S2Shape represents full or empty loops as chains
      // with no vertices.
      int numVertices = polygon.loops.get(chainId).numVertices();
      return numVertices == 1 ? 0 : numVertices;
    }

    @Override
    public void getChainEdge(int chainId, int offset, MutableEdge result) {
      Preconditions.checkElementIndex(offset, getChainLength(chainId));
      S2Loop loop = polygon.loop(chainId);
      result.set(loop.orientedVertex(offset), loop.orientedVertex(offset + 1));
    }

    @Override
    public S2Point getChainVertex(int chainId, int edgeOffset) {
      Preconditions.checkElementIndex(edgeOffset, getChainLength(chainId) + 1);
      S2Loop loop = polygon.loop(chainId);
      return loop.orientedVertex(edgeOffset);
    }

    @Override
    public void getChainPosition(int edgeId, ChainPosition result) {
      Preconditions.checkArgument(edgeId < numEdges());

      // When the number of loops is small, linear search is faster than a more complex approach,
      // and most often there is exactly one loop and the code below executes zero times. Anyway,
      // due to the way S2Polygon stores its data, a faster search along the lines of S2ShapeAspect
      // .Multi is not straightforward.
      int chainId = 0;
      for (; edgeId >= polygon.loop(chainId).numVertices(); ++chainId) {
        edgeId -= polygon.loop(chainId).numVertices();
      }
      result.set(chainId, edgeId);
    }

    @Override
    public int dimension() {
      return 2;
    }
  }
}
