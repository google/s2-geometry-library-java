package com.google.common.geometry;

/**
 * A kind of query target.
 */
interface S2MinDistanceTarget {
    /**
     * Returns the approximate center of the target.
     */
    S2Point center();

    /**
     * Returns the distance between this target and the given cell.
     */
    S1ChordAngle getDistance(S2Cell cell);

    /**
     * Returns the radian radius of an angular cap that encloses this target.
     */
    double radius();

    /**
     * Returns the smaller of {@code distance} and a new distance from target to {@code point}.
     */
    S1ChordAngle getMinDistance(S2Point point, S1ChordAngle minDist);

    S1ChordAngle getMinDistance(S2Point a, S2Point b, S1ChordAngle minDist);

    default S1ChordAngle getMinDistance(S2Cell cell, S1ChordAngle minDist) {
        S1ChordAngle angle = getDistance(cell);
        return angle.compareTo(minDist) > 0 ? minDist : angle;
    }

    default boolean visitContainingShapes(S2ShapeIndex index, S2ContainsPointQuery.ShapeVisitor visitor) {
        return new S2ContainsPointQuery(index, S2ContainsPointQuery.Options.SEMI_OPEN).visitContainingShapes(center(), visitor);
    }

    default int maxBruteForceIndexSize() {
        return -1;
    }

    default S2Cap getCapBound() {
        return S2Cap.fromAxisAngle(center(), S1Angle.radians(radius()));
    }
}


/** A point query, used to find the closest points to a query point. */
abstract class S2MinDistancePointTarget implements S2MinDistanceTarget {
    private final S2Point point;

    public S2MinDistancePointTarget(S2Point point) {
        this.point = point;
    }

    @Override
    public S2Point center() {
        return point;
    }

    @Override
    public double radius() {
        return 0;
    }

    @Override
    public S1ChordAngle getMinDistance(S2Point x, S1ChordAngle minDist) {
        S1ChordAngle angle = new S1ChordAngle(x, point);
        // See comment regarding ">=" in the findClosestPoints() main loop.
        return angle.compareTo(minDist) > 0 ? minDist : angle;
    }

    @Override
    public S1ChordAngle getMinDistance(S2Point a, S2Point b, S1ChordAngle minDist) {
        return S2EdgeUtil.updateMinDistance(point, a, b, minDist);
    }

    @Override
    public S1ChordAngle getDistance(S2Cell cell) {
        return cell.getDistance(point);
    }
}


/** An edge query, used to find the closest points to a query edge. */
abstract class S2MinDistanceEdgeTarget implements S2MinDistanceTarget {
    private S2Point a;
    private S2Point b;

    public S2MinDistanceEdgeTarget(S2Point a, S2Point b) {
        this.a = a;
        this.b = b;
    }

    @Override
    public S2Point center() {
        return S2Point.normalize(S2Point.add(a, b));
    }

    @Override
    public double radius() {
        return 0.5 * a.angle(b);
    }

    @Override
    public S1ChordAngle getMinDistance(S2Point x, S1ChordAngle minDist) {
        return S2EdgeUtil.updateMinDistance(x, a, b, minDist);
    }

    @Override
    public S1ChordAngle getMinDistance(S2Point c, S2Point d, S1ChordAngle minDist) {
        return S2EdgeUtil.getEdgePairMinDistance(a, b, c, d, minDist);
    }

    @Override
    public S1ChordAngle getDistance(S2Cell cell) {
        return cell.getDistanceToEdge(a, b);
    }
}
