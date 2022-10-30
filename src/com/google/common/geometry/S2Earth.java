package com.google.common.geometry;

import java.io.Serializable;

public class S2Earth implements Serializable {

    /**
     * The Earth's mean radius in meters, which is the radius of the equivalent sphere with the same surface area.
     * According to NASA, this value is 6371.01 +/- 0.02 km.  The equatorial radius is 6378.136 km, and the polar
     * radius is 6356.752 km.  They differ by one part in 298.257.
     * <p>
     * Reference: http://ssd.jpl.nasa.gov/phys_props_earth.html, which quotes
     * Yoder, C.F. 1995. "Astrometric and Geodetic Properties of Earth and the
     * Solar System" in Global Earth Physics, A Handbook of Physical Constants,
     * AGU Reference Shelf 1, American Geophysical Union, Table 2.
     */
    public static final double RADIUS_METERS = 6371010.0;

    /**
     * The Earth's mean radius in kilometers.
     */
    public static final double RADIUS_KM = 0.001 * RADIUS_METERS;
    public S1Angle metersToAngle(double meters){
       return S1Angle.radians(metersToRadians(meters));
    }

    public static double metersToRadians(double meters) {
        return meters / RADIUS_METERS;
    }
}
