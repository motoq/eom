/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef PHY_CONST_H
#define PHY_CONST_H

#include <cmath>

#include <cal_const.h>

/**
 * Constants and functions related physical properties.  Units defined
 * here are dependent on such properties, such as an ER (earth radius),
 * DU (distance unit), and the TU (time unit).
 *
 * @author  Kurt Motekew
 */
namespace phy_const {

  //
  // Base constants defining the earth and rotation rate
  //

  /** Ellipsoid radius, GRS80/WGS 84, km */
  constexpr double km_per_er {6378.1370};
  /** Ellipsoid flattening, WGS 84 */
  constexpr double flat {1.0/298.257223563};
  constexpr double j2 {1.08262668355e-3};
  constexpr double j3 {-2.5327e-6};
  /** Gravitational parameter, EGM96/EGM2008, TN 36 TT compatible, km^3/s^2 */
  constexpr double gm_km3_sec2 {398600.4415};
  /** Gravitational scaling radius, EGM96/EGM2008, TN 36 TT compatible, km */
  constexpr double km_per_du {6378.1363};
  /** Nominal mean angular velocity of earth w.r.t ECI, GRS80/WGS 84, rad/sec */
  constexpr double we_rad_sec {7292115.0e-11};
  
  //
  // Derived
  //
  
    // Square of ellipsoid eccentricity and eccentricity
  constexpr double ecc2 {flat*(2.0 - flat)};
  constexpr double ecc {std::sqrt(ecc2)};

    // distance units
  constexpr double er_per_km {1.0/km_per_er};
  constexpr double m_per_er {1000.0*km_per_er};
  constexpr double er_per_m {1.0/m_per_er};

  constexpr double du_per_km {1.0/km_per_du};
  constexpr double m_per_du {1000.0*km_per_du};
  constexpr double du_per_m {1.0/m_per_du};

    // Earth ellipsoid in DU
  constexpr double er_per_du {km_per_du/km_per_er};
  constexpr double du_per_er {1.0/er_per_du};
  constexpr double earth_smaj {er_per_du};
  constexpr double earth_smin {earth_smaj*(1.0 - flat)};

  /** Time unit definition */
  constexpr double sec_per_tu {std::sqrt(km_per_du*
                                         ((km_per_du*km_per_du)/gm_km3_sec2))};
  constexpr double min_per_tu {cal_const::min_per_sec*sec_per_tu};
  constexpr double day_per_tu {cal_const::day_per_sec*sec_per_tu};
    //
  constexpr double tu_per_sec {1.0/sec_per_tu};
  constexpr double tu_per_min {1.0/min_per_tu};
  constexpr double tu_per_day {1.0/day_per_tu};

  /** Gravitation parameter and scaling */
  constexpr double gm {1.0};
  constexpr double re {1.0};

  /** 1 mm movement at a sea level orbital altitude */
  constexpr double epsdt {1.0e-6/km_per_du};

  /**
   * Computes the angular velocity of the earth w.r.t. inertial space.
   * This is the correct angular velocity to use when accounting for the
   * Coriolis effect and/or centripetal acceleration.  Supplying this
   * function with a literal "0.0" when LOD is not available provides a
   * more accurate angular velocity than we_rad_sec above while incuring
   * no additional computational burden since it will be evaluated at
   * compile time.
   *
   * @param  lod  Length of day, TU
   *
   * @return  Earth angular velocity, radians/TU
   */
  constexpr double earth_angular_velocity(double lod)
  {
    return 7.292115146706979e-5*sec_per_tu*(1.0 - lod/tu_per_day);
  }
}

#endif
