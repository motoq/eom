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
  /** EGM96 Zonals */
  constexpr double j2 {1.082626173852223e-03};
  constexpr double j3 {-2.532410518567722e-06};
  constexpr double j4 {-1.619897599916973e-06};
  constexpr double j5 {-2.277535907308362e-07};
  constexpr double j6 {5.406665762838132e-07};
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

  /** Solar gravitational parameter, TN36, TCB-compatible value, DU^3/TU^2 */
  constexpr double gm_sun {1.32712442099e20*du_per_m*du_per_m*du_per_m*
                                            sec_per_tu*sec_per_tu};
  /** Astronomical unit, TN36, DU */
  constexpr double du_per_au {1.49597870700e11*du_per_m};

  /** TN36 moon/earth mass ratio */
  constexpr double moon_per_earth_mass {1.0/81.3005690699};
  /** Lunar gravitational parameter, TT-compatible value, DU^3/TU^2 */
  constexpr double gm_moon {moon_per_earth_mass*gm};

  /** Planetary mass, TN36 */
  constexpr double sun_per_earth_mass {332946.048166};
    //
  constexpr double mercury_per_earth_mass {sun_per_earth_mass/6023597.400017};
  constexpr double gm_mercury {mercury_per_earth_mass*gm};
    //
  constexpr double venus_per_earth_mass {sun_per_earth_mass/408523.718655};
  constexpr double gm_venus {venus_per_earth_mass*gm};
    //
  constexpr double mars_per_earth_mass {sun_per_earth_mass/3098703.590267};
  constexpr double gm_mars {mars_per_earth_mass*gm};
    //
  constexpr double jupiter_per_earth_mass {sun_per_earth_mass/1047.348625};
  constexpr double gm_jupiter {jupiter_per_earth_mass*gm};
    //
  constexpr double saturn_per_earth_mass {sun_per_earth_mass/3497.901768};
  constexpr double gm_saturn {saturn_per_earth_mass*gm};
    //
  constexpr double uranus_per_earth_mass {sun_per_earth_mass/22902.981613};
  constexpr double gm_uranus {uranus_per_earth_mass*gm};
    //
  constexpr double neptune_per_earth_mass {sun_per_earth_mass/19412.237346};
  constexpr double gm_neptune {neptune_per_earth_mass*gm};
    //
  constexpr double pluto_per_earth_mass {sun_per_earth_mass/135836683.767599};
  constexpr double gm_pluto {pluto_per_earth_mass*gm};

  /** 1 mm movement at a sea level orbital altitude */
  constexpr double epsdt {1.0e-6*du_per_km};
  constexpr double epsdt_days {epsdt*phy_const::day_per_tu};

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
  constexpr double earth_angular_velocity(double lod) noexcept
  {
    return 7.292115146706979e-5*sec_per_tu*(1.0 - lod/tu_per_day);
  }

  /**
   * Compute the magnitude of the velocity of the surface of the earth
   * at the equator.
   *
   * @param  lod  Length of day, TU
   *
   * @return  Speed of earth surface at equator, DU/TU
   */
  constexpr double earth_equatorial_speed(double lod = 0.0) noexcept
  {
    return earth_smaj*earth_angular_velocity(lod);
  }
  
}

#endif
