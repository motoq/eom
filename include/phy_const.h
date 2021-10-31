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

namespace phy_const {

  //
  // Base constants defining the earth and rotation rate
  //

  /** Ellipsoid radius, GRS80/WGS 84, km */
  constexpr double km_per_er {6378.1370};
  /** Ellipsoid flattening, WGS 84 */
  constexpr double flat {1.0/298.257223563};
  /** Gravitational parameter, EGM96/EGM2008, TN 36 TT compatible, km^3/s^2 */
  constexpr double gm_km3_sec2 {398600.4415};
  /** Gravitational scaling radius, EGM96/EGM2008, TN 36 TT compatible, km */
  constexpr double km_per_du {6378.1363};
  /** Nominal mean angular velocity of earth w.r.t ECI, GRS80/WGS 84, rad/sec */
  constexpr double we_rad_sec {7292115.0e-11};
  
  //
  // Derived
  //
  
  /** Time unit definition */
  constexpr double sec_per_tu {std::sqrt(km_per_du*
                                         ((km_per_du*km_per_du)/gm_km3_sec2))};
  constexpr double min_per_tu {cal_const::min_per_sec*sec_per_tu};
  constexpr double day_per_tu {cal_const::day_per_sec*sec_per_tu};
    //
  constexpr double tu_per_sec {1.0/sec_per_tu};
  constexpr double tu_per_min {1.0/min_per_tu};
  constexpr double tu_per_day {1.0/day_per_tu};

  /** Gravitation parameter */
  constexpr double gm {1.0};

  /** 1 mm movement at a sea level orbital altitude */
  constexpr double epsdt {1.0e-6/km_per_du};
}

#endif
