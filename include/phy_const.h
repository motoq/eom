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
  constexpr double KM_PER_ER {6378.1370};
  /** Ellipsoid flattening, WGS 84 */
  constexpr double FLAT {1.0/298.257223563};
  /** Gravitational parameter, EGM96/EGM2008, TN 36 TT compatible, km^3/s^2 */
  constexpr double GM_KM3_SEC2 {398600.4415};
  /** Gravitational scaling radius, EGM96/EGM2008, TN 36 TT compatible, km */
  constexpr double KM_PER_DU {6378.1363};
  /** Nominal mean angular velocity of earth w.r.t ECI, GRS80/WGS 84, rad/sec */
  constexpr double WE_RAD_SEC {7292115.0e-11};
  
  
  //
  // Derived
  //

  
  /** Time unit definition */
  constexpr double SEC_PER_TU {std::sqrt(KM_PER_DU*
                                         ((KM_PER_DU*KM_PER_DU)/GM_KM3_SEC2))};
  constexpr double MIN_PER_TU {cal_const::MIN_PER_SEC*SEC_PER_TU};
  constexpr double DAY_PER_TU {cal_const::DAY_PER_SEC*SEC_PER_TU};
    //
  constexpr double TU_PER_SEC {1.0/SEC_PER_TU};
  constexpr double TU_PER_MIN {1.0/MIN_PER_TU};
  constexpr double TU_PER_DAY {1.0/DAY_PER_TU};
  /** Gravitation parameter */
  constexpr double GM {1.0};
}

#endif
