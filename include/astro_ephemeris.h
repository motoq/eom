/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_EPHEMERIS_H
#define ASTRO_EPHEMERIS_H

#include <Eigen/Dense>

#include <cal_julian_date.h>

namespace eom {

/**
 * The GCRF and ITRF are the two reference frames supported throughout
 * the codebase.
 */
enum class EphemFrame {
  eci,                            ///< GCRF (IAU 2000A/2006)
  ecf                             ///< ITRF (~WGS 84)
};

/**
 * Interface defining an ephemeris provider.  Ephemeris is typically
 * assumed to be from an orbit, but can be anything for which ECF or ECI
 * Cartesian postion and velocity make sense.
 */
class Ephemeris {
public:
  virtual ~Ephemeris() {}

  /**
   * @param  jd     UTC time for which to return a state vector
   * @param  frame  Reference frame of returned state vector
   *
   * @return  Cartesian position and velocity state vector, DU and DU/TU
   */
  virtual Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate& jd,
                                                     EphemFrame frame) const=0;
};


}

#endif
