/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_ECFECI_SYS_H
#define ASTRO_ECFECI_SYS_H

#include <vector>

#include <Eigen/Geometry>

#include <mth_quaternion_interp.h>
#include <cal_julian_date.h>
#include <cal_duration.h>

namespace eom {

/* Enable after adding parsing of IERS EOP
struct raw_eop {
  long mjd {0L};
  double xp {0.0};
  double yp {0.0};
  double ut1mutc {0.0};
  double lod {0.0};
  double dx {0.0};
  double dy {0.0};
};
*/

  // The transformation directoin is ECF to ECI
struct ecf_eci {
  double mjd2000 {0.0};
  double ut1mutc {0.0};
  double lod {0.0};
  Eigen::Quaterniond pm;     ///< Polar motion; ITRF to TIRF
  Eigen::Quaterniond bpn;    ///< Frame bias, precession, nutation; TIRF to CIRF
};

/**
 * System resource utility for ECF/ECI conversions.
 *
 * @author  Kurt Motekew
 * @date    20211020
 */
class EcfEciSys {
public:
  /**
   * This constructor creates an ECF to ECI conversion utility that
   * makes use of a list of generated precession and nutation data.  It
   * allows for the most accurate conversions with some upfront
   * computational burden.  Full benefit of this option comes with
   * enabling interpolation.  Interpolation also comes at a price since
   * it requires interpolating ECF to ECI data vs. directly referencing
   * an existing structure of data in a list.
   *
   * @param  startTime    Earliest UTC time for which ECF to ECI data is
   *                      to be generated or available.
   * @param  startTime    Latest UTC time for which ECF to ECI data is
   *                      to be generated or available.
   * @param  dt           Rate at which to generate ECF to ECI data
   * @param  interpolate  If true, ECF to ECI data will be interpolated.
   *                      Note that UT1-UTC will always be interpolated.
   *                      Defaults to true.
   */
  EcfEciSys(const JulianDate& startTime, const JulianDate& stopTime,
            const Duration& dt, bool interpolate = true);

  ecf_eci getEcfEciData(JulianDate& utc) const;

private:
  JulianDate jdStart;
  JulianDate jdStop;
  double rate_days {0.0};
  unsigned long nfi {0UL};
  bool interpolate_bpnpm {true};
  std::vector<ecf_eci> f2iData;
};


}

#endif
