/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_ECFECI_SYS_H
#define ASTRO_ECFECI_SYS_H

#include <memory>

#include <Eigen/Geometry>

#include <cal_julian_date.h>
#include <cal_duration.h>

namespace eom {

struct raw_eop {
  long mjd {0L};
  double xp {0.0};
  double yp {0.0};
  double ut1mutc {0.0};
  double lod {0.0};
  double dx {0.0};
  double dy {0.0};
};

  // The transformation directoin is ECF to ECI
struct ecf_eci {
  double mjd {0.0};
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
   * Default zero EcfEciSys
   */
  EcfEciSys(JulianDate startTime, JulianDate stopTime, Duration dt);

private:
  JulianDate jdStart;
  JulianDate jdStop;
  double dt_days {0.0};
  long nfi {0L};
  std::unique_ptr<ecf_eci> f2iData {nullptr};
};


}

#endif
