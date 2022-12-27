/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_sp_ephemeris.h>

#include <iostream>

#include <string>
#include <vector>
#include <memory>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_duration.h>
#include <cal_julian_date.h>
#include <mth_hermite2.h>
#include <mth_ode_solver.h>
#include <astro_ephemeris.h>

namespace eom {

/*
 * Work in progress...
 * Need to add interpolator finder, exceptions for out of range,
 * boundary/edge resolution
 */
SpEphemeris::SpEphemeris(const std::string& name,
                         const JulianDate& epoch,
                         const Eigen::Matrix<double, 6, 1>& xeci,
                         const JulianDate& jdStart,
                         const JulianDate& jdStop,
                         const std::shared_ptr<const EcfEciSys>& ecfeciSys,
                         std::unique_ptr<OdeSolver<JulianDate, double, 3>> sp)
{
  m_name = name;
  m_jdEpoch = epoch;
  m_jdStart = jdStart;
  m_jdStop = jdStop;
  nullState = xeci;
  m_ecfeciSys = ecfeciSys;

    // SP propagator
  std::unique_ptr<OdeSolver<JulianDate, double, 3>> c_sp = std::move(sp);

    // Pad stop time
  JulianDate jdEndProp {m_jdStop + utl_const::day_per_min};
  JulianDate jdNow = epoch;

    // Forward ephemeris
  std::vector<eph_record> fwd_eph;
  Eigen::Matrix<double, 6, 1> x0 = xeci;
  Eigen::Matrix<double, 3, 1> a0;
  Eigen::Matrix<double, 6, 1> x1;
  Eigen::Matrix<double, 3, 1> a1;
    // Need acceleration vector for first time
  c_sp->step(jdNow, x0, a0, x1, a1);
  fwd_eph.emplace_back(jdNow, x0.block<3, 1>(0, 0), x0.block<3, 1>(3, 0), a0);
  while (jdNow < jdEndProp) {
    jdNow = c_sp->step(jdNow, x0, a0, x1, a1);
    fwd_eph.emplace_back(jdNow, x1.block<3, 1>(0, 0), x1.block<3, 1>(3, 0), a1);
    x0 = x1;
  }

    // Generate and store Hermite interpolation objects
  for (unsigned int ii=1U; ii<fwd_eph.size(); ++ii) {
    eph_record& r1 = fwd_eph[ii-1U];
    eph_record& r2 = fwd_eph[ii];
    double dt_tu {phy_const::tu_per_day*(r2.t - r1.t)};
    Hermite2<double, 3> hItp(dt_tu, r1.p, r1.v, r1.a, r2.p, r2.v, r2.a);
    m_eph_interpolators.emplace_back(r1.t, r2.t, hItp);
  }
}

Eigen::Matrix<double, 6, 1> SpEphemeris::getStateVector(const JulianDate& jd,
                                                        EphemFrame frame) const
{
    // Simple retrieval for now
  Eigen::Matrix<double, 6, 1> xeci = nullState;
  bool found {false};
  for (const auto& interp_record : m_eph_interpolators) {
    if (interp_record.jd1 <= jd  &&  jd <= interp_record.jd2) {
      double dt_tu {phy_const::tu_per_day*(jd - interp_record.jd1)};
      xeci.block<3,1>(0,0) = interp_record.hItp.getPosition(dt_tu);
      xeci.block<3,1>(3,0) = interp_record.hItp.getVelocity(dt_tu);
      found = true;
      break;
    }
  }
  if (!found) {
    std::cout << "\n\nError - Didn't Find it\n\n";
  }

  if (frame == EphemFrame::ecf) {
    return m_ecfeciSys->eci2ecf(jd, xeci.block<3,1>(0,0), xeci.block<3,1>(3,0));
  }

  return xeci;
}

}
