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
#include <utility>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_duration.h>
#include <cal_julian_date.h>
#include <cal_time_indexer.h>
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
                         const JulianDate& jdStart,
                         const JulianDate& jdStop,
                         const std::shared_ptr<const EcfEciSys>& ecfeciSys,
                         std::unique_ptr<OdeSolver<JulianDate, double, 6>> sp)
{
  m_name = name;
  m_jdStart = jdStart;
  m_jdStop = jdStop;
  m_ecfeciSys = ecfeciSys;

    // SP propagator
  std::unique_ptr<OdeSolver<JulianDate, double, 6>> c_sp = std::move(sp);
  m_jdEpoch = c_sp->getT();
  nullState = c_sp->getX();

    // Pad stop time
  JulianDate jdEndProp {m_jdStop + utl_const::day_per_min};
  JulianDate jdNow = c_sp->getT();

    // Forward ephemeris
  std::vector<eph_record> fwd_eph;
  Eigen::Matrix<double, 6, 1> c_x = c_sp->getX();
  Eigen::Matrix<double, 6, 1> c_dx = c_sp->getXdot();
  fwd_eph.emplace_back(jdNow, c_x.block<3, 1>(0, 0),
                              c_x.block<3, 1>(3, 0),
                              c_dx.block<3, 1>(3, 0));
  while (jdNow < jdEndProp) {
    jdNow = c_sp->step();
    c_x = c_sp->getX();
    c_dx = c_sp->getXdot();
    fwd_eph.emplace_back(jdNow, c_x.block<3, 1>(0, 0),
                                c_x.block<3, 1>(3, 0),
                                c_dx.block<3, 1>(3, 0));
  }

  std::vector<std::pair<JulianDate, JulianDate>> times;
    // Generate and store Hermite interpolation objects
  for (unsigned int ii=1U; ii<fwd_eph.size(); ++ii) {
    eph_record& r1 = fwd_eph[ii-1U];
    eph_record& r2 = fwd_eph[ii];
    double dt_tu {phy_const::tu_per_day*(r2.t - r1.t)};
    Hermite2<double, 3> hItp(dt_tu,
                             r1.p, r1.v, r1.a,
                             r2.p, r2.v, r2.a,
                             phy_const::epsdt);
    m_eph_interpolators.emplace_back(r1.t, r2.t, hItp);
    times.emplace_back(r1.t, r2.t);
  }
  m_ndxr = std::make_unique<TimeIndexer>(times);
}

Eigen::Matrix<double, 6, 1> SpEphemeris::getStateVector(const JulianDate& jd,
                                                        EphemFrame frame) const
{
    // Simple retrieval for now
  Eigen::Matrix<double, 6, 1> xeci = nullState;
/*
  bool found {false};
  for (const auto& irec : m_eph_interpolators) {
    if (irec.jd1 <= jd  &&  jd <= irec.jd2) {
      double dt_tu {phy_const::tu_per_day*(jd - irec.jd1)};
      xeci.block<3,1>(0,0) = irec.hItp.getPosition(dt_tu);
      xeci.block<3,1>(3,0) = irec.hItp.getVelocity(dt_tu);
      found = true;
      break;
    }
  }
  if (!found) {
    std::cout << "\n\nError - Didn't Find it\n\n";
  }
*/

  unsigned long ndx {m_ndxr->getIndex(jd)};
  if (ndx >= 0) {
    const auto& irec = m_eph_interpolators[ndx];
    double dt_tu {phy_const::tu_per_day*(jd - irec.jd1)};
    xeci.block<3,1>(0,0) = irec.hItp.getPosition(dt_tu);
    xeci.block<3,1>(3,0) = irec.hItp.getVelocity(dt_tu);
  } else {
    std::cout << "\n\nError - Didn't Find it\n\n";
  }

  if (frame == EphemFrame::ecf) {
    return m_ecfeciSys->eci2ecf(jd, xeci.block<3,1>(0,0), xeci.block<3,1>(3,0));
  }

  return xeci;
}

}
