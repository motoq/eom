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

#include <algorithm>

#include <cal_duration.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_propagator_config.h>

#include <utl_const.h>
#include <phy_const.h>
#include <mth_hermite.h>
#include <astro_gravity.h>
#include <astro_gravity_jn.h>
#include <astro_deq.h>
#include <mth_ode_solver.h>
#include <astro_rk4.h>

namespace eom {

/*
 * Work in progress - testing implementation of the tools needed for SP
 * methods (ODE, integrator, gravity model, interpolation) before
 * setting up factory like process using PropagatorConfig settings.
 */
SpEphemeris::SpEphemeris(const std::string& name,
                         const JulianDate& epoch,
                         const Eigen::Matrix<double, 6, 1>& xeci,
                         const std::shared_ptr<const EcfEciSys>& ecfeciSys,
                         const PropagatorConfig& propCfg)
{
  m_name = name;
  m_jdEpoch = epoch;
  m_jdStart = propCfg.getStartTime();
  m_jdStop = propCfg.getStopTime();
  nullState = xeci;
  m_ecfeciSys = ecfeciSys;

    // SP propagator
  std::unique_ptr<OdeSolver<JulianDate, double, 3>> sp {nullptr};
  Duration dt(1.0, phy_const::tu_per_min);
  {
    std::unique_ptr<Gravity> forceModel {nullptr};
    forceModel = std::make_unique<GravityJn>(2);
    auto deq = std::make_unique<Deq>(std::move(forceModel), ecfeciSys);
    sp = std::make_unique<Rk4>(std::move(deq), dt);
  }


  double dt_days {dt.getDays()};
  JulianDate jdStop = propCfg.getStopTime() + 2.0*dt_days;
  JulianDate jdNow = epoch;

    // Forward ephemeris
  std::vector<eph_record> fwd_eph;
  Eigen::Matrix<double, 6, 1> x0 = xeci;
  Eigen::Matrix<double, 3, 1> a0;
  Eigen::Matrix<double, 6, 1> x1;
  Eigen::Matrix<double, 3, 1> a1;
    // Need acceleration vector for first time
  sp->step(jdNow, x0, a0, x1, a1);
  fwd_eph.emplace_back(jdNow, x0.block<3, 1>(0, 0), x0.block<3, 1>(3, 0), a0);
  while (jdNow < jdStop) {
    jdNow = sp->step(jdNow, x0, a0, x1, a1);
    fwd_eph.emplace_back(jdNow, x1.block<3, 1>(0, 0), x1.block<3, 1>(3, 0), a1);
    x0 = x1;
  }

  m_dt_factor = dt_days;
  for (unsigned int ii=1U; ii<fwd_eph.size(); ++ii) {
    eph_record& r1 = fwd_eph[ii-1U];
    eph_record& r2 = fwd_eph[ii];
    Hermite<double, 3> hItp(1.0, r1.p, r1.v, r1.a, r2.p, r2.v, r2.a);
    m_eph_interpolators.emplace_back(r1.t, r2.t, hItp);
  }
  
}

Eigen::Matrix<double, 6, 1> SpEphemeris::getStateVector(const JulianDate& jd,
                                                        EphemFrame frame) const
{
  Eigen::Matrix<double, 6, 1> xeci = nullState;
  bool found {false};
  for (const auto& interp_record : m_eph_interpolators) {
    if (interp_record.jd1 <= jd  &&  jd <= interp_record.jd2) {
      double dt {m_dt_factor*(jd - interp_record.jd1)};
      xeci.block<3,1>(0,0) = interp_record.hItp.getX(dt);
      xeci.block<3,1>(3,0) = interp_record.hItp.getdX(dt);
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
