/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_sp_ephemeris.h>

#include <string>
#include <vector>
#include <memory>

#include <cal_duration.h>
#include <cal_julian_date.h>
#include <astro_propagator_config.h>

#include <utl_const.h>
#include <phy_const.h>
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
  m_jdStart = epoch;
  m_jdStop = propCfg.getStopTime();
  nullState = xeci;
  m_ecfeciSys = ecfeciSys;

  std::unique_ptr<Gravity> forceModel {nullptr};
  forceModel = std::make_unique<GravityJn>(2);
  auto deq = std::make_unique<Deq>(std::move(forceModel), ecfeciSys);

  std::unique_ptr<OdeSolver<JulianDate, double, 3>> sp {nullptr};
  Duration dt(1.0, phy_const::tu_per_min);
  sp = std::make_unique<Rk4>(std::move(deq), dt);

  double dt_days {dt.getDays()};
  //JulianDate jdStart = propCfg.getStartTime();
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
  }

  double dt_sec {utl_const::sec_per_day*dt_days};
  //for (const auto& feph : fwd_eph) {
    
  
}

}
