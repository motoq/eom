/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_sp_ephemeris.h>

#include <string>
#include <memory>

#include <cal_duration.h>
#include <cal_julian_date.h>
#include <astro_propagator_config.h>

#include <phy_const.h>
#include <astro_gravity.h>
#include <astro_gravity_jn.h>
#include <astro_deq.h>
#include <mth_ode_solver.h>
#include <astro_rk4.h>

namespace eom {

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

  std::unique_ptr<OdeSolver<JulianDate, double, 6>> sp {nullptr};
  Duration dt(1.0, phy_const::tu_per_min);
  sp = std::make_unique<Rk4>(std::move(deq), dt);
}

}
