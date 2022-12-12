/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_sp_ephemeris.h>

#include <cal_julian_date.h>
#include <astro_propagator_config.h>

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
  

}

}
