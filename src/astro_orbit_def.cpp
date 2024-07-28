/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_orbit_def.h>

#include <string>
#include <array>

#include <utl_const.h>
#include <cal_julian_date.h>
#include <astro_propagator_config.h>
#include <astro_tle.h>

namespace eom {

OrbitDef::OrbitDef(const std::string& orbit_name,
                   const PropagatorConfig& propConfig,
                   const JulianDate& epoch,
                   const std::array<double, 6>& state,
                   CoordType coord_type,
                   FrameType frame_type)
{
  m_name = orbit_name;
  m_propCfg = propConfig,
  m_jd0 = epoch;
  m_x0 = state;
  m_coord = coord_type;
  m_frame = frame_type;
}


OrbitDef::OrbitDef(const std::string& orbit_name,
                   const std::string& tle1,
                   const std::string& tle2)
{
  m_name = orbit_name;

  PropagatorConfig cfg(PropagatorType::sgp4);
  m_propCfg = cfg;

  m_coord = CoordType::tle;
  m_frame = FrameType::teme;

  m_tle.set(tle1, tle2);
  m_jd0 = m_tle.getEpoch();
    // TLE values but in units of radians and minutes
    // n e i o w m
  m_x0[0] = utl_const::tpi*utl_const::day_per_min*m_tle.getMeanMotion();
  m_x0[1] = m_tle.getEccentricity();
  m_x0[2] = utl_const::rad_per_deg*m_tle.getInclination();
  m_x0[3] = utl_const::rad_per_deg*m_tle.getRaan();
  m_x0[4] = utl_const::rad_per_deg*m_tle.getArgumentOfPerigee();
  m_x0[5] = utl_const::rad_per_deg*m_tle.getMeanAnomaly();
}


Tle OrbitDef::getTle() const noexcept
{
  return m_tle;
}


}
