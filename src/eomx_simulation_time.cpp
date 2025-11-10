/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <utility>

#include <phy_const.h>
#include <astro_orbit_def.h>
#include <cal_julian_date.h>

#include <eom_config.h>
#include <eomx.h>

/**
 * See eomx.h
 */
std::pair<eom::JulianDate, eom::JulianDate>
eomx_simulation_time(const eom_app::EomConfig& cfg,
                     const std::vector<eom::OrbitDef>& orbit_defs)
{
    // Determine time span that must be supported by the simulation
    // based on the input scenario time and orbit epoch times.
  eom::JulianDate minJd = cfg.getStartTime();
  eom::JulianDate maxJd = cfg.getStopTime();
  for (const auto& orbit : orbit_defs) {
    if (orbit.getEpoch() < minJd) {
      minJd = orbit.getEpoch();
    }
    if (maxJd < orbit.getEpoch()) {
      maxJd = orbit.getEpoch();
    }
  }

  return std::make_pair(minJd, maxJd);
}
