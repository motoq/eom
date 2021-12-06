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

#include <cal_julian_date.h>
#include <astro_propagator_config.h>

namespace eom {

OrbitDef::OrbitDef(std::string orbit_name, const PropagatorConfig& propConfig,
                   const JulianDate& epoch, const std::array<double, 6>& state,
                                    CoordType coord_type, FrameType frame_type)
{
  name = orbit_name;
  propCfg = propConfig,
  jd0 = epoch;
  x0 = state;
  coord = coord_type;
  frame = frame_type;
}


}
