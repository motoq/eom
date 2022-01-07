/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_rel_orbit_def.h>

#include <string>
#include <array>

namespace eom {

RelOrbitDef::RelOrbitDef(std::string orbit_name, std::string template_name,
               const std::array<double, 6>& rel_state, RelCoordType coord_type)
{
  name = orbit_name;
  ref_name = template_name;
  dx0 = rel_state;
  coord = coord_type;
}


}
