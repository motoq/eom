/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_rel_orbit_def.h>

#include <array>
#include <string>

namespace eom {

RelOrbitDef::RelOrbitDef(const std::string& orbit_name, 
                         const std::string& template_name,
                         const std::array<double, 6>& rel_state, 
                         RelCoordType coord_type)
{
  m_name = orbit_name;
  m_template_name = template_name;
  m_dx0 = rel_state;
  m_coord = coord_type;
}


std::string RelOrbitDef::getOrbitName() const noexcept
{
  return m_name;
}


std::string RelOrbitDef::getTemplateOrbitName() const noexcept
{
  return m_template_name;
}


std::array<double, 6> RelOrbitDef::getInitialState() const noexcept
{
  return m_dx0;
}


RelCoordType RelOrbitDef::getRelCoordinateType() const noexcept 
{
  return m_coord;
}


}
