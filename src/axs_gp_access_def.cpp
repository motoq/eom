/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_gp_access_def.h>

#include <string>

namespace eom {

std::string GpAccessDef::getOrbitName() const noexcept
{
  return m_orbit_name;
}


std::string GpAccessDef::getGpName() const noexcept
{
  return m_gp_name;
}


AccessModel GpAccessDef::getAccessModel() const noexcept
{
  return m_model;
}


GpConstraints GpAccessDef::getConstraints() const noexcept
{
  return m_xcs;
}


// Bang together all constraint activation flags
bool GpAccessDef::useAuxConstraints() const noexcept
{
  return m_aux_constraints.use_max_sun_el;
}


aux_gp_constraints GpAccessDef::getAuxConstraints() const noexcept
{
  return m_aux_constraints;
}


}
