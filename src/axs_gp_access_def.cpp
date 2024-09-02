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

GpAccessDef::GpAccessDef(const std::string& orbit_name,
                         const std::string& gp_name,
                         const GpConstraints& xcs,
                         const aux_gp_constraints& axcs,
                         AccessModel  mdl)
{
  m_orbit_name = orbit_name;
  m_gp_name = gp_name;
  m_xcs = xcs;
  m_axcs = axcs;
  m_model = mdl;
}


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
  return m_axcs.use_max_sun_el;
}


aux_gp_constraints GpAccessDef::getAuxConstraints() const noexcept
{
  return m_axcs;
}


}
