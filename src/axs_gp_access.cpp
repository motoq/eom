/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_gp_access.h>

#include <string>
#include <memory>

#include <cal_julian_date.h>
#include <astro_ground_point.h>
#include <astro_ephemeris.h>
#include <axs_gp_constraints.h>

namespace eom {

static bool is_visible(const JulianDate& jd);

GpAccess::GpAccess(const JulianDate& jdStart,
                   const JulianDate& jdStop,
                   const GroundPoint& gp,
                   const GpConstraints& xcs,
                   std::shared_ptr<const Ephemeris> eph) : m_jdStart {jdStart},
                                                           m_jdStop {jdStop},
                                                           m_gp {gp},
                                                           m_xcs {xcs},
                                                           m_eph {
                                                             std::move(eph)
                                                           }
{
  is_visible(m_jdStart);
}


std::string GpAccess::getGpName() const noexcept
{
  return m_gp.getName();
}


std::string GpAccess::getOrbitName() const
{
  return (*m_eph).getName();
}


static bool is_visible(const JulianDate& jd)
{
  return true;
}


}
