/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_gp_access.h>

#include <memory>

#include <cal_julian_date.h>
#include <astro_ground_point.h>
#include <astro_ephemeris.h>
#include <axs_gp_constraints.h>

namespace eom {

static bool is_visible(const JulianDate& jd, 
                       const GroundPoint& gp,
                       const Ephemeris& eph,
                       const GpConstraints& xcs);

GpAccess::GpAccess(const JulianDate& jdStart,
                   const JulianDate& jdStop,
                   const GroundPoint& gp,
                   const Ephemeris& eph,
                   const GpConstraints& xcs)
{
  m_gp_name = gp.getName();
  m_eph_name = eph.getName();
  
  is_visible(jdStart, gp, eph, xcs);
}


static bool is_visible(const JulianDate& jd, 
                       const GroundPoint& gp,
                       const Ephemeris& eph,
                       const GpConstraints& xcs)
{
  return true;
}


}
