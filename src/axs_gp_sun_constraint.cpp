/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_gp_sun_constraint.h>

#include <cmath>
#include <memory>

#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ground_point.h>
#include <astro_sun_meeus.h>

namespace eom {

GpSunConstraint::GpSunConstraint(const GroundPoint& gp,
                                 std::shared_ptr<const EcfEciSys> ecfeciSys) :
                                 m_gp {gp}
{
  m_sunEph = std::make_shared<eom::SunMeeus>(ecfeciSys);
}


bool GpSunConstraint::isSatisfied(JulianDate utc) const
{
    // Check for failures - otherwise, success
  if (m_use_max_el  &&
      m_gp.getSinElevation(m_sunEph->getPosition(utc,
                                     EphemFrame::ecf)) > m_max_sin_el) {
    return false;
  }

  return true;
}


void GpSunConstraint::setMaxElevation(double max_el)
{
  m_use_max_el = true;
  m_max_sin_el = std::sin(max_el);
}


}
