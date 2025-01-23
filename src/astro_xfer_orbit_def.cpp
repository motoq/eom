/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_xfer_orbit_def.h>

#include <string>

#include <cal_duration.h>
#include <cal_julian_date.h>

namespace eom {

XferOrbitDef::XferOrbitDef(const std::string& orbit_name,
                           const std::string& orbit1_name,
                           const std::string& orbit2_name,
                           const JulianDate& xferStart,
                           const Duration& dur)
{
  m_name = orbit_name;
  m_start_name = orbit1_name;
  m_end_name = orbit2_name;
  m_start = xferStart;
  m_dur = dur;
}


std::string XferOrbitDef::getOrbitName() const noexcept
{
  return m_name;
}


std::string XferOrbitDef::getStartOrbitName() const noexcept
{
  return m_start_name;
}


std::string XferOrbitDef::getEndOrbitName() const noexcept
{
  return m_end_name;
}


JulianDate XferOrbitDef::getXferStartTime() const noexcept
{
  return m_start;
}


Duration XferOrbitDef::getXferDuration() const noexcept
{
  return m_dur;
}



}
