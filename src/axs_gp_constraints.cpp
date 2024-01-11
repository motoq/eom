/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_gp_constraints.h>

#include <string>
#include <cmath>
#include <stdexcept>

#include <Eigen/Dense>

#include <utl_const.h>
#include <cal_julian_date.h>
#include <astro_ground_point.h>

namespace eom {

void GpConstraints::setMinEl(double min_el)
{
  using namespace std::string_literals;
  if (min_el < -utl_const::pio2  ||  min_el > utl_const::pio2) {
    throw std::invalid_argument("(GpConstraints::setMinEl() "s +
                                "invalid Minimum Elevation "s +
                                std::to_string(min_el) + " radians"s);
  }
  m_sin_min_el = std::sin(min_el);
}


void GpConstraints::setMaxEl(double max_el)
{
  using namespace std::string_literals;
  if (max_el < -utl_const::pio2  ||  max_el > utl_const::pio2) {
    throw std::invalid_argument("(GpConstraints::setMaxEl() "s +
                                "invalid Maximum Elevation "s +
                                std::to_string(max_el) + " radians"s);
  }
  m_sin_max_el = std::sin(max_el);
  m_only_min_el = false;
}


void GpConstraints::setMinMaxAz(double min_az, double max_az)
{
  using namespace std::string_literals;
  if (min_az < 0.0  ||  min_az > utl_const::tpi) {
    throw std::invalid_argument("(GpConstraints::setMinMaxAz() "s +
                                "invalid Minimum Azimuth "s +
                                std::to_string(min_az) + " radians"s);
  }
  if (max_az < 0.0  ||  max_az > utl_const::tpi) {
    throw std::invalid_argument("(GpConstraints::setMinMaxAz() "s +
                                "invalid Maximum Azimuth "s +
                                std::to_string(max_az) + " radians"s);
  }
  m_min_az = min_az;
  m_max_az = max_az;
  m_only_min_el = false;
  m_check_az = true;

  if (m_min_az < m_max_az) {
    m_az_shift = 0.0;
    m_min_az_shifted = m_min_az;
    m_max_az_shifted = m_max_az;
  } else {
    m_az_shift =  m_max_az;
    m_max_az_shifted = utl_const::tpi;
    m_min_az_shifted = m_min_az - m_az_shift;
  }
}


bool GpConstraints::isVisible(const JulianDate& jd, 
                              const GroundPoint& gp,
                              const Eigen::Matrix<double, 3, 1>& pos) const
{
  double sel {gp.getSinElevation(pos)};

  return sel >= m_sin_min_el  &&  sel <= m_sin_max_el;
}


}
