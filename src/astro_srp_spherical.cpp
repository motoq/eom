/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_srp_spherical.h>

#include <memory>
#include <utility>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>

namespace eom {

SrpSpherical::SrpSpherical(double cr,
                           double aom,
                           std::unique_ptr<Ephemeris> sun_eph)
{
  m_cr = cr;
  m_aom = aom;
  m_sun_eph = std::move(sun_eph);
}


Eigen::Matrix<double, 3, 1>
SrpSpherical::getAcceleration(const JulianDate& jd,
                              const Eigen::Matrix<double, 6, 1>& state)
{
  Eigen::Matrix<double, 3, 1> r_sat_o = state.block<3,1>(0,0);
  Eigen::Matrix<double, 3, 1> r_sun_o = m_sun_eph->getPosition(jd,
                                                               EphemFrame::eci);
  Eigen::Matrix<double, 3, 1> rhat_sat_sun = (r_sat_o - r_sun_o).normalized();

  return 4.57e-6*m_cr*m_aom*rhat_sat_sun;
}


}
