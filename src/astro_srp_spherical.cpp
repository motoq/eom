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

#include <mth_angle.h>
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

  Eigen::Matrix<double, 3, 1> rhat_sat_o = r_sat_o.normalized();
  Eigen::Matrix<double, 3, 1> rhat_sun_o = r_sun_o.normalized();
  if (r_sat_o.dot(r_sun_o) > 0.0) {
    auto angle = mth_angle::unit_vec_angle<double>(rhat_sat_o, rhat_sun_o);
    if (angle*r_sat_o.norm() < 1.0) {
      return Eigen::Matrix<double, 3, 1>::Zero();
    }
  }

  Eigen::Matrix<double, 3, 1> rhat_sat_sun = (r_sat_o - r_sun_o).normalized();

  double cr_i {0.0};
  double cr_r {m_cr - 1.0};
  if (cr_r > 0.0) {
    cr_i = m_cr - cr_r;
  } else {
    cr_r = m_cr;
    cr_i = 0.0;
  }

  return 4.57e-6*m_aom*(cr_i*rhat_sat_sun - 
                        cr_r*rhat_sat_o);
}


}
