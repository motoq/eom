/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_third_body_gravity.h>

#include <memory>
#include <utility>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>

namespace eom {

ThirdBodyGravity::ThirdBodyGravity(double gm,
                                   std::unique_ptr<Ephemeris> eph)
{
  m_gm = gm;
  m_eph = std::move(eph);
}


Eigen::Matrix<double, 3, 1>
    ThirdBodyGravity::getAcceleration(const JulianDate& jd,
                                      const Eigen::Matrix<double, 6, 1>& state)
{
  Eigen::Matrix<double, 3, 1> r_sat_o = state.block<3,1>(0,0);
  Eigen::Matrix<double, 3, 1> r_3rd_o = m_eph->getPosition(jd,
                                                           EphemFrame::eci);
  Eigen::Matrix<double, 3, 1> r_sat_3rd = r_sat_o - r_3rd_o;

  double r_sat_3rd3 {r_sat_3rd.norm()};
  r_sat_3rd3 *= r_sat_3rd3*r_sat_3rd3;
  double r_3rd_o3 {r_3rd_o.norm()};
  r_3rd_o3 *= r_3rd_o3*r_3rd_o3;

  return -1.0*m_gm*(r_sat_3rd/r_sat_3rd3 + r_3rd_o/r_3rd_o3);
}


}
