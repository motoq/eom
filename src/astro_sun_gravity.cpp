/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_sun_gravity.h>

#include <memory>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>

namespace eom {

SunGravity::SunGravity(std::unique_ptr<Ephemeris> sunEph)
{
  m_sunEph = std::move(sunEph);
}


Eigen::Matrix<double, 3, 1>
        SunGravity::getAcceleration(const JulianDate& jd,
                                    const Eigen::Matrix<double, 6, 1>& state)
{
  Eigen::Matrix<double, 3, 1> r_sat_o = state.block<3,1>(0,0);
  Eigen::Matrix<double, 3, 1> r_sun_o = m_sunEph->getPosition(jd,
                                                              EphemFrame::eci);
  Eigen::Matrix<double, 3, 1> r_sat_sun = r_sat_o - r_sun_o;

  double r_sat_sun3 {r_sat_sun.norm()};
  r_sat_sun3 *= r_sat_sun3*r_sat_sun3;
  double r_sun_o3 {r_sun_o.norm()};
  r_sun_o3 *= r_sun_o3*r_sun_o3;

  return -1.0*phy_const::gm_sun*(r_sat_sun/r_sat_sun3 + r_sun_o/r_sun_o3);
}


}
