/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_kepler.h>

#include <array>
#include <memory>
#include <string>
#include <utility>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>

#include <Vinti.h>

namespace eom {


Kepler::Kepler(const std::string& orbit_name,
               const JulianDate& epoch,
               const Eigen::Matrix<double, 6, 1>& xeci,
               std::shared_ptr<const EcfEciSys> ecfeciSys)
{
  m_name = orbit_name;
  m_jd0 = epoch;
  m_ecfeci = std::move(ecfeciSys);
  m_x0[0] = phy_const::km_per_du*xeci(0);
  m_x0[1] = phy_const::km_per_du*xeci(1);
  m_x0[2] = phy_const::km_per_du*xeci(2);
  m_x0[3] = phy_const::km_per_du*xeci(3)*phy_const::tu_per_sec;
  m_x0[4] = phy_const::km_per_du*xeci(4)*phy_const::tu_per_sec;
  m_x0[5] = phy_const::km_per_du*xeci(5)*phy_const::tu_per_sec;
}


Eigen::Matrix<double, 6, 1> Kepler::getStateVector(const JulianDate& jd,
                                                   EphemFrame frame) const 
{
  double x {0.0};
  std::array<double, 6> x1;
  double t1 {utl_const::sec_per_day*(jd - m_jd0)};
  Kepler1(m_planet.data(), 0.0, m_x0.data(), t1, x1.data(), &x);

  Eigen::Matrix<double, 6, 1> xeci;
  xeci(0) = phy_const::du_per_km*x1[0];
  xeci(1) = phy_const::du_per_km*x1[1];
  xeci(2) = phy_const::du_per_km*x1[2];
  xeci(3) = phy_const::du_per_km*x1[3]*phy_const::sec_per_tu;
  xeci(4) = phy_const::du_per_km*x1[4]*phy_const::sec_per_tu;
  xeci(5) = phy_const::du_per_km*x1[5]*phy_const::sec_per_tu;

  if (frame == EphemFrame::ecf) {
    return m_ecfeci->eci2ecf(jd, xeci.block<3,1>(0,0), xeci.block<3,1>(3,0));
  }
  return xeci;
}


Eigen::Matrix<double, 3, 1> Kepler::getPosition(const JulianDate& jd,
                                                EphemFrame frame) const
{
  double x {0.0};
  std::array<double, 6> x1;
  double t1 {utl_const::sec_per_day*(jd - m_jd0)};
  Kepler1(m_planet.data(), 0.0, m_x0.data(), t1, x1.data(), &x);

  Eigen::Matrix<double, 3, 1> xeci;
  xeci(0) = phy_const::du_per_km*x1[0];
  xeci(1) = phy_const::du_per_km*x1[1];
  xeci(2) = phy_const::du_per_km*x1[2];

  if (frame == EphemFrame::ecf) {
    return m_ecfeci->eci2ecf(jd, xeci);
  }
  return xeci;
}


}
