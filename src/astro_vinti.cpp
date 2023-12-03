/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_vinti.h>

#include <string>
#include <array>
#include <utility>
#include <memory>

#include <Eigen/Dense>

#include <Vinti.h>

#include <cal_const.h>
#include <cal_julian_date.h>
#include <phy_const.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>

namespace eom {


Vinti::Vinti(const std::string& orbit_name,
             const JulianDate& epoch,
             const Eigen::Matrix<double, 6, 1>& xeci,
             std::shared_ptr<const EcfEciSys> ecfeciSys,
             VintiPertModel pertModel)
{
  if (pertModel == VintiPertModel::J2_ONLY) {
    m_planet[3] = 0.0;
  }
  m_name = orbit_name;
  m_jd0 = epoch;
  m_ecfeci = std::move(ecfeciSys);
    // A true equator ECI frame is required for propagation
  Eigen::Matrix<double, 6, 1> xecf = m_ecfeci->eci2ecf(m_jd0,
                                                       xeci.block<3,1>(0,0),
                                                       xeci.block<3,1>(3,0));
  Eigen::Matrix<double, 6, 1> xteme = m_ecfeci->ecf2teme(m_jd0,
                                                         xecf.block<3,1>(0,0),
                                                         xecf.block<3,1>(3,0));
  m_x0[0] = phy_const::km_per_du*xteme(0);
  m_x0[1] = phy_const::km_per_du*xteme(1);
  m_x0[2] = phy_const::km_per_du*xteme(2);
  m_x0[3] = phy_const::km_per_du*xteme(3)*phy_const::tu_per_sec;
  m_x0[4] = phy_const::km_per_du*xteme(4)*phy_const::tu_per_sec;
  m_x0[5] = phy_const::km_per_du*xteme(5)*phy_const::tu_per_sec;
}


Eigen::Matrix<double, 6, 1> Vinti::getStateVector(const JulianDate& jd,
                                                  EphemFrame frame) const
{
  std::array<double, 6> oe;
  std::array<double, 6> x1;
  double t1 {cal_const::sec_per_day*(jd - m_jd0)};
  Vinti6(m_planet.data(), 0.0, m_x0.data(), t1, x1.data(), oe.data());

  Eigen::Matrix<double, 6, 1> xteme;
  xteme(0) = phy_const::du_per_km*x1[0];
  xteme(1) = phy_const::du_per_km*x1[1];
  xteme(2) = phy_const::du_per_km*x1[2];
  xteme(3) = phy_const::du_per_km*x1[3]*phy_const::sec_per_tu;
  xteme(4) = phy_const::du_per_km*x1[4]*phy_const::sec_per_tu;
  xteme(5) = phy_const::du_per_km*x1[5]*phy_const::sec_per_tu;
  Eigen::Matrix<double, 6, 1> xecf = m_ecfeci->teme2ecf(jd,
                                                        xteme.block<3,1>(0,0),
                                                        xteme.block<3,1>(3,0));

  if (frame == EphemFrame::eci) {
    return m_ecfeci->ecf2eci(jd, xecf.block<3,1>(0,0), xecf.block<3,1>(3,0));
  }
  return xecf;
}


Eigen::Matrix<double, 3, 1> Vinti::getPosition(const JulianDate& jd,
                                               EphemFrame frame) const
{
  std::array<double, 6> oe;
  std::array<double, 6> x1;
  double t1 {cal_const::sec_per_day*(jd - m_jd0)};
  Vinti6(m_planet.data(), 0.0, m_x0.data(), t1, x1.data(), oe.data());

  Eigen::Matrix<double, 3, 1> xteme;
  xteme(0) = phy_const::du_per_km*x1[0];
  xteme(1) = phy_const::du_per_km*x1[1];
  xteme(2) = phy_const::du_per_km*x1[2];
  Eigen::Matrix<double, 3, 1> xecf = m_ecfeci->teme2ecf(jd, xteme);

  if (frame == EphemFrame::eci) {
    return m_ecfeci->ecf2eci(jd, xecf);
  }
  return xecf;
}


}
