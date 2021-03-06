/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_vinti.h>

#include <array>
#include <memory>

#include <Eigen/Dense>

#include <Vinti.h>

#include <cal_const.h>
#include <cal_julian_date.h>
#include <phy_const.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>

namespace eom {


Vinti::Vinti(const JulianDate& epoch,
             const Eigen::Matrix<double, 6, 1>& xeci,
             const std::shared_ptr<const EcfEciSys>& ecfeciSys,
             VintiPertModel pertModel)
{
  if (pertModel == VintiPertModel::J2_ONLY) {
    planet[3] = 0.0;
  }
  jd0 = epoch;
  ecfeci = ecfeciSys;
    // A true equator ECI frame is required for propagation
  Eigen::Matrix<double, 6, 1> xecf = ecfeci->eci2ecf(jd0,
                                                     xeci.block<3,1>(0,0),
                                                     xeci.block<3,1>(3,0));
  Eigen::Matrix<double, 6, 1> xteme = ecfeci->ecf2teme(jd0,
                                                       xecf.block<3,1>(0,0),
                                                       xecf.block<3,1>(3,0));
  x0[0] = phy_const::km_per_du*xteme(0,0);
  x0[1] = phy_const::km_per_du*xteme(1,0);
  x0[2] = phy_const::km_per_du*xteme(2,0);
  x0[3] = phy_const::km_per_du*xteme(3,0)*phy_const::tu_per_sec;
  x0[4] = phy_const::km_per_du*xteme(4,0)*phy_const::tu_per_sec;
  x0[5] = phy_const::km_per_du*xteme(5,0)*phy_const::tu_per_sec;
}


Eigen::Matrix<double, 6, 1> Vinti::getStateVector(const JulianDate& jd,
                                                  EphemFrame frame) const
{
  std::array<double, 6> oe;
  std::array<double, 6> x1;
  double t1 {cal_const::sec_per_day*(jd - jd0)};
  Vinti6(planet.data(), 0.0, x0.data(), t1, x1.data(), oe.data());

  Eigen::Matrix<double, 6, 1> xteme;
  xteme(0,0) = phy_const::du_per_km*x1[0];
  xteme(1,0) = phy_const::du_per_km*x1[1];
  xteme(2,0) = phy_const::du_per_km*x1[2];
  xteme(3,0) = phy_const::du_per_km*x1[3]*phy_const::sec_per_tu;
  xteme(4,0) = phy_const::du_per_km*x1[4]*phy_const::sec_per_tu;
  xteme(5,0) = phy_const::du_per_km*x1[5]*phy_const::sec_per_tu;
  Eigen::Matrix<double, 6, 1> xecf = ecfeci->teme2ecf(jd,
                                                      xteme.block<3,1>(0,0),
                                                      xteme.block<3,1>(3,0));

  if (frame == EphemFrame::eci) {
    Eigen::Matrix<double, 6, 1> xeci = ecfeci->ecf2eci(jd,
                                                       xecf.block<3,1>(0,0),
                                                       xecf.block<3,1>(3,0));
    return xeci;
  }
  return xecf;
}


}
