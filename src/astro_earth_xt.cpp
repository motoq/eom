/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_earth_xt.h>

#include <cmath>

#include <Eigen/Dense>

#include <phy_const.h>
#include <mth_unit_circle.h>
#include <utl_no_solution_exception.h>
#include <astro_earth_surf.h>

namespace {
  constexpr double ainv {1.0/phy_const::earth_smaj};
  constexpr double binv {1.0/phy_const::earth_smin};
}

namespace eom {

EarthXt::EarthXt(HorizonMode mode, double alt)
{
  m_mode = mode;

    // Initialize transformation matrices
  if (alt == 0) {
    t_ac(0,0) = ainv;
    t_ac(1,1) = ainv;
    t_ac(2,2) = binv;
    t_ca(0,0) = phy_const::earth_smaj;
    t_ca(1,1) = phy_const::earth_smaj;
    t_ca(2,2) = phy_const::earth_smin;
  } else {
    double a_alt {alt + phy_const::earth_smaj};
    double b_alt {alt + phy_const::earth_smin};
    t_ac(0,0) = 1.0/a_alt;
    t_ac(1,1) = 1.0/a_alt;
    t_ac(2,2) = 1.0/b_alt;
    t_ca(0,0) = a_alt;
    t_ca(1,1) = a_alt;
    t_ca(2,2) = b_alt;
  }
}


// Based on A Fine Way to Kiss the Hyperellipsoid
// Kurt A. Motekew, May 30, 2022
void EarthXt::setEarthXt(const Eigen::Matrix<double, 3, 1>& pos,
                         const Eigen::Matrix<double, 3, 1>& pnt)
{
  Eigen::Matrix<double, 3, 1> posAff {t_ac*pos};
  Eigen::Matrix<double, 3, 1> pntAff {t_ac*pnt};
    // 3D to 2D
  Eigen::Matrix<double, 3, 1> yhat {posAff};
  Eigen::Matrix<double, 3, 1> zhat {posAff.cross(pntAff)};
  Eigen::Matrix<double, 3, 1> xhat {yhat.cross(zhat)};
  xhat.normalize();
  yhat.normalize();
  zhat.normalize();
  Eigen::Matrix<double, 3, 3> dcm;
  dcm.row(0) = xhat;
  dcm.row(1) = yhat;
  dcm.row(2) = zhat;

  Eigen::Matrix<double, 3, 1> pos2D = dcm*posAff;
  Eigen::Matrix<double, 3, 1> pnt2D = dcm*pntAff;

  Eigen::Matrix<double, 2, 1> r = pos2D.block<2,1>(0,0);
  Eigen::Matrix<double, 2, 1> p = pnt2D.block<2,1>(0,0);
  Eigen::Matrix<double, 2, 1> xy;
  if (m_mode == HorizonMode::horizon_always) {
    xy = unit_circle::tangent(r, p);
    earth_x = false;
  } else {
    try {
      xy = unit_circle::intersect(r, p);
      earth_x = true;
    } catch (const NoSolutionException& nsl) {
      earth_x = false;
      if (m_mode == HorizonMode::horizon_miss) {
        xy = unit_circle::tangent(r, p);
      } else if (m_mode == HorizonMode::horizon_never) {
        m_xyz.setZero();
        return;
      }
    }
  }
  Eigen::Matrix<double, 3, 1> xyz2D = {xy(0), xy(1), 0.0};
  dcm.transposeInPlace();
  m_xyz = t_ca*dcm*xyz2D;
}


}

