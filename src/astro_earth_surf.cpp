/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_earth_surf.h>

#include <cmath>

#include <Eigen/Dense>

#include <phy_const.h>
#include <utl_no_solution_exception.h>


namespace eom {

EarthSurf::EarthSurf(HorizonMode mode, double alt)
{
  m_mode = mode;

  a2 = phy_const::earth_smaj + alt;
  a2 *= a2;
  b2 = phy_const::earth_smin + alt;
  b2 *= b2;
}


// Based on "Pointing Vector to Oblate Spheroid",
// Kurt A. Motekew, July 24, 2014
void EarthSurf::setEarthSurf(const Eigen::Matrix<double, 3, 1>& pos,
                             const Eigen::Matrix<double, 3, 1>& pnt)
{
  Eigen::Matrix<double, 3, 1> uhat = pnt.normalized();
  
  double alpha {(uhat[0]*uhat[0] + uhat[1]*uhat[1])/a2 + uhat[2]*uhat[2]/b2};
  double beta  {(uhat[0]*pos[0] + uhat[1]*pos[1])/a2 + uhat[2]*pos[2]/b2};
  double gamma {(pos[0]*pos[0] + pos[1]*pos[1])/a2 + pos[2]*pos[2]/b2};

  double d {beta*beta - alpha*(gamma - 1.0)};

  if (m_mode == HorizonMode::horizon_always  ||
      m_mode == HorizonMode::horizon_miss  &&  d < 0.0) {
    double a1 {std::sqrt((gamma - 1.0)/(gamma*alpha - beta*beta))};
    double a0 {(1.0 - beta*a1)/gamma};
    m_xyz = a0*pos + a1*uhat;
    earth_x = false;
  } else if (m_mode == HorizonMode::horizon_never  &&  d < 0.0) {
    m_xyz.setZero();
    earth_x = false;
  } else {
    double s {-1.0*(beta + std::sqrt(d))/alpha};
    m_xyz = pos + s*uhat;
    earth_x = true;
  }
}


}

