/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_gravity_jn.h>

#include <cmath>
#include <string>
#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>
#include <mth_ode.h>

namespace eom
{

GravityJn::GravityJn(int degree) : nterm {degree}
{
  if (degree > GravityJn::getMaxDegree()) {
    throw std::invalid_argument("GravityJn::GravityJn() Unsupported Degree: " +
                                std::to_string(degree));
  }
}


/*
 * Based on Vallado's "Fundamentals of Astrodynamics and Applications",
 * 3rd edition, Section 8.7.1 "Simplified Acceleration Model", modified
 * for use with GM = 1 DU^3/TU^2 and Re (gravitational scaling
 * parameter) = 1 DU
 */
Eigen::Matrix<double, 3, 1>
    GravityJn::getAcceleration(const Eigen::Matrix<double, 3, 1>& pos,
                               OdeEvalMethod)
{
    // Always used
  double ri {pos(0)};
  double rj {pos(1)};
  double rk {pos(2)};
  double rk2 {rk*rk};
    //
  double rmag2 {ri*ri + rj*rj + rk*rk};
  double invr2 {1.0/rmag2};
  double invr {std::sqrt(invr2)};
    // Always used unless 2-body
  double invr5 {0.0};
  if (nterm > 1) {
    invr5 = invr*invr2*invr2;
  }
  
    // 2-body contribution
  Eigen::Matrix<double, 3, 1> acc {-invr2*(invr*pos)};

    // Intentional fall through - include contributions
    // smallest to largest
  switch (nterm) {
    case 4:
    {
      double invr7 {invr2*invr5};
      double c1 {15.0*phy_const::j4*invr7/8.0};
      double c2 {7.0*rk2*invr2};
      double c3 {3.0*rk2*invr2};
      acc(0) += c1*ri*(1.0 - c2*(2.0 - c3));
      acc(1) += c1*rj*(1.0 - c2*(2.0 - c3));
      acc(2) += c1*rk*(5.0 - c2*(10.0/3.0 - c3));
    }
    case 3:
    {
      double invr7 {invr2*invr5};
      double c1 {2.5*phy_const::j3*invr7};
      double c2 {7.0*rk2*invr2};
      acc(0) -= c1*ri*rk*(3.0 - c2);
      acc(1) -= c1*rj*rk*(3.0 - c2);
      acc(2) -= c1*(rk2*(6.0 - c2) - 3.0*rmag2/5.0);
    }
    case 2:
    {
      double c1 {1.5*phy_const::j2*invr5};
      double c2 {5.0*rk2*invr2};
      acc(0) -= c1*ri*(1.0 - c2);
      acc(1) -= c1*rj*(1.0 - c2);
      acc(2) -= c1*rk*(3.0 - c2);
    }
  }

  return acc;
}


}
