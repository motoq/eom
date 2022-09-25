/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_gravity_jn.h>

#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>

namespace {
  constexpr int max_deg {3};
}

namespace eom
{

GravityJn::GravityJn(int degree) : nterm {degree}
{
  if (degree > max_deg) {
    throw std::runtime_error("GravityJn::GravityJn Unsupported Degree");
  }
}

/*
 * GM = 1, Re = 1
 */
Eigen::Matrix<double, 3, 1>
    GravityJn::getAcceleration(const Eigen::Matrix<double, 3, 1>& pos,
                               EvalMethod)
{
  double ri {pos(0)};
  double rj {pos(1)};
  double rk {pos(2)};
  double rk2 {rk*rk};

  double rmag {pos.norm()};
  double invr {1.0/rmag};
  double invr2 {invr*invr};

  double invr5 {0.0};
  if (nterm > 1) {
    invr5 = invr*invr2*invr2;
  }
  
  Eigen::Matrix<double, 3, 1> acc {-invr2*(invr*pos)};

    // Intentional fall through - include contributions
    // smallest to largest
  switch (nterm) {
    case 3:
    {
      double invr7 {invr2*invr5};
      double c1 {2.5*phy_const::j3*invr7};
      double c2 {7.0*rk2*invr2};
      acc(0) -= c1*ri*rk*(3.0 - c2);
      acc(1) -= c1*rj*rk*(3.0 - c2);
      acc(2) -= c1*(rk2*(6.0 - c2) - 3.0*rmag*rmag/5.0);
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
