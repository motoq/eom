/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_gravity_std.h>

#include <string>
#include <memory>
#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>
#include <mth_ode.h>
#include <astro_egm_coeff.h>
#include <astro_gravity.h>
#include <mth_legendre_af.h>

namespace eom {

/*
 * Coefficients in astro_egm_coeff.h are expected to be unormalized and
 * continuous, sorted by order first, then degree.
 */
GravityStd::GravityStd(int max_degree, int max_order)
{
  if (max_order > max_degree) {
    throw std::invalid_argument("GravityStd::GravityStd Order > Degree");
  }
  if (max_degree < 0  ||  max_degree > egm_coeff::degree) {
    throw std::invalid_argument(
        "GravityStd::GravityStd() Unsupported Degree: " +
        std::to_string(max_degree)
    );
  }
  if (max_order < 0  ||  max_order > egm_coeff::order) {
    throw std::invalid_argument(
        "GravityStd::GravityStd() Unsupported Order: " +
        std::to_string(max_order)
    );
  }

  m_degree = max_degree;
  m_order = max_order;
  m_alf = std::make_unique<LegendreAf>(m_degree, m_order);
}


Eigen::Matrix<double, 3, 1>
    GravityStd::getAcceleration(const Eigen::Matrix<double, 3, 1>& pos,
                           OdeEvalMethod entry)
{
  double rmag {pos.norm()};
  double invr {1.0/rmag};
  //double rxy {std::sqrt(pos(0)*pos(0) + pos(1)*pos(1))};
  //double invrxy {1.0/rxy};
    
  return -1.0*phy_const::gm*invr*invr*invr*pos;
}


}
