/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_deq.h>

#include <utility>
#include <memory>
#include <vector>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <mth_ode.h>
#include <astro_gravity.h>

namespace eom {


Deq::Deq(std::unique_ptr<Gravity> grav, std::shared_ptr<const EcfEciSys> ecfeci)
{
  m_ecfeci = std::move(ecfeci);
  m_grav = std::move(grav);
}


Eigen::Matrix<double, 6, 1> Deq::getXdot(const JulianDate& utc,
                                         const Eigen::Matrix<double, 6, 1>& x,
                                         OdeEvalMethod method)
{
    // Central body gravity model
  Eigen::Matrix<double, 3, 1> posf = m_ecfeci->eci2ecf(utc, x.block<3,1>(0,0));
  Eigen::Matrix<double, 3, 1> a_i_f = m_grav->getAcceleration(posf, method);

  Eigen::Matrix<double, 6, 1> xd;
    // Velocity is derivative of position
  xd.block<3,1>(0,0) = x.block<3,1>(3,0);
    // Acceleration derivative is w.r.t. ECI, but need to transform
    // vector components to ECI frame
  xd.block<3,1>(3,0) = m_ecfeci->ecf2eci(utc, a_i_f);

    // Add non-central body accelerations
  if (m_fmodels.size() > 0) {
    Eigen::Matrix<double, 3, 1> a_i_i = Eigen::Matrix<double, 3, 1>::Zero();
    for (auto& fm : m_fmodels) {
      a_i_i += fm->getAcceleration(utc, x);
    }
    xd.block<3,1>(3,0) += a_i_i;
  }

  return xd;
}


void Deq::addForceModel(std::unique_ptr<ForceModel> fm)
{
  m_fmodels.push_back(std::move(fm));
}


}
