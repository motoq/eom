/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_deq.h>

#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_gravity.h>

namespace eom {


Deq::Deq(std::unique_ptr<Gravity> grav,
         const std::shared_ptr<const EcfEciSys>& ecfeci)
{
  m_ecfeci = ecfeci;
  m_grav = std::move(grav);
}


void Deq::getXDot(const JulianDate& utc,
                  const Eigen::Matrix<double, 6, 1>& x,
                        Eigen::Matrix<double, 6, 1>& xd)
{
  Eigen::Matrix<double, 3, 1> posi = x.block<3,1>(0,0);
  Eigen::Matrix<double, 3, 1> veli = x.block<3,1>(3,0);

    // Derivative of position is just the velocity
  xd.block<3,1>(0,0) = veli;
  
    // Force model
  Eigen::Matrix<double, 3, 1> posf = m_ecfeci->eci2ecf(utc, posi);
  Eigen::Matrix<double, 3, 1> a_i_f =
      m_grav->getAcceleration(posf, eom::EvalMethod::predictor);
    // Derivative in inertial, convert components from fixed to inertial
  Eigen::Matrix<double, 3, 1> a_i_i = m_ecfeci->ecf2eci(utc, a_i_f);
  xd.block<3,1>(3,0) = a_i_i;
}


}
