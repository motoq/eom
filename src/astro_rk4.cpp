/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_rk4.h>

#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <cal_duration.h>
#include <mth_ode.h>

namespace eom {


Rk4::Rk4(std::unique_ptr<Ode<JulianDate, double, 6>> deq, const Duration& dt)
{
  m_deq = std::move(deq);
  m_dt = dt;
}


JulianDate Rk4::step(const JulianDate& utc,
                     const Eigen::Matrix<double, 6, 1>& x0,
                           Eigen::Matrix<double, 3, 1>& a0,
                           Eigen::Matrix<double, 6, 1>& x,
                           Eigen::Matrix<double, 3, 1>& a)
{
  auto dt = m_dt.getDays();

  Eigen::Matrix<double, 6, 1> xd;
  Eigen::Matrix<double, 6, 1> xx;
  Eigen::Matrix<double, 6, 1> xa;

    // first
  auto jdNow = utc;
  m_deq->getXDot(jdNow, x0, xd);
  a0 = xd.block<3, 1>(3, 0);
  for (int ii=0; ii<6; ++ii) {
    xa(ii) = dt*xd(ii);
    xx(ii) = 0.5*xa(ii) + x0(ii);
  }
    // second
  jdNow += 0.5*dt;
  m_deq->getXDot(jdNow, xx, xd);
  for (int ii=0; ii<6; ++ii) {
    auto q = dt*xd(ii);
    xx(ii) =  x0(ii) + 0.5*q;
    xa(ii) += q + q;
  }
    // third
  m_deq->getXDot(jdNow, xx, xd);
  for (int ii=0; ii<6; ++ii) {
    auto q = dt*xd(ii);
    xx(ii) = x0(ii) + q;
    xa(ii) += q + q;
  }
    // forth
  jdNow = utc + dt;
  m_deq->getXDot(jdNow, xx, xd);
  a = xd.block<3, 1>(3, 0);
  for (int ii=0; ii<6; ++ii) {
    x(ii) = x0(ii) + (xa(ii) + dt*xd(ii))/6.0;
  }

  return jdNow;
}


}
