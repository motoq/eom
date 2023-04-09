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


Rk4::Rk4(std::unique_ptr<Ode<JulianDate, double, 6>> deq,
         const Duration& dt,
         const JulianDate& jd,
         const Eigen::Matrix<double, 6, 1>& x)
{
  m_deq = std::move(deq);
  m_dt = dt;
  m_jd = jd;
  m_x = x;
    // Default integration step size if not explicitly set
  if (m_dt.getTu() == 0.0) {
    Duration dt_default(0.3, phy_const::tu_per_min);
    m_dt = dt_default;
  }
    // Acceleration at jd
  m_dx = m_deq->getXdot(m_jd, m_x);
}


JulianDate Rk4::getT() const noexcept
{
  return m_jd;
}


Eigen::Matrix<double, 6, 1> Rk4::getX() const noexcept
{
  return m_x;
}


Eigen::Matrix<double, 6, 1> Rk4::getXdot() const noexcept
{
  return m_dx;
}


/*
 * RK4 algorithm adapted from "Aircraft Control and Simulation" by
 * Brian L. Stevens and Frank L. Lewis, 1st ed.
 */
JulianDate Rk4::step()
{
  auto dt = m_dt.getTu();
  auto dt_days = m_dt.getDays();

  Eigen::Matrix<double, 6, 1> x0 = m_x;
  Eigen::Matrix<double, 6, 1> xd;
  Eigen::Matrix<double, 6, 1> xx;
  Eigen::Matrix<double, 6, 1> xa;
  Eigen::Matrix<double, 6, 1> q;

    // first
  auto jdNow = m_jd;
  xd = m_deq->getXdot(jdNow, x0);
  xa = dt*xd;
  xx = 0.5*xa + x0;
    // second
  jdNow += 0.5*dt_days;
  xd = m_deq->getXdot(jdNow, xx);
  q = dt*xd;
  xx = x0 + 0.5*q;
  xa += q + q;
    // third
  xd = m_deq->getXdot(jdNow, xx);
  q = dt*xd;
  xx = x0 + q;
  xa += q + q;
    // forth - update member variables vs. locals
  m_jd += dt_days;
  m_dx = m_deq->getXdot(m_jd, xx);
  m_x = x0 + (xa + dt*m_dx)/6.0;

  return m_jd;
}


std::unique_ptr<Ode<JulianDate, double, 6>> Rk4::returnDeq()
{
  return std::move(m_deq);
}


}
