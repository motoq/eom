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
  rk4_step(m_deq.get(), m_dt, m_jd, m_x, m_dx, OdeEvalMethod::corrector);

  return m_jd;
}


void rk4_step(Ode<JulianDate, double, 6>* deq,
              const Duration& dt,
              JulianDate& jd,
              Eigen::Matrix<double, 6, 1>& x,
              Eigen::Matrix<double, 6, 1>& dx,
              OdeEvalMethod dx_method)
{
  auto dt_tu = dt.getTu();
    // No forward integration - just populate derivative
  if (dt_tu == 0.0) {
    dx = deq->getXdot(jd, x);
  }
  auto dt_days = dt.getDays();

  Eigen::Matrix<double, 6, 1> x0 = x;
  Eigen::Matrix<double, 6, 1> xd;
  Eigen::Matrix<double, 6, 1> xx;
  Eigen::Matrix<double, 6, 1> xa;
  Eigen::Matrix<double, 6, 1> q;

    // first
  auto jdNow = jd;
  xd = deq->getXdot(jdNow, x0);
  xa = dt_tu*xd;
  xx = 0.5*xa + x0;
    // second
  jdNow += 0.5*dt_days;
  xd = deq->getXdot(jdNow, xx);
  q = dt_tu*xd;
  xx = x0 + 0.5*q;
  xa += q + q;
    // third
  xd = deq->getXdot(jdNow, xx);
  q = dt_tu*xd;
  xx = x0 + q;
  xa += q + q;
    // forth - update member variables vs. locals
  jd += dt;
  dx = deq->getXdot(jd, xx);
  x = x0 + (xa + dt_tu*dx)/6.0;
  dx = deq->getXdot(jd, x, dx_method);
}


}
