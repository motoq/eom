/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_adams_4th.h>

#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <cal_duration.h>
#include <mth_ode.h>
#include <astro_rk4.h>

namespace eom {


Adams4th::Adams4th(std::unique_ptr<Ode<JulianDate, double, 6>> deq,
                   const Duration& dt,
                   const JulianDate& jd,
                   const Eigen::Matrix<double, 6, 1>& x)
{
  m_dt = dt;
    // Default integration step size if not explicitly set
  if (m_dt.getTu() == 0.0) {
    Duration dt_default(0.3, phy_const::tu_per_min);
    m_dt = dt_default;
  }

  Rk4 rk(std::move(deq), m_dt, jd, x);
  m_jdW[0] = rk.getT();
  m_w[0] = rk.getX();
  m_dw[0] = rk.getXdot();
  for (int ii=1; ii<4; ++ii) {
    rk.step();
    m_jdW[ii] = rk.getT();
    m_w[ii] = rk.getX();
    m_dw[ii] = rk.getXdot();
  }

  m_jd = m_jdW[0];
  m_x = m_w[0];
  m_dx = m_dw[0];
    // Retrieve eom - done with rk
  m_deq = rk.returnDeq();
}


JulianDate Adams4th::getT() const noexcept
{
  return m_jd;
}


Eigen::Matrix<double, 6, 1> Adams4th::getX() const noexcept
{
  return m_x;
}


Eigen::Matrix<double, 6, 1> Adams4th::getXdot() const noexcept
{
  return m_dx;
}


JulianDate Adams4th::step()
{

  auto dt = m_dt.getTu();
//  auto dt_days = m_dt.getDays();

  Eigen::Matrix<double, 6, 1> wNow = m_w[3] + dt*(55.0*m_dw[3] -
                                                  59.0*m_dw[2] +
                                                  37.0*m_dw[1] -
                                                  9.0*m_dw[0])/24.0;
  JulianDate jdNow = m_jdW[3] + m_dt;
  Eigen::Matrix<double, 6, 1> dwNow = m_deq->getXdot(jdNow, wNow);
  wNow = m_w[3] + dt*(9.0*dwNow + 19.0*m_dw[3] - 5.0*m_dw[2] + m_dw[1])/24.0;

  for (int ii=0; ii<3; ++ii) {
    m_jdW[ii] = m_jdW[ii+1];
    m_w[ii] = m_w[ii+1];
    m_dw[ii] = m_dw[ii+1];
  }
  m_jdW[3] = jdNow;
  m_w[3] = wNow;
  m_dw[3] = m_deq->getXdot(jdNow, wNow);

  m_jd = m_jdW[2];
  m_x = m_w[2];
  m_dx = m_dw[2];

  return m_jd;
}


}
