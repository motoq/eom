/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_rk4.h>

#include <utility>
#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <cal_duration.h>
#include <mth_ode.h>
#include <mth_ode_solver.h>
#include <mth_rk4.h>

namespace eom {


Rk4::Rk4(std::unique_ptr<Ode<JulianDate, double, 6>> deq,
         const Duration& dt,
         const JulianDate& jd,
         const Eigen::Matrix<double, 6, 1>& x) : m_deq {std::move(deq)},
                                                 m_dt {dt},
                                                 m_jd0 {jd},
                                                 m_x0 {x}
{
  reset();
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


JulianDate Rk4::step()
{
  mth_rk4::rk4_step(m_deq.get(), m_dt, m_jd,
                    m_x, m_dx, OdeEvalMethod::corrector);

  return m_jd;
}


void Rk4::reset()
{
    // Reset state vector
  m_jd = m_jd0;
  m_x = m_x0;

    // Default integration step size if not explicitly set
  if (m_dt.getTu() == 0.0) {
    Duration dt_default(0.3, phy_const::tu_per_min);
    m_dt = dt_default;
  }
    // Update integration direction
  if (m_dt.getTu() < 0.0) {
    m_prop_dir = StepDirection::reverse;
  } else {
    m_prop_dir = StepDirection::forward;
  }

    // Acceleration at jd
  m_dx = m_deq->getXdot(m_jd, m_x);
}


void Rk4::resetAndReverse()
{
  m_dt = m_dt*(-1.0);
  reset();
}


StepDirection Rk4::getStepDirection() const noexcept
{
  return m_prop_dir;
}


}
