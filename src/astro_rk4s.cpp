/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_rk4s.h>

#include <utility>
#include <memory>

#include <Eigen/Dense>

#include <utl_units.h>
#include <cal_julian_date.h>
#include <mth_ode.h>
#include <astro_regularize.h>

namespace eom {

Rk4s::Rk4s(std::unique_ptr<Ode<JulianDate, double, 6>> deq,
           const JulianDate& jd,
           const Eigen::Matrix<double, 6, 1>& x) : m_deq {std::move(deq)},
                                                   m_jd00 {jd},
                                                   m_x00 {x}
{
  reset();
}


/**
 * RK4 algorithm adapted from "Aircraft Control and Simulation" by
 * Brian L. Stevens and Frank L. Lewis, 1st ed.
 *
 * Regularization still a bit experimental...
 */
JulianDate Rk4s::step()
{
  using namespace utl_units;

  const double ds {m_ds/16.0};
  const JulianDate jd0 {m_jd0};

    // first
  Eigen::Matrix<double, 8, 1> y0 = m_reg->getY();
  Eigen::Matrix<double, 8, 1> yd = m_reg->getYdot();
  Eigen::Matrix<double, 8, 1> ya = ds*yd;
  Eigen::Matrix<double, 8, 1> yy = 0.5*ya + y0;
    // second
  yd.block<4,1>(0,0) = yy.block<4,1>(4,0);
  m_reg->setRegularizedState(yy, yd);                      // yd ??
  JulianDate jd = jd0 + m_reg->getTime().getDays();
  Eigen::Matrix<double, 6, 1> x = m_reg->getX();
  Eigen::Matrix<double, 6, 1> dx = m_deq->getXdot(jd, x);
  m_reg->setTimeState({jd - jd0, 1.0_min}, x, dx);
  yd = m_reg->getYdot();
  Eigen::Matrix<double, 8, 1> q = ds*yd;
  yy = y0 + 0.5*q;
  ya += q + q;
    // third
  yd.block<4,1>(0,0) = yy.block<4,1>(4,0);
  m_reg->setRegularizedState(yy, yd);                      // yd ??
  jd = jd0 + m_reg->getTime().getDays();
  x = m_reg->getX();
  dx = m_deq->getXdot(jd, x);
  m_reg->setTimeState({jd - jd0, 1.0_min}, x, dx);
  yd = m_reg->getYdot();
  q = ds*yd;
  yy = y0 + q;
  ya += q + q;
    // forth - update member variables vs. locals
  yd.block<4,1>(0,0) = yy.block<4,1>(4,0);
  m_reg->setRegularizedState(yy, yd);                      // yd ??
  jd = jd0 + m_reg->getTime().getDays();
  x = m_reg->getX();
  dx = m_deq->getXdot(jd, x);
  m_reg->setTimeState({jd - jd0, 1.0_min}, x, dx);
  yd = m_reg->getYdot();
  yy = y0 + (ya + ds*yd)/6.0;
  yd.block<4,1>(0,0) = yy.block<4,1>(4,0);
  m_reg->setRegularizedState(yy, yd);                      // yd ??

  return jd0 + m_reg->getTime().getDays();
}


void Rk4s::reset()
{
  m_jd0 = m_jd00;
  Eigen::Matrix<double, 6, 1> dx = m_deq->getXdot(m_jd0, m_x00,
                                                  OdeEvalMethod::predictor);
    // Initialize regularization
  m_reg = std::make_unique<Regularize>(m_x00, dx);
  if (m_prop_dir == StepDirection::forward) {
    m_ds = std::fabs(m_reg->getDsMax());
  } else {
    m_ds = -1.0*std::fabs(m_reg->getDsMax());
  }
  m_reg->setTimeState({0.0, phy_const::tu_per_min}, m_x00, dx);
}


void Rk4s::resetAndReverse()
{
  if (m_prop_dir == StepDirection::forward) {
    m_prop_dir = StepDirection::reverse;
  } else {
    m_prop_dir = StepDirection::forward;
  }
  reset();
}


StepDirection Rk4s::getStepDirection() const noexcept
{
  return m_prop_dir;
}


}
