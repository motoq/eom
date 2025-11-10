/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_RK4_H
#define MTH_RK4_H

#include <Eigen/Dense>

#include <mth_ode.h>

namespace mth_rk4 {

/**
 * Integrates a 1st order ordinary differential equation using
 * RK4 algorithm adapted from "Aircraft Control and Simulation" by
 * Brian L. Stevens and Frank L. Lewis, 1st ed.
 *
 * @tparam  T    Time type
 * @tparam  DT   Duration compatible with T type and time units (TU)
 * @tparam  F    Data type of state vector
 * @tparam  DIM  State vector dimension
 *
 * @param  deq        Equations of motion
 * @param  dt         Integration step size
 * @param  time       Input state vector epoch
 *                    Output time of output state vector, time + dt
 * @param  x          Input initial conditions - state vector at epoch
 *                    Output propagated state vector
 * @param  dx         Derivative of state vector at output time
 * @param  dx_method  When computing dx to be returned, predictor
 *                    (default) will reevaluate the full spherical
 *                    harmonics using the final update to x.  The
 *                    corrector option will reuse past accumulated
 *                    partials and update only the 2-body portion.
 *                    This option saves time when dx is needed as
 *                    an output for Hermite interpolation.
 */
template<typename T, typename DT, typename F, int DIM>
void rk4_step(eom::Ode<T, F, DIM>* deq,
              const DT& dt,
              T& time,
              Eigen::Matrix<F, DIM, 1>& x,
              Eigen::Matrix<F, DIM, 1>& dx,
              eom::OdeEvalMethod dx_method = eom::OdeEvalMethod::predictor)
{
  constexpr F half {static_cast<F>(0.5)};
  constexpr F sixth {static_cast<F>(1.0/6.0)};

  auto dt_tu = dt.getTu();
    // No forward integration - just populate derivative
  if (dt_tu == 0) {
    dx = deq->getXdot(time, x);
    return;
  }

  Eigen::Matrix<F, DIM, 1> x0 = x;
  Eigen::Matrix<F, DIM, 1> xd;
  Eigen::Matrix<F, DIM, 1> xx;
  Eigen::Matrix<F, DIM, 1> xa;
  Eigen::Matrix<F, DIM, 1> q;

    // first
  auto tmpTime = time;
  xd = deq->getXdot(tmpTime, x0);
  xa = dt_tu*xd;
  xx = half*xa + x0;
    // second
  tmpTime += dt*half;
  xd = deq->getXdot(tmpTime, xx);
  q = dt_tu*xd;
  xx = x0 + half*q;
  xa += q + q;
    // third
  xd = deq->getXdot(tmpTime, xx);
  q = dt_tu*xd;
  xx = x0 + q;
  xa += q + q;
    // forth - update member variables vs. locals
  time += dt;
  dx = deq->getXdot(time, xx);
  x = x0 + (xa + dt_tu*dx)*sixth;
  dx = deq->getXdot(time, x, dx_method);
}

}

#endif
