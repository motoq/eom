/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_RK4_H
#define MTH_RK4_H

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <cal_duration.h>
#include <mth_ode.h>

namespace eom {

  /**
   * Propagate equations of motion forward in time by dt.
   *
   * @param  deq        Equations of motion
   * @param  dt         Integration step size
   * @param  jd         Input state vector epoch
   *                    Output time of output state vector, jd + dt
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
template<typename T, int DIM>
void rk4_step(Ode<JulianDate, T, DIM>* deq,
              const Duration& dt,
              JulianDate& jd,
              Eigen::Matrix<T, DIM, 1>& x,
              Eigen::Matrix<T, DIM, 1>& dx,
              OdeEvalMethod dx_method = OdeEvalMethod::predictor)
{
  auto dt_tu = dt.getTu();
    // No forward integration - just populate derivative
  if (dt_tu == 0.0) {
    dx = deq->getXdot(jd, x);
  }
  auto dt_days = dt.getDays();

  Eigen::Matrix<T, DIM, 1> x0 = x;
  Eigen::Matrix<T, DIM, 1> xd;
  Eigen::Matrix<T, DIM, 1> xx;
  Eigen::Matrix<T, DIM, 1> xa;
  Eigen::Matrix<T, DIM, 1> q;

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

#endif
