/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_RK4_H
#define ASTRO_RK4_H

#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <cal_duration.h>
#include <mth_ode.h>

#include <mth_ode_solver.h>

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
void rk4_step(Ode<JulianDate, double, 6>* deq,
              const Duration& dt,
              JulianDate& jd,
              Eigen::Matrix<double, 6, 1>& x,
              Eigen::Matrix<double, 6, 1>& dx,
              OdeEvalMethod dx_method = OdeEvalMethod::predictor);

/**
 * Propagates astrodynamics equations of motion using an RK4 integrator.
 *
 * @author  Kurt Motekew
 * @date    2022/09/11
 */
class Rk4 : public OdeSolver<JulianDate, double, 6> {
public:
  ~Rk4() = default;
  Rk4(const Rk4&) = delete;
  Rk4& operator=(const Rk4&) = delete;
  Rk4(Rk4&&) = default;
  Rk4& operator=(Rk4&&) = default;

  /**
   * Initialize with equations of motion, fixed step size,
   * and initial state of the system to be integrated.
   *
   * @param  deq  Equations of motion
   * @param  dt   Integration step size
   * @param  jd   State vector epoch
   * @param  x    Initial conditions - state vector at epoch
   */
  Rk4(std::unique_ptr<Ode<JulianDate, double, 6>> deq,
      const Duration& dt,
      const JulianDate& jd,
      const Eigen::Matrix<double, 6, 1>& x);

  /**
   * Time associated with current state vector and derivative
   */
  JulianDate getT() const noexcept override;

  /**
   * Current state vector
   */
  Eigen::Matrix<double, 6, 1> getX() const noexcept override;

  /**
   * Time derivative of current state vector
   */
  Eigen::Matrix<double, 6, 1> getXdot() const noexcept override;

 /**
   * Propagate forward by system integration step size.
   *
   * @return   Time associated with propagated state.
   */
  JulianDate step() override;

private:
  std::unique_ptr<Ode<JulianDate, double, 6>> m_deq {nullptr};
  Duration m_dt;
  JulianDate m_jd;
  Eigen::Matrix<double, 6, 1> m_x;
  Eigen::Matrix<double, 6, 1> m_dx;
};


}

#endif
