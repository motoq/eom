/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_RK4S_H
#define ASTRO_RK4S_H

#include <memory>
#include <array>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <mth_ode.h>
#include <astro_regularize.h>

#include <mth_ode_solver.h>

namespace eom {

/**
 * Propagates astrodynamics equations of motion using time
 * regularization with a fixed step-size RK4 integrator.
 *
 * @author  Kurt Motekew
 * @date    2023/05/29
 */
class Rk4s : public OdeSolver<JulianDate, double, 6> {
public:
  /**
   * Initialize with equations of motion, fixed step size,
   * and initial state of the system to be integrated.
   *
   * @param  deq  Equations of motion
   * @param  jd   State vector epoch
   * @param  x    Initial conditions - state vector at epoch
   */
  Rk4s(std::unique_ptr<Ode<JulianDate, double, 6>> deq,
            const JulianDate& jd,
            const Eigen::Matrix<double, 6, 1>& x);

  /**
   * Time associated with current state vector and derivative
   */
  JulianDate getT() const noexcept override
  {
    return m_jd0 + m_reg->getTime().getDays();
  }

  /**
   * Current state vector
   */
  Eigen::Matrix<double, 6, 1> getX() const noexcept override
  {
    return m_reg->getX();
  }

  /**
   * Time derivative of current state vector
   */
  Eigen::Matrix<double, 6, 1> getXdot() const noexcept override
  {
    return m_reg->getXdot();
  }

  /**
   * Propagate forward by system integration step size.
   *
   * @return   Time associated with propagated state.
   */
  JulianDate step() override;

  /**
   * Reset the integrator to the orginal state (the state it was
   * upon instantiation) with the current integration direction
   * (StepDirection remains the same)
   */
  void reset() override;

  /**
   * Reset the integrator and reverse the direction of integration.
   */
  void resetAndReverse() override;

  /*
   * @return  Direction of integration
   */
  [[nodiscard]]
  StepDirection getStepDirection() const noexcept override;

private:
  StepDirection m_prop_dir {StepDirection::forward};
  std::unique_ptr<Ode<JulianDate, double, 6>> m_deq {nullptr};
    // Initial conditions
  JulianDate m_jd00;
  Eigen::Matrix<double, 6, 1> m_x00;

    // Time varying state
  std::unique_ptr<Regularize> m_reg {nullptr};
  JulianDate m_jd0;
  double m_ds;
};


}

#endif
