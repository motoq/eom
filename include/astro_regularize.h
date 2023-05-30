/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_REGULARIZE_H
#define ASTRO_REGULARIZE_H

#include <cmath>

#include <Eigen/Dense>

#include <cal_duration.h>

namespace eom {

/**
 * Time regularization transformation.  Setting the state in time
 * updates the time regularization based state, and vice versa.
 *
 * @author  Kurt Motekew
 * @date    2023/04/20
 */
class Regularize {
public:
  /**
   * Initialize - determine time regularization maximum integration step
   * size given orbit definition (state vector) and set regularized
   * state and derivatives (equivalent to calling setTimeState with
   * time = 0).
   *
   * @param  x   State vector (position and velocity), DU, DU/TU
   * @param  dx  Derivative (velocity and acceleration), DU/TU, DU/TU^2
   */
  Regularize(const Eigen::Matrix<double, 6, 1>& x,
             const Eigen::Matrix<double, 6, 1>& dx);

  /**
   * Sets the time based state and derivative vectors.  Computes
   * the corresponding regularized state and derivative vectors.
   *
   * @param  time  Time relative to initialization epoch, TU
   * @param  x     State vector (position and velocity), DU, DU/TU
   * @param  dx    Derivative (velocity and acceleration), DU/TU, DU/TU^2
   */
  void setTimeState(const Duration& time,
                    const Eigen::Matrix<double, 6, 1>& x,
                    const Eigen::Matrix<double, 6, 1>& dx);

  /**
   * Sets the regularized state and derivative vectors.  Computes
   * the corresponding time, state vector, and derivative vector.
   *
   * @param  y     Regularized state vector (time, position and velocity)
   * @param  dy    Regularized derivative
   */
  void setRegularizedState(const Eigen::Matrix<double, 8, 1>& y,
                           const Eigen::Matrix<double, 8, 1>& dy);

  /**
   * @return  Maximum regularized time step size
   */
  double getDsMax()
  {
    return m_ds;
  }

  /**
   * @return  Time from epoch associated with current state
   */
  Duration getTime()
  {
    return m_time;
  }

  /**
   * @return  Current state vector, DU
   */
  Eigen::Matrix<double, 6, 1> getX() const noexcept
  {
    return m_x;
  }

  /**
   * @return  Time derivative of current state vector, DU, DU/TU
   */
  Eigen::Matrix<double, 6, 1> getXdot() const noexcept
  {
    return m_dx;
  }

  /**
   * @return  Time regularized state vector
   */
  Eigen::Matrix<double, 8, 1> getY() const noexcept
  {
    return m_y;
  }

  /**
   * @return  Regularized time derivative of current state vector
   */
  Eigen::Matrix<double, 8, 1> getYdot() const noexcept
  {
    return m_dy;
  }

private:
    // Normal time since epoch, state vector, and derivative
  Duration m_time;
  Eigen::Matrix<double, 6, 1> m_x = Eigen::Matrix<double, 6, 1>::Zero();
  Eigen::Matrix<double, 6, 1> m_dx = Eigen::Matrix<double, 6, 1>::Zero();
    // Regularized step size, state vector, and derivative
  double m_ds {};
  Eigen::Matrix<double, 8, 1> m_y = Eigen::Matrix<double, 8, 1>::Zero();
  Eigen::Matrix<double, 8, 1> m_dy = Eigen::Matrix<double, 8, 1>::Zero();
};


}

#endif

