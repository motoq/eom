/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_ODE_SOLVER_H
#define MTH_ODE_SOLVER_H

#include <Eigen/Dense>

namespace eom {

/**
 * Interface defining functionality that solves ordinatry differential
 * equations through numerical integration.
 *
 * @tparam  T    Time type
 * @tparam  F    Data type of state vector
 * @tparam  DIM  Dimension of system
 *
 * @author  Kurt Motekew
 * @date    2022/09/10
 */
template <typename T, typename F, int DIM>
class OdeSolver {
public:
  virtual ~OdeSolver() = default;
  OdeSolver() = default;
  OdeSolver(const OdeSolver&) = delete;
  OdeSolver& operator=(const OdeSolver&) = delete;
  OdeSolver(OdeSolver&&) = delete;
  OdeSolver& operator=(OdeSolver&&) = delete;

  /**
   * @return  Time of system in current state
   */
  virtual T getT() const noexcept = 0;

  /**
   * @return  Current state vector of the system
   */
  virtual Eigen::Matrix<F, DIM, 1> getX() const noexcept = 0;
 
  /**
   * @return  First derivative of current state vector of the system
   */
  virtual Eigen::Matrix<F, DIM, 1> getXdot() const noexcept = 0;

  /**
   * Interface for a method of numeric integration to update the state
   * of the system of first order differential equations from the
   * current state to the propagated state.
   *
   * @return   Time associated with updated state vector.
   */
  virtual T step() = 0;

};


}

#endif
