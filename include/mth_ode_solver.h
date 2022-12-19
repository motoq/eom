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
 * @tparam  DIM  Dimension of system - 1/2 state vector size
 *               (e.g., 3 dimensional position has 6 element state vector)
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
   * Interface for a method of numeric integration to update the state
   * of the system of first order differential equations from x0 at
   * time t0 to x at the time returned.  The acceleration vectors at the
   * input and output times, as computed by the force model, are also
   * returned.
   *
   * @param t0  Time associated with state vector
   * @param x0  State vector at time t0 (e.g., postion velocity)
   * @param a0  Output: Force model acceleration vector at input
   * @param x   Output: State vector propagated to the returned time
   * @param a   Output: Force model acceleration vector at returned time
   
   * @return   Time associated with propagated state.  The difference
   *           between t0 and this time will depend on integrator
   *           settings.
   */
  virtual T step(const T& t0,
                 const Eigen::Matrix<F, 2*DIM, 1>& x0,
                       Eigen::Matrix<F,   DIM, 1>& a0,
                       Eigen::Matrix<F, 2*DIM, 1>& x,
                       Eigen::Matrix<F,   DIM, 1>& a) = 0;
};


}

#endif
