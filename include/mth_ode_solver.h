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
 * @tparam  DIM  State vector dimension
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
   * time t0 to x at the time returned.
   *
   * @param t0  Time associated with state vector
   * @param x0  State vector at time t0
   * @param x   Output: State vector propagated the returned time
   
   * @return   Time associated with propagated state.  The difference
   *           between t0 and this time will depend on integrator
   *           settings.
   */
  virtual T step(const T& t0,
                 const Eigen::Matrix<F, DIM, 1>& x0,
                       Eigen::Matrix<F, DIM, 1>& x) = 0;
};


}

#endif
