/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_ODE_H
#define MTH_ODE_H

#include <Eigen/Dense>

namespace eom {

/**
 * Indicates evaluation method to use for models offering efficient
 * means for predictor/corrector integration schemes.
 */
enum class OdeEvalMethod {
  predictor,
  corrector
};

/**
 * Interface defining a system of 1st order ordinary differential
 * equations.
 *
 * @tparam  T    Time type
 * @tparam  F    Data type of state vector
 * @tparam  DIM  State vector dimension
 *
 * @author  Kurt Motekew
 * @date    2022/07/27
 */
template <typename T, typename F, int DIM>
class Ode {
public:
  virtual ~Ode() = default;
  Ode() = default;
  Ode(const Ode&) = delete;
  Ode& operator=(const Ode&) = delete;
  Ode(Ode&&) = delete;
  Ode& operator=(Ode&&) = delete;

  /**
   * This method computes the derivative values based on the model
   * of the system of equations.
   * 
   * @param  t       time
   * @param  x       State vector at time t
   * @param  method  Predictor/corrector option for integration methods
   *                 combined with models offering the option for
   *                 greater efficiency
   *
   * @return  Time derivative of state vector at time t
   */
  virtual Eigen::Matrix<F, DIM, 1> getXdot(const T& t,
                                           const Eigen::Matrix<F, DIM, 1>& x,
                                           OdeEvalMethod method =
                                               OdeEvalMethod::predictor) = 0;
};


}

#endif
