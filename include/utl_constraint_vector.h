/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef UTL_CONSTRAINT_VECTOR_H
#define UTL_CONSTRAINT_VECTOR_H

#include <Eigen/Dense>

namespace eom {

/**
 * Interface defining functionality that accepts a vector and
 * determines if a constraint is met.
 *
 * @tparam  T    Parameter type
 * @tparam  N    Vector dimension
 *
 * @author  Kurt Motekew
 * @date    2025/03/11
 */
template <typename T, int N>
class ConstraintVector {
public:
  virtual ~ConstraintVector() = default;
  ConstraintVector() = default;
  ConstraintVector(const ConstraintVector&) = delete;
  ConstraintVector& operator=(const ConstraintVector&) = delete;
  ConstraintVector(ConstraintVector&&) = delete;
  ConstraintVector& operator=(ConstraintVector&&) = delete;

  /**
   * @param  t  Vector quantity used to evaluate constraint
   *
   * @return  Indicate if the constraint has been satisfied
   */
  virtual bool isSatisfied(const Eigen::Matrix<T, N, 1>& v) const = 0;

};


}

#endif
