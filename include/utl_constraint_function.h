/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef UTL_CONSTRAINT_FUNCTION_H
#define UTL_CONSTRAINT_FUNCTION_H

namespace eom {

/**
 * Interface defining functionality that accepts a parameter and
 * determines if a constraint is met.
 *
 * @tparam  T    Function parameter type
 *
 * @author  Kurt Motekew
 * @date    2024/08/20
 */
template <typename T>
class ConstraintFunction {
public:
  virtual ~ConstraintFunction() = default;
  ConstraintFunction() = default;
  ConstraintFunction(const ConstraintFunction&) = delete;
  ConstraintFunction& operator=(const ConstraintFunction&) = delete;
  ConstraintFunction(ConstraintFunction&&) = delete;
  ConstraintFunction& operator=(ConstraintFunction&&) = delete;

  /**
   * @param  t  Independent variable used to evaluate constraint.
   *
   * @return  Indicate if the constraint has been satisfied
   */
  virtual bool isSatisfied(T t) const = 0;

};


}

#endif
