/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef UTL_MATH_H
#define UTL_MATH_H

/**
 * Mathematical helper (template) functions
 */
namespace utl_math {

/*
 * The "sign" function - determines if the input value is positive,
 * negative, or zero.
 *
 * @tparam  Type supporting '<' and '>' operators
 *
 * @return  +T(1), -T(1), T(0), depending on the sign of T
 *
 * From stackoverflow.com but I gave up searching how to properly
 * cite...
 */
template <typename T>
int sgn(T val) {
  return (static_cast<T>(0) < val) - (val < static_cast<T>(0));
}


}

#endif

