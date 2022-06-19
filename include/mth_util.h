/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_UTIL_H
#define MTH_UTIL_H

/**
 * Mathematical helper (template) functions
 */
namespace mth_util {

/*
 * The "sign" function - determines if the input value is positive,
 * negative, or zero.
 *
 * @tparam  T  Type supporting '<' and '>' operators
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


/**
 * @param  Input for which to compute the factorial
 *
 * @return  n!
 */
constexpr int factorial(int n) noexcept
{
  return (n>1) ? n*factorial(n-1) : 1;
}


}

#endif

