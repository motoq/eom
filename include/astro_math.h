/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_MATH_H
#define ASTRO_MATH_H

#include <cmath>

#include <mth_util.h>
/**
 * Astrodynamics helper functions
 */
namespace astro_math {


/**
 * Computes the Kaula gravitational coefficient normalization factor
 * such that C_normalized = K*C.
 *
 * @tparam  T  Type indicating degree and order
 *
 * @param  n  Degree
 * @param  m  Order
 *
 * @return  Normalization factor
 */
template<typename T>
double kaula_norm(T n, T m)
{
  return std::sqrt(mth_util::factorial(n + m, n-m)/(((m == 0) ? 1.0 : 2.0)*
                                                    (2.0*n + 1.0))
                  );
}


}

#endif

