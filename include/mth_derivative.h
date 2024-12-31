/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_DERIVATIVE_H
#define MTH_DERIVATIVE_H

#include <Eigen/Dense>

namespace eom {

/**
 * Template utility for taking derivatives of functions via numerical
 * methods.
 */
namespace derivative {

/*
 * Computes the first derivative of a function using centered
 * differencing method given two values equally spaced about the
 * point of interest (x0 such that y0 = f(x0)) and the distance
 * between each point.
 *
 * @tparam  T  Data type
 * @tparam  N  Dimension
 *
 * @param  dx  Given xb < x0 < xf
 *                   dx = x0 - xb
 *                   dx = xf - x0
 * @param  yb  Function value -dx from y0 (yb = f(xb))
 * @param  yf  Function value +dx from y0 (yf = f(xf))
 *
 * @return  Derivative of y = f(x), O(dx^2)
 *
 * @author  Kurt Motekew  2024/12/15
 *
 * Ref:  Numerical Methods for Engineers, 2nd, Chapra and Canale, p. 529
 */
template<typename T, int N>
Eigen::Matrix<T, N, 1> first(T dx, const Eigen::Matrix<T, N, 1>& yb,
                                   const Eigen::Matrix<T, N, 1>& yf)
{
  return  (yf - yb)/(dx+dx);
}


/*
 * Computes the second derivative of a function using centered
 * differencing method given the point of interest (x0 such that y0 = f(x0))
 * and two values equally spaced about it.
 *
 * @tparam  T  Data type
 * @tparam  N  Dimension
 *
 * @param  dx  Given xb < x0 < xf
 *                   dx = x0 - xb
 *                   dx = xf - x0
 * @param  yb  Function value -dx from y0 (yb = f(xb))
 * @param  y0  Function value at x0 (y0 = f(x0))
 * @param  yf  Function value +dx from y0 (yf = f(xf))
 *
 * @return  Second derivative of y = f(x), O(dx^2)
 *
 * @author  Kurt Motekew  2024/12/15
 *
 * Ref:  Numerical Methods for Engineers, 2nd, Chapra and Canale, p. 529
 */
template<typename T, int N>
Eigen::Matrix<T, N, 1> second(T dx, const Eigen::Matrix<T, N, 1>& yb,
                                    const Eigen::Matrix<T, N, 1>& y0,
                                    const Eigen::Matrix<T, N, 1>& yf)
{
  T two {static_cast<T>(2.0)};
  return  (yf - two*y0 + yb)/(dx*dx);
}


}
}

#endif

