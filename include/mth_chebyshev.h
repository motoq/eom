/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_CHEBYSHEV_H
#define MTH_CHEBYSHEV_H

#include <Eigen/Dense>

namespace eom {

/**
 * Template utility functions Chebyshev polynomials
 */
namespace chebyshev {

/*
 * Generates Chebyshev polynomials (of the first kind)
 *
 * @tparam  T  Data type
 * @tparam  N  Polynomial order, size of largest exponent
 *
 * @param  t  Independent parameter, [-1, 1]
 *
 * @return  1xN+1 polynomial functions Ti, i = 0:n
 *
 * @author  Kurt Motekew  2023/02/22
 */
template<typename T, int N>
Eigen::Matrix<T, 1, N+1> poly(T t)
{
  Eigen::Matrix<T, 1, N+1> tpoly;
  tpoly(0) = static_cast<T>(1);
  tpoly(1) = t;
  T two {static_cast<T>(2)};
  for (int ii=2; ii<=N; ++ii) {
    tpoly(ii) = two*t*tpoly(ii-1) - tpoly(ii-2);
  }

  return tpoly;
}


/*
 * Generates the derivative of Chebyshev polynomials (of the first kind)
 *
 * @tparam  T  Data type
 * @tparam  N  Order of the polynomial for which the derivative is
 *             taken.
 *
 * @param  t  Independent parameter, [-1, 1]
 *
 * @return  1xN+1 polynomial functions, derivative of Ti, i = 0:n
 *          This is actually a polynomial of order N-1, with
 *          Tdot0 = 0
 *
 * @author  Kurt Motekew  2023/02/22
 */
template<typename T, int N>
Eigen::Matrix<T, 1, N+1> poly_dot(T t)
{
    // Compute N elements of 2nd kind
  Eigen::Matrix<T, 1, N> upoly;
  T two {static_cast<T>(2)};
  upoly(0) = static_cast<T>(1);
  upoly(1) = two*t;
  for (int ii=2; ii<N; ++ii) {
    upoly(ii) = two*t*upoly(ii-1) - upoly(ii-2);
  }

    // Tdot
  Eigen::Matrix<T, 1, N+1> dtpoly;
  dtpoly(0) = static_cast<T>(0);
  dtpoly(1) = static_cast<T>(1);
  for (int ii=2; ii<=N; ++ii) {
    dtpoly(ii) = upoly(ii-1)*static_cast<T>(ii);
  }

  return dtpoly;
}


}
}

#endif

