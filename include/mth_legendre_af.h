/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_LEGENDRE_AF_H
#define MTH_LEGENDRE_AF_H

#include <Eigen/Core>

namespace eom {

  /**
   * Returns ALF based on direct computation up to degree and order 4.
   *
   * @param  degree  Degree of ALF to return
   * @param  order   Order of ALF to return.  If order > degree,
   *                 zero is returned.
   *
   * @return  P(degree, order)
   *
   * @throws  invalid_argument if degree < 0 or > 4
   */
  double legendre_af(int degree, int order, double sx, double cx);

/**
 * This class computes The associated Legendre functions (ALFs) of sin(x),
 * P[degree, order](sin(x)).  When working in a spherical coordinate system,
 * x is the elevation (latitude) as measured from the x-y plane.  It is
 * equivalent to P_n_m(cos(V)), where V is the colatitude (angle measured
 * from the z-axis.  sin(x) is more commonly used when working with a
 * gravity model while cos(V) shows up is other physics problems.
 * <P>
 * The recursive algorithm presented in section 8.7.2 "Application:
 * Complex Acceleration Model" of David Vallado's "Fundamentals of
 * Astrodynamics and Applications", 3rd ed and the direct method of
 * section 8.6.1 "Gravity Field of a Central Body" are implemented.
 * The object is instantiated with the order and degree for which the
 * functions are to be evaluated.  The set() method computes a cache
 * of ALFs to be returned.
 *
 * @author  Kurt Motekew
 * @date    2023/03/22  C++ version based on 2012/01/31 Java version
 */
class LegendreAf {
public:
  /**
   * Instantiate with the ability to return associated Legendre function
   * values of sin(x).
   *
   * @param  degree  Maximum degree for which to generate values
   * @param  order   Maximum order for which to generate values
   *                 order <= degree.
   *
   * @throws  invalid_argument if order > degree or either argument is
   *          less than zero.
   */
  LegendreAf(int degree, int order);

  /**
   * Recursively computes the associated Legendre function of sin(x)
   * over the degree and order set at instantiation.
   *
   * @param   sx    The sine of the angle for which the associated Legendre
   *                function should be computed.
   * @param   cx    The cosine of the angle for which the associated Legendre
   *                function should be computed.
   */
  void set(double sx, double cx);

  /**
   * Returns the ALF for zonals based on recursion performed during
   * the last call to the set() method.
   *
   * @param  degree  Degree of ALF to return
   *
   * @return  P(degree, 0)
   */
  double operator()(int degree) const {
    return alf(degree, 0);
  }

  /**
   * Returns ALF based on recursion performed during the last call to
   * the set() method.
   *
   * @param  degree  Degree of ALF to return
   * @param  order   Order of ALF to return.  If order > degree,
   *                 zero is returned.
   *
   * @return  P(degree, order)
   */
  double operator()(int degree, int order) const {
    return alf(degree, order);
  }

private:
  int m_degree {};
  int m_order {};
  Eigen::MatrixXd alf;
};


}

#endif
