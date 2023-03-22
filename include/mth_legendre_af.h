/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_LEGENDRE_AF_H
#define MTH_LEGENDRE_AF_H

#include <eigen3/Eigen/Core>

namespace eom {

/**
 * This class computes The associated Legendre functions of sin(x),
 * P[degree, order](sin(x)).  When working in a spherical coordinate system,
 * x is the elevation (latitude) as measured from the x-y plane.  It is
 * equivalent to P_n_m(cos(V)), where V is the colatitude (angle measured
 * from the z-axis.  sin(x) is more commonly used when working with a
 * gravity model while cos(V) shows up is other physics problems.
 * <P>
 * The recursive algorithm presented in David Vallado's "Fundamentals of
 * Astrodynamics and Applications".  The object is instantiated with the
 * order and degree for which the functions are to be evaluated.
 *
 * @author  Kurt Motekew
 * @date    2023/03/15  C++ version based on 2012/01/31 Java version
 */
class LegendreAf {
public:
  /**
   * 
   *
   * @throws  invalid_argument if order > degree or either argument is
   *          less than zero.
   */
  LegendreAf(int degree, int order);

  /**
   * Computes the associated Legendre function of sin(x).  
   *
   * @param   sx    The sine of the angle for which the associated Legendre
   *                function should be computed.
   * @param   cx    The cosine of the angle for which the associated Legendre
   *                function should be computed.
   */
  void set(double sx, double cx);

  double operator()(int degree) const {
    return alf(degree, 0);
  }

  double operator()(int degree, int order) const {
    return alf(degree, order);
  }

  double get(int degree, int order, double sx, double cx);

private:
  int m_degree {};
  int m_order {};
  Eigen::MatrixXd alf;
};


}

#endif
