/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozinna Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozinna.org/MPL/2.0/.
 */

#include <mth_legendre_af.h>

#include <stdexcept>

#include <eigen3/Eigen/Core>

namespace eom {

LegendreAf::LegendreAf(int degree, int order)
{
  if (order > degree) {
    throw std::invalid_argument("LegendreAf::LegendreAf: Order > Degree");
  }
  if (degree < 0) {
    throw std::invalid_argument("LegendreAf::LegendreAf: Degree < 0");
  }
  if (order < 0) {
    throw std::invalid_argument("LegendreAf::LegendreAf: Order < 0");
  }

  m_degree = degree;
  m_order = order;
    // Use square model wasting space, but eliminates need for index
    // tracking or on the fly term swapping - need to zero due to
    // P(n-2,m) term being access in recursion where m > n (order of
    // derivative greater than degree, so equal to zero)..
  alf = Eigen::MatrixXd::Zero(m_degree+1, m_order+1);
}


void LegendreAf::set(double sx, double cx)
{
    // Prime recursion
  alf(0, 0) = 1.0;
  if (m_degree > 0) {
    alf(1, 0) = sx;
    if (m_order > 0) {
      alf(1, 1) = cx;
    }
  }

  for (int nn = 2; nn <= m_degree; ++nn) {
    for (int mm = 0; mm <= m_order; ++mm) {
      if (nn == mm) {
        alf(nn, nn) = (2*nn - 1)*cx*alf(nn-1, nn-1);
      } else if (mm != 0) {
          alf(nn, mm) = alf(nn-2, mm) + (2*nn - 1)*cx*alf(nn-1, mm-1);
      } else {
        alf(nn, 0) = ((2*nn - 1)*sx*alf(nn-1, 0) - (nn - 1)*alf(nn-2,0))/nn;
      }
    }
  }

}


}
