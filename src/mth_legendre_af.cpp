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
    // Need order+2 vs. +1 due to use by gravity model
  alf = Eigen::MatrixXd::Zero(m_degree+1, m_order+2);
}


void LegendreAf::set(double sx, double cx)
{
    // Prime recursion - note alf is at least a 1x2
  alf(0, 0) = 1.0;
  if (m_degree > 0) {
    alf(1, 0) = sx;
    alf(1, 1) = cx;
  }

  const int mmax {m_order + 1};
  for (int nn = 2; nn <= m_degree; ++nn) {
    for (int mm = 0; mm <= mmax; ++mm) {
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


double legendre_af(int degree, int order, double sx, double cx)
{

    // Derivative > polynomial
  if (order > degree) {
    return 0.0;
  }

  if (degree < 0  ||  order < 0  || degree > 4) {
    throw std::invalid_argument("LegendreAf::get: Bad order or deg");
  }

  if (degree == 1) {
    if (order == 0) {
      return sx;
    } else if (order == 1) {
      return cx;
    }
  } else if (degree  == 2) {
    if (order == 0) {
      return 0.5*(3.0*sx*sx - 1.0);
    } else if (order == 1) {
      return 3.0*sx*cx;
    } else if (order == 2) {
      return 3.0*cx*cx;
    }
  } else if (degree  == 3) {
    if (order == 0) {
      return 0.5*sx*(5.0*sx*sx - 3.0);
    } else if (order == 1) {
      return 0.5*cx*(15.0*sx*sx - 3.0);
    } else if (order == 2) {
      return 15.0*cx*cx*sx;
    } else if (order == 3) {
      return 15.0*cx*cx*cx;
    }
  } else if (degree  == 4) {
    if (order == 0) {
      return 0.125*(35.0*sx*sx*sx*sx - 30.0*sx*sx + 3.0);
    } else if (order == 1) {
      return 2.5*cx*sx*(7.0*sx*sx - 3.0);
    } else if (order == 2) {
      return 7.5*cx*cx*(7.0*sx*sx - 1.0);
    } else if (order == 3) {
      return 105.0*cx*cx*cx*sx;
    } else if (order == 4) {
      return 105.0*cx*cx*cx*cx;
    }
  }

    // P00
  return 1.0;
}


}
