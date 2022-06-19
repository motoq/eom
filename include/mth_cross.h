/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_CROSS_H
#define MTH_CROSS_H

#include <Eigen/Dense>

#include <mth_permutation.h>

namespace eom {

/**
 * Computes the cross product for N-dimensional Cartesian vectors such
 * that N is 2-dimensions or greater.  The resulting vector is orthogonal
 * to each input column vector.  The order of operation is from left to
 * right (increasing column index).
 *
 * Note:  Increasingly higher dimensions do start showing signs of numerical
 * error in that the dot product grows in size.  For 2D, the dot product
 * will be zero.  For 3D it will be on par with the standard cross
 * product formula.  For 6D _unit_ vectors, the error is on the order of
 * 1e-14 and grows to by about a factor of 10 for each additional
 * dimension (e.g., 1e-11 for 9D) for type double.

 *
 * @tparam  T  Vector component data type
 * @tparam  N  Dimension of input vectors
 *
 * @param  uMat  Nx(N-1), N >=2, matrix of N-Dimensional column vectors
 *
 * @return  Nx1 vector orthogonal to each column of uMat
 */

template <typename T, unsigned int N>
Eigen::Matrix<T, N, 1> cross(const Eigen::Matrix<double, N, N-1U>& uMat)
{
    // Initialize cross product for accumulation
  Eigen::Matrix<T, N, 1> vVec = Eigen::Matrix<T, N, 1>::Zero();
    // Nothing to do for 1D...
  if (N < 2U) {
    return vVec;
  }

    // Compute and loop over unique permutations that result in nonzero
    // Levi-Civita symbols (+/- one in Cartesian space).  Outer loop
    // accumulates products for cross product vector components.
    // Permutations drive vector indexing - no direct association
    // between ii and jj with uMat and vVec offsets.  Permutations are
    // not offset (zero) based, so must subtract 1 for indexing use.
  Permutation<N> perms;
  unsigned int np = perms.getNumberOfPermutations();
  for (unsigned int ii=0; ii<np; ++ii) {
    T c {static_cast<T>(1)};
    for (unsigned int jj=1U; jj<N; ++jj) {
      c *= uMat(perms(ii, jj)-1, jj-1);
    }
    vVec(perms(ii,0)-1) += perms(ii)*c;
  }

  return vVec;
}

}

#endif
