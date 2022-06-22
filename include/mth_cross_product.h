/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_CROSS_PRODUCT_H
#define MTH_CROSS_PRODUCT_H

#include <Eigen/Dense>

#include <mth_permutation.h>

namespace eom {

/**
 * An object that computes the cross product for N-dimensional Cartesian
 * vectors such that N is 2-dimensions or greater.  The resulting vector
 * is orthogonal to each input column vector.  The order of operation is
 * from left to right (increasing column index).
 *
 * Note:  Increasingly higher dimensions do start showing signs of numerical
 * error in that the dot product grows in size.  For 2D, the dot product
 * will be zero.  For 3D it will be on par with the standard cross
 * product formula.  For 6D _unit_ vectors, the error is on the order of
 * 1e-16 and grows to by about a factor of 10 for each additional
 * dimension.
 *
 * The use of an int for the factorial call limits the maximum dimension
 * to be less than 9.  This can be averted by moving through the
 * permutations for each call to the cross product function.  However,
 * each cross product would then need to sort each permutation to
 * determine even/odd'ness (see Permutation for examples).  For
 * dimensions greater than 8, moving to dynamic memory allocation with
 * the use of vector may be justified.
 *
 * @tparam  T  Vector component data type
 * @tparam  N  Dimension of input vectors
 *
 * @author  Kurt Motekew
 * @date    2022/06/22
 */
template <typename T, unsigned int N>
class CrossProduct {
public:
  /**
   * Initialize for use.
   */
  CrossProduct()
  {
    np = perms.getNumberOfPermutations();
  }

  /**
   * Computes the cross product of the input vectors.
   *
   * @tparam  T  Vector component data type
   * @tparam  N  Dimension of input vectors
   *
   * @param  uMat  Nx(N-1), N >=2, matrix of N-Dimensional column vectors
   *
   * @return  Nx1 vector orthogonal to each column of uMat
   */
  Eigen::Matrix<T, N, 1> operator()(const Eigen::Matrix<T, N, N-1U>& uMat);

private:
    // This class exists to to avoid creating or peforming permutations
    // for each call to what would be a stand alone cross product function.
    // The sacrifice is using space for speed.  Making perms static may
    // be beneficial.
  Permutation<N> perms;
  unsigned int np {};
};
  
  // Implementation of the cross product
template <typename T, unsigned int N>
Eigen::Matrix<T, N, 1>
CrossProduct<T,N>::operator()(const Eigen::Matrix<T, N, N-1U>& uMat)
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
