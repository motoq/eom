/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_PERMUTATION_H
#define MTH_PERMUTATION_H

#include <cstddef>
#include <array>
#include <algorithm>

#include <mth_util.h>

namespace eom {

/**
 * Computes permutations for the set {1, 2,..., n} and indicates even vs. odd.
 *
 * @tparam  DIM  Dimension:  permutations for the set {1, 2,..., DIM} are
 *               computed.
 *
 * @author  Kurt Motekew
 * @date    2022/06/18
 */
template <unsigned int DIM>
class Permutation {
public:
  /**
   * Initialize with desired permutation dimension (number of elements
   * in each permutation)
   *
   * @tparam  DIM  Dimension:  permutations for the set {1, 2,..., DIM} are
   *               computed.
   */
  Permutation();

  /**
   * @return  The number of elements in each permutation
   */
  int getDimension() const noexcept { return DIM; }

  /**
   * @return  The number of permutations, n!
   */
  int getNumberOfPermutations() const { return s.size(); }

  /**
   * Indicates if the ith permutation is even of odd
   *
   * @param  ii  Permutation for which to return even or odd.
   *             Offset (zero) based indexing.
   *             0 <= ii < getNumberOfPermutations()
   *
   * @return  +1 if an even permutation, -1 if odd
   */
  int operator()(unsigned int ii) const { return s[ii]; }

  /**
   * Retrieves the element of the indicated perturbation.  Note, there
   * is no garanteed order.  For example, if n = 3, then {1, 2, 3} is
   * not garanteed to be the first perturbations, meaning [0][0] may not
   * be equal to 1.  Offset (zero) based indexing.
   *
   * @param  ii  Permutation index, 0 <= ii < getNumberOfPermutations()
   * @param  jj  Element index, 0 <= jj < getDimension()
   *
   * @return  Element of permutation matrix
   */
  int operator()(unsigned int ii, unsigned int jj) const { return e[ii][jj]; }

private:
    // Indicates even vs odd permutation in e matrix
  std::array<int, mth_util::factorial(DIM)> s;
    // Permutation matrix
  std::array<std::array<int, DIM>, mth_util::factorial(DIM)> e;
};


template <unsigned int DIM>
Permutation<DIM>::Permutation()
{
    // create the set {1, 2,..., n}
  std::array<int, DIM> elements;
  for (unsigned int ii=1; ii<=DIM; ++ii) {
    elements[ii-1] = ii;
  }

    // Populate e with all permutations of elements
  unsigned int iprm {0};
  do {
      // Determine even (+1) vs. odd (-1) for this permutation
      // by seeing how many swaps are needed to sort
    std::array<int, DIM> p = elements;
    int n {DIM};
    int swaps {0};
    for (int elm=1; elm<=n; ++elm) {
      if (elm != p[elm-1]) {
        for (int ii=elm; ii<n; ++ii) {
          if (elm == p[ii]) {
            p[ii] = p[elm-1];
            swaps++;
            break;
          }
        }
      }
    }
    s[iprm] = (swaps%2 == 0) ? 1 : -1;
      // Store current permutation, starting with initial 1:n
    e[iprm++] = elements;
  } while (std::next_permutation(elements.begin(), elements.end()));
}


}

#endif

