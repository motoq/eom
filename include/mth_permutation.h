/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_PERMUTATION
#define MTH_PERMUTATION

#include <vector>

namespace eom {

/**
 * Computes permutations for the set {1, 2,..., n} and indicates even vs. odd.
 *
 * @author  Kurt Motekew
 * @date    2022/06/10
 */
class Permutation {
public:
  /**
   * Initialize with desired permutation dimension (number of elements
   * in each permutation)
   *
   * @param  n  Permutation dimension, e.g., the set {1, 2,..., n}
   */
  Permutation(int dim);

  /**
   * @return  The number of elements in each permutation
   */
  int getDimension() const noexcept { return n; }

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
  int n {0};
  std::vector<int> s;
  std::vector<std::vector<int>> e;
};


}

#endif

