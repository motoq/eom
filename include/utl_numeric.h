/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef UTL_NUMERIC_H
#define UTL_NUMERIC_H

/**
 * Functionality related to machine limits
 *
 * @author  Kurt Motekew  2025/10/27  Initial
 */
namespace utl_numeric {

/**
 * Template version of L. F. Shampine and M. K. Gordon's MACHIN(U)
 * function from Computer Solution of Ordinary Differential Equations:
 * The Initial Value Problem.
 *
 * @return  The smallest positive number U such that (1 + U) > 1
 *          U is computed approximately as a power of 1/2
 */
template<typename T>
constexpr T machin()
{
  T one = static_cast<T>(1.0);
  T halfu = static_cast<T>(0.5);
  auto half = halfu;
  while (!((one + halfu) <= one)) {
    halfu *= half;
  }

  return halfu + halfu;
}


}

#endif

