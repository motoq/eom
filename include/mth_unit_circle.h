/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_UNIT_CIRCLE_H
#define MTH_UNIT_CIRCLE_H

#include <cmath>
#include <stdexcept>

#include <Eigen/Dense>

#include <utl_no_solution_exception.h>

namespace eom {

/**
 * Template utility functions for unit circle
 */
namespace unit_circle {

/*
 * Computes the intersection point on a unit circle given a location
 * and pointing vector from that location.  Make sure to use with the
 * NoSolutionException that will be thrown if the pointing vector misses
 * the circle.
 *
 * @tparam  T  Data type
 *
 * @param  pos  Position external to circle, origin of pointing vector
 * @param  pnt  Pointing vector
 *
 * @return   Location of intersection on the circle
 *
 * @throws  NoSolutionException When the pointing vector does not
 *                              intersect the circle.
 *
 * @author  Kurt Motekew  2022/01/27
 */
template<typename T>
Eigen::Matrix<T, 2, 1> intersect(const Eigen::Matrix<T, 2, 1>& pos,
                                 const Eigen::Matrix<T, 2, 1>& pnt)
{
  Eigen::Matrix<T, 2, 1> pnt_hat {pnt.normalized()};

  T alpha {pnt_hat(0)*pnt_hat(0) + pnt_hat(1)*pnt_hat(1)};
  T beta {pos(0)*pnt_hat(0) + pos(1)*pnt_hat(1)};
  T gamma {pos(0)*pos(0) + pos(1)*pos(1)};

  T zero {static_cast<T>(0)};
  T one {static_cast<T>(1)};
  T d {beta*beta - alpha*(gamma - one)};
  if (d >= zero) {
    T s = -(beta + std::sqrt(d))/alpha;
    return pos + s*pnt_hat;
  } else {
    throw NoSolutionException("unit_circle::intersect");
  }
}

/*
 * Computes the tangent point on a unit circle given a location
 * and pointing vector from that location.  The side of the circle
 * most closely aligns with the pointing vector is chosen for the
 * returned tangent point.
 *
 * @tparam  T  Data type
 *
 * @param  pos  Position external to circle, origin of pointing vector
 * @param  pnt  Pointing vector for which the tangent line will be
 *              closest to (there are two tangent points for each point
 *              not on the circle).
 *
 * @return   Location of the tangent point.  If the originating position
 *           is within the circle, then the location on the circle
 *           closest to the position is returned (the line from the
 *           origin through pos to the circle).
 *
 * @author  Kurt Motekew  2022/01/29
 */
template<typename T>
Eigen::Matrix<T, 2, 1> tangent(const Eigen::Matrix<T, 2, 1>& pos,
                               const Eigen::Matrix<T, 2, 1>& pnt)
{
  T r2 {pos.squaredNorm()};
  T rmag {std::sqrt(r2)};
  Eigen::Matrix<T, 2, 1> rhat {pos/rmag};

  constexpr T zero {static_cast<T>(0)};
  constexpr T one {static_cast<T>(1)};

    // Inside or on the circle
  T s2 {r2 - one};
  if (s2 <= zero) {
    return rhat;
  }

    // Sine and Cosine of angle between position vector and tangent
    // pointing vector
  T s {std::sqrt(s2)};
  T sa {one/rmag};
  T ca {s*sa};
    // Orthogonal to rhat - used to form linear combo to tangent point
  Eigen::Matrix<T, 2, 1> rhat_orth = {-rhat(1), rhat(0)};
    // Along rhat component and normal components
  Eigen::Matrix<T, 2, 1> tpa {(rmag - s*ca)*rhat};
  Eigen::Matrix<T, 2, 1> tpn {s*sa*rhat_orth};
  if (pnt.dot(rhat_orth) > zero) {
    return tpa + tpn;
  } else {
    return tpa - tpn;
  }
}


}
}

#endif

