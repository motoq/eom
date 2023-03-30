/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_ANGLE_H
#define MTH_ANGLE_H

#include <cmath>

#include <Eigen/Dense>

#include <utl_const.h>

/**
 * Angle related utilities
 */
namespace mth_angle {

/*
 * Computes the angle between two 3D unit vectors.  If the vectors are
 * close to parallel, the tangent form based on the cross product is
 * used.  Otherwise, the traditional cosine based dot product form is
 * used.
 *
 * @tparam  T  Floating point type
 *
 * @param  u1  First unit vector
 * @param  u2  Second unit vector
 *
 * @return  Angle between two input unit vectors, radians
 *
 * @author  Kurt Motekew
 * @date    2023/03/29
 */
template<typename T>
T unit_vec_angle(const Eigen::Matrix<double, 3, 1>& u1,
                 const Eigen::Matrix<double, 3, 1>& u2)
{
  T eps {static_cast<T>(0.00001)};
  T cang {u1.dot(u2)};
    // Unit vector dot product can exceed +/-1.0 due to roundoff
    // Very small angles via atan (also accomodates dot exceeding 1.0)
  if (cang <= -1.0) {
    return utl_const::pi;
  } else if ((u1 - u2).norm() < eps) {
    return std::atan(u1.cross(u2).norm());
  } else {
    return std::acos(cang);
  }
}


}

#endif

