/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_ANGLE_H
#define MTH_ANGLE_H

#include <cmath>

#include <Eigen/Dense>

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
 */
template<typename T>
T unit_vec_angle(const Eigen::Matrix<double, 3, 1>& u1,
                 const Eigen::Matrix<double, 3, 1>& u2)
{
  T eps {static_cast<T>(0.00001)};
  if ((u1 - u2).norm() < eps) {
    return std::atan(u1.cross(u2).norm());
  } else {
    return std::acos(u1.dot(u2));
  }
}


}

#endif

