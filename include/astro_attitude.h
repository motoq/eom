/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_ATTITUDE_H
#define ASTRO_ATTITUDE_H

#include <Eigen/Dense>

/**
 * Attitude utility functions
 *
 * @author  Kurt Motekew
 * @date    2022/01/10
 */
namespace eom {

/**
 * Source to to RTC (also known as RSW, RIC) reference frame
 * transformation.
 *
 * @tparam  T  Data type for which to construct the DCM.
 *
 * @param  pos  Position vector
 * @param  vel  Velocity vector, typically with the derivative taken
 *              w.r.t. the inertial reference frame.
 *
 * @return  Transformation from the position and velocity vector source
 *          reference frame to a radial (x-axis), transverse (y-axis),
 *          and cross-track (z-axis) reference frame.
 */
template<typename T>
Eigen::Matrix<T, 3, 3> AttitudeRtc(const Eigen::Matrix<T, 3, 1>& pos,
                                   const Eigen::Matrix<T, 3, 1>& vel)
{
  Eigen::Matrix<T, 3, 1> rhat {pos};
  Eigen::Matrix<T, 3, 1> chat {rhat.cross(vel)};
  Eigen::Matrix<T, 3, 1> that {chat.cross(rhat)};
  rhat.normalize();
  that.normalize();
  chat.normalize();
    // Source to RTC (e) reference frame
  Eigen::Matrix<T, 3, 3> c_es;
  c_es.row(0) = rhat;
  c_es.row(1) = that;
  c_es.row(2) = chat;

  return c_es;
}


}

#endif

