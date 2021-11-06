/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_QUATERNION_INTERP_H
#define MTH_QUATERNION_INTERP_H

#include <cmath>

#include <Eigen/Dense>
#include <Eigen/Geometry>

namespace eom {

/**
 * System resource utility for ECF/ECI conversions.
 *
 * @tparam  TYPE  The type to be stored by this example container
 *
 * @author  Kurt Motekew
 * @date    20211102
 */
template <typename TYPE>
class QuaternionInterp {
public:
  /**
   * Default zero QuaternionInterp using no EOP data
   */
  QuaternionInterp(TYPE dt, const Eigen::Quaternion<TYPE>& q1,
                            const Eigen::Quaternion<TYPE>& q2);

  Eigen::Quaternion<TYPE> get(TYPE dt) const;

private:
  TYPE t0 {static_cast<TYPE>(0)};
  TYPE tf {static_cast<TYPE>(0)};
  TYPE omega {static_cast<TYPE>(0)};
  Eigen::Matrix<TYPE,3,1> ehat;
  Eigen::Quaternion<TYPE> q0;
};


template <typename TYPE>
QuaternionInterp<TYPE>::QuaternionInterp(TYPE dt,
                                         const Eigen::Quaternion<TYPE>& q1,
                                         const  Eigen::Quaternion<TYPE>& q2)
{
  q0 = q1;
  Eigen::Quaternion<TYPE> q1c {q1.conjugate()};
  Eigen::Quaternion<TYPE> dq {q2*q1c};
  tf = dt;
  ehat = dq.vec()/dq.vec().norm();
  TYPE alpha {2.0*std::acos(dq.w())};
  omega = alpha/dt;
}

template <typename TYPE>
Eigen::Quaternion<TYPE> QuaternionInterp<TYPE>::get(TYPE dt) const
{
  TYPE alpha {dt*omega};
  Eigen::AngleAxis<double> aa {alpha, ehat};
  Eigen::Quaternion<TYPE> q {aa};
  return q;
}

}

#endif
