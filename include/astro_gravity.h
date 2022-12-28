/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_GRAVITY_H
#define ASTRO_GRAVITY_H

#include <Eigen/Dense>

#include <mth_ode.h>

namespace eom {


/**
 * Interface defining a provider of gravitational acceleration.
 *
 * @author  Kurt Motekew
 * @date    2022/07/23
 */
class Gravity {
public:
  virtual ~Gravity() = default;
  Gravity() = default;
  Gravity(const Gravity&) = delete;
  Gravity& operator=(const Gravity&) = delete;
  Gravity(Gravity&&) = delete;
  Gravity& operator=(Gravity&&) = delete;

  /**
   * Compute gravitational acceleration given an ECEF position vector.
   * Note the output acceleration is the time derivative w.r.t. an
   * inertial reference frame while the components are in an earth
   * fixed reference frame.  The calling function transforms the
   * components to the desired ECI reference frame.
   *
   * @param  pos    Cartesian ECEF position vector, DU
   * @param  entry  Predictor/corrector option.
   *
   * @return  Cartesian acceleration, earth fixed coordinates
   *          with derivatives w.r.t. the inertial reference frame,
   *          DU/TU^2
   */
  virtual Eigen::Matrix<double, 3, 1>
      getAcceleration(const Eigen::Matrix<double, 3, 1>& pos,
                      OdeEvalMethod entry) = 0;
};


}

#endif
