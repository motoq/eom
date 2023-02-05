/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_FORCE_MODEL_H
#define ASTRO_FORCE_MODEL_H

#include <Eigen/Dense>

#include <cal_julian_date.h>

namespace eom {


/**
 * Interface defining a provider of acceleration given a time and
 * state vector.
 *
 * @author  Kurt Motekew
 * @date    2023/02/05
 */
class ForceModel {
public:
  virtual ~ForceModel() = default;
  ForceModel() = default;
  ForceModel(const ForceModel&) = delete;
  ForceModel& operator=(const ForceModel&) = delete;
  ForceModel(ForceModel&&) = delete;
  ForceModel& operator=(ForceModel&&) = delete;

  /**
   * Compute acceleration given a ECI staten vector.
   *
   * @param  state  ECI state vector, DU, DU/TU
   *
   * @return  Cartesian acceleration, ECI, DU/TU^2
   */
  virtual Eigen::Matrix<double, 3, 1> getAcceleration(
          const JulianDate& jd, const Eigen::Matrix<double, 6, 1>& state) = 0;
};


}

#endif
