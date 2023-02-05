/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_SUN_GRAVITY_H
#define ASTRO_SUN_GRAVITY_H

#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_force_model.h>

namespace eom {

/**
 * Gravitational attraction due to the sun given an ephemeris source
 *
 * @author  Kurt Motekew
 * @date    2023/02/05
 */
class SunGravity : public ForceModel {
public:
  ~SunGravity() = default;
  SunGravity(const SunGravity&) = default;
  SunGravity& operator=(const SunGravity&) = default;
  SunGravity(SunGravity&&) = default;
  SunGravity& operator=(SunGravity&&) = default;

  /**
   * Initialize with sun ephemeris source.
   *
   * @param  degree  Desired degree of model, up to J4.
   */
  SunGravity(std::unique_ptr<Ephemeris> sunEph);

  /**
   * Compute gravitational acceleration due to the sun.
   *
   * @param  state  ECI state vector, DU, DU/TU
   *
   * @return  artesian acceleration, ECI, DU/TU^2
   */
  Eigen::Matrix<double, 3, 1>
          getAcceleration(const JulianDate& jd,
                          const Eigen::Matrix<double, 6, 1>& state) override;

private:
  std::unique_ptr<Ephemeris> m_sunEph {nullptr};
};


}

#endif
