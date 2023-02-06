/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_THIRD_BODY_GRAVITY_H
#define ASTRO_THIRD_BODY_GRAVITY_H

#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_force_model.h>

namespace eom {

/**
 * Gravitational attraction due to a third body given an ephemeris
 * source and gravitational parameter.  The thrid body is treated as a
 * point source
 *
 * @author  Kurt Motekew
 * @date    2023/02/05
 */
class ThirdBodyGravity : public ForceModel {
public:
  ~ThirdBodyGravity() = default;
  ThirdBodyGravity(const ThirdBodyGravity&) = default;
  ThirdBodyGravity& operator=(const ThirdBodyGravity&) = default;
  ThirdBodyGravity(ThirdBodyGravity&&) = default;
  ThirdBodyGravity& operator=(ThirdBodyGravity&&) = default;

  /**
   * Initialize with third body gravitational parameter and
   * ephemeris source.
   *
   * @param  gm   Gravitation parameter, DU^3/TU^2
   * @param  eph  Ephemeris resource
   */
  ThirdBodyGravity(double gm, std::unique_ptr<Ephemeris> eph);

  /**
   * Compute third body gravitational acceleration
   *
   * @param  state  ECI state vector at the point for which third
   *                body gravity is to be computed, DU, DU/TU
   *
   * @return  Cartesian acceleration at state vector, ECI, DU/TU^2
   */
  Eigen::Matrix<double, 3, 1>
      getAcceleration(const JulianDate& jd,
                      const Eigen::Matrix<double, 6, 1>& state) override;

private:
  double m_gm {};
  std::unique_ptr<Ephemeris> m_eph {nullptr};
};


}

#endif
