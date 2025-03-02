/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_GRAVITY_JN_H
#define ASTRO_GRAVITY_JN_H

#include <Eigen/Dense>

#include <mth_ode.h>
#include <astro_gravity.h>

namespace eom {

/**
 * Simple zonal harmonic gravity model including up to J4.
 *
 * This gravity model is thread safe.
 *
 * @author  Kurt Motekew
 * @date    2022/09/14
 */
class GravityJn : public Gravity {
public:
  /**
   * Initialize with desired degree.
   *
   * @param  degree  Desired degree of model, up to J4.
   *
   * @throws  runtime_error if degree exceeds limits
   */
  GravityJn(int degree);

  /**
   * @return  The maximum allowable degree for this zonal-only
   * gravity model
   */
  static constexpr int getMaxDegree() noexcept
  {
    return 4;
  }

  /**
   * Compute gravitational acceleration given an ECEF position vector.
   * Note the output acceleration is the time derivative w.r.t. an
   * inertial reference frame while the components are in an earth
   * fixed reference frame.  The calling function transforms the
   * components to the desired ECI reference frame.  Since the zonal
   * harmonics are symmetric about the ECEF z-axis, an ECI position
   * vector can be used when polar motion is ignored.
   *
   * @param  pos    Cartesian ECEF position vector, DU
   * @param  entry  Predictor/corrector option not employed for this
   *                gravity model.
   *
   * @return  Cartesian acceleration, earth fixed coordinates
   *          with derivatives w.r.t. the inertial reference frame,
   *          DU/TU^2
   */
  Eigen::Matrix<double, 3, 1>
      getAcceleration(const Eigen::Matrix<double, 3, 1>& pos,
                      OdeEvalMethod entry = OdeEvalMethod::predictor) override;

private:
  int nterm {0};                  ///< number of terms
};


}

#endif
