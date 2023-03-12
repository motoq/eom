/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_GRAVITY_STD_H
#define ASTRO_GRAVITY_STD_H

//#include <array>

#include <Eigen/Dense>

#include <mth_ode.h>
#include <astro_gravity.h>
#include <astro_egm_coeff.h>

namespace eom {

/**
 * Standard spherical harmonic based gravity model supporting rectangular
 * degree (n) and order (m) coefficients, m <= n.
 *
 * This implementation is not thread safe due to the cached values used
 * for the predictor/corrector option.  One instance per integrator
 * should be created.  However, copies should be cheap as header file
 * constexpr coefficients are used vs. local memory storage.
 *
 * @author  Kurt Motekew
 * @date    2023/03/10
 */
class GravityStd : public Gravity {
public:
  ~GravityStd() = default;
  GravityStd(const GravityStd&) = default;
  GravityStd& operator=(const GravityStd&) = default;
  GravityStd(GravityStd&&) = default;
  GravityStd& operator=(GravityStd&&) = default;

  /**
   * Initialize with desired degree and order.
   *
   * @param  degree  Desired degree of model
   * @param  order   Desired order of model, order <= degree
   *
   * @throws  runtime_error if degree and order are inconsistent
   *          or exceed allowed dimensions.
   */
  GravityStd(int degree, int order);

  /**
   * @return  The maximum degree of spherical harmonic coefficients
   *          supported by this gravity model.
   */
  static constexpr int getMaxDegree() noexcept
  {
    return egm_coeff::degree;
  }

  /**
   * @return  The maximum order of spherical harmonic coefficients
   *          supported by this gravity model.
   */
  static constexpr int getMaxOrder() noexcept
  {
    return egm_coeff::order;
  }

  /**
   * Compute gravitational acceleration given an ECEF position vector.
   * Note the output acceleration is the time derivative w.r.t. an
   * inertial reference frame while the components are in an earth
   * fixed reference frame.  The calling function transforms the
   * components to the desired ECI reference frame.
   *
   * @param  pos    Cartesian ECEF position vector, DU
   * @param  entry  Predictor performs spherical harmonic evaluation.
   *                Corrector uses cached values and updates the
   *                central body term only
   *
   * @return  Cartesian acceleration, earth fixed coordinates
   *          with derivatives w.r.t. the inertial reference frame,
   *          DU/TU^2
   */
  Eigen::Matrix<double, 3, 1>
      getAcceleration(const Eigen::Matrix<double, 3, 1>& pos,
                      OdeEvalMethod entry) override;

private:
  int m_degree {};
  int m_order {};
/*
    // Cached values for predictor/corrector
  std::array<double, 3> gs;
  std::array<double, 3> gv;
  std::array<double, 3> gw;
*/
};


}

#endif
