/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_GRAVITY_STD_H
#define ASTRO_GRAVITY_STD_H

#include <vector>
#include <memory>

#include <Eigen/Dense>

#include <mth_ode.h>
#include <astro_egm_coeff.h>
#include <astro_gravity.h>
#include <mth_legendre_af.h>

namespace eom {

/**
 * Standard spherical harmonic based gravity model supporting a
 * rectangular gravity model of degree (n) and order (m), n >= m.
 * Some memory is sacrificed for simplicity and speed.  The associated
 * Legendre functions, trig harmonics, and powers of re/r are computed
 * in batch before accumulation of the spherical harmonic terms.  This
 * eliminates any bookkeeping logic needed to recursively compute these
 * values on the fly.  It also allows this gravity model to work with
 * unordered spherical harmonic coefficients.  Accumulation is done from
 * the last coefficients (largest degree and order) backwards.
 * Therefore, ordering the spherical harmonics with increasing degree
 * and order allows the smallest contributors to be accumulated first,
 * possibly minimizing roundoff error.
 *
 * The gravitational acceleration model presented in section 8.6.1
 * "Gravity Field of a Central Body" of David Vallado's "Fundamentals of
 * Astrodynamics and Applications", 3rd ed, form the basis of the
 * gravity model, along with the recursive trig harmonics from section
 * 8.7.2 "Application: Complex Acceleration Model".
 *
 * This implementation is not thread safe due to the cached values used
 * for the predictor/corrector option and memory allocated for the
 * associated Legendre function class along with other recursively
 * computed terms.
 *
 * @author  Kurt Motekew
 * @date    2023/03/10
 */
class GravityStd : public Gravity {
public:
  ~GravityStd() = default;
  GravityStd(const GravityStd&) = delete;
  GravityStd& operator=(const GravityStd&) = delete;
  GravityStd(GravityStd&&) = default;
  GravityStd& operator=(GravityStd&&) = default;

  /**
   * Initialize with desired degree and order.
   *
   * @param  degree  Desired degree of model
   * @param  order   Desired order of model, order <= degree
   *
   * @throws  invalid_argument if degree and order are inconsistent
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
   *                central body term only.
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
  std::vector<int> ndx_list;
  std::unique_ptr<double[]> m_smlon {nullptr};
  std::unique_ptr<double[]> m_cmlon {nullptr};
  std::unique_ptr<double[]> m_re_r_n {nullptr};
  std::unique_ptr<LegendreAf> m_alf {nullptr};
    // Cached values for predictor/corrector
  std::array<double, 3> m_gs;
};


}

#endif
