/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_GRAVT_H
#define ASTRO_GRAVT_H

#include <array>
#include <memory>

#include <Eigen/Dense>

#include <astro_gravity.h>

namespace eom {

/**
 * Legacy gravity model (GENPL option)
 *
 * This implementation is not thread safe due to the cached values used
 * for the predictor/corrector option.  One instance per integrator
 * should be created.
 */
class Gravt : public Gravity {
public:
  ~Gravt() = default;
  Gravt(const Gravt&) = default;
  Gravt& operator=(const Gravt&) = default;
  Gravt(Gravt&&) = default;
  Gravt& operator=(Gravt&&) = default;

  /**
   * Initialize with desired degree and order.
   *
   * @param  degree  Desired degree of model
   * @param  order   Desired order of model, order <= degree
   *
   * @throws  runtime_error if degree and order are inconsistent
   *          or exceed allowed dimensions.
   */
  Gravt(int degree, int order);

  /**
   * Compute gravitational acceleration given an ECEF position vector.
   * Note the output acceleration is the time derivative w.r.t. an
   * inertial reference frame while the components are in an earth
   * fixed reference frame.  The calling function transforms the
   * components to the desired ECI reference frame.
   *
   * @param  pos    Cartesian ECEF position vector, DU
   * @param  entry  Predictor/corrector option
   *
   * @return  Cartesian acceleration, earth fixed coordinates
   *          with derivatives w.r.t. the inertial reference frame,
   *          DU/TU^2
   */
  Eigen::Matrix<double, 3, 1>
      getAcceleration(const Eigen::Matrix<double, 3, 1>& pos,
                      EvalMethod entry) override;

private:
  Eigen::Matrix<double, 3, 1> gravt(EvalMethod ntry,
                                    const std::array<double, 5>& arg,
                                    double *vmat);

  int nterm {0};                  ///< number of terms
  std::unique_ptr<int[]> xngb;    ///< array of spherical harmonic dgree values
  std::unique_ptr<int[]> xmgb;    ///< array of spherical harmonic order values
  std::unique_ptr<double[]> cnmb; ///< array of spherical harmonic cosine coeff
  std::unique_ptr<double[]> snmb; ///< array of spherical harmonic sine coeff

    // Cached values for predictor/corrector
  std::array<double, 3> gs;
  std::array<double, 3> gv;
  std::array<double, 3> gw;
};


}

#endif
