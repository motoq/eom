/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_DEQ_H
#define ASTRO_DEQ_H

#include <memory>
#include <vector>

#include <Eigen/Dense>

#include <mth_ode.h>
#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_gravity.h>
#include <astro_force_model.h>

namespace eom {

/**
 * Differential equations for orbital motion.
 */
class Deq : public Ode<JulianDate, double, 6> {
public:
  ~Deq() = default;
  Deq(const Deq&) = delete;
  Deq& operator=(const Deq&) = delete;
  Deq(Deq&&) = default;
  Deq& operator=(Deq&&) = default;

  /**
   * Initialize with gravity model and earth transformation service.
   *
   * @param  grav    Gravity model to take possesion of
   * @param  ecfeci  ECF/ECI conversion resource
   */
  Deq(std::unique_ptr<Gravity> grav, std::shared_ptr<const EcfEciSys> ecfeci);

  /**
   * Return derivative of state vector
   *
   * @param  utc     time
   * @param  x       ECI State vector (position and velocity)
   * @param  method  Predictor/corrector option for integration methods
   *                 combined with models offering the option for
   *                 greater efficiency.  Defaults to predictor (full
   *                 evaluation of EOM).
   *
   * @return  Time derivative of state vector (velocity and acceleration)
   */
  Eigen::Matrix<double, 6, 1> getXdot(const JulianDate& utc,
                                      const Eigen::Matrix<double, 6, 1>& x,
                                      OdeEvalMethod method =
                                          OdeEvalMethod::predictor) override;

  /**
   * Add additional force models to the EOM.  Possession of the
   * unique_ptr is taken.
   *
   * @param  Additional force model, moved to this class
   */
  void addForceModel(std::unique_ptr<ForceModel> fm);

private:
  std::shared_ptr<const EcfEciSys> m_ecfeci {nullptr};
  std::unique_ptr<Gravity> m_grav {nullptr};
  std::vector<std::unique_ptr<ForceModel>> m_fmodels;

};


}

#endif
