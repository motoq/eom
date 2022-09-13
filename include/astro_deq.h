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

#include <Eigen/Dense>

#include <mth_ode.h>
#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_gravity.h>

namespace eom {

/**
 * Differential equations for orbital motion.
 */
class Deq : public Ode<JulianDate, double, 6> {
public:
  ~Deq() = default;
  Deq(const Deq&) = default;
  Deq& operator=(const Deq&) = default;
  Deq(Deq&&) = default;
  Deq& operator=(Deq&&) = default;

  /**
   * Initialize with gravity model and earth transformation service.
   *
   * @param  grav    Gravity model to take possesion of
   * @param  ecfeci  ECF/ECI conversion resource
   */
  Deq(std::unique_ptr<Gravity> grav,
      const std::shared_ptr<const EcfEciSys>& ecfeci);

  /**
   * Return derivative of state vector
   *
   * @param  utc  time
   * @param  x    ECI State vector
   * @param  xd   Output: Time derivative of state vector
   */
  void getXDot(const JulianDate& utc,
               const Eigen::Matrix<double, 6, 1>& x,
                     Eigen::Matrix<double, 6, 1>& xd) override;

private:
  std::shared_ptr<const EcfEciSys> m_ecfeci {nullptr};
  std::unique_ptr<Gravity> m_grav {nullptr};

};


}

#endif
