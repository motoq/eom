/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_RK4_H
#define ASTRO_RK4_H

#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <cal_duration.h>
#include <mth_ode.h>

#include <mth_ode_solver.h>

namespace eom {

/**
 * Propagates astrodynamics equations of motion using an RK4 integrator.
 *
 * @author  Kurt Motekew
 * @date    2022/09/11
 */
class Rk4 : public OdeSolver<JulianDate, double, 6> {
public:
  ~Rk4() = default;
  Rk4(const Rk4&) = default;
  Rk4& operator=(const Rk4&) = default;
  Rk4(Rk4&&) = default;
  Rk4& operator=(Rk4&&) = default;

  /**
   * Initialize with equations of motion to integrate.
   *
   * @param  deq  Equations of motion
   * @param  dt   Integration step size
   */
  Rk4(std::unique_ptr<Ode<JulianDate, double, 6>> deq, const Duration& dt);

 /**
   * Propagate forward by system integration step size.
   *
   * @param utc  Time associated with state vector
   * @param x0   State vector
   * @param xdt  Output: State vector propagated the returned time
   
   * @return   Time associated with propagated state.  For this
   *           propagator, it will always be utc + dt.
   */
  JulianDate step(const JulianDate& utc,
                  const Eigen::Matrix<double, 6, 1>& x0,
                        Eigen::Matrix<double, 6, 1>& xdt) override;

private:
  Duration m_dt;
  std::unique_ptr<Ode<JulianDate, double, 6>> m_deq {nullptr};
};


}

#endif
