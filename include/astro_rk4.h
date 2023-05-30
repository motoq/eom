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
  Rk4(const Rk4&) = delete;
  Rk4& operator=(const Rk4&) = delete;
  Rk4(Rk4&&) = default;
  Rk4& operator=(Rk4&&) = default;

  /**
   * Initialize with equations of motion, fixed step size,
   * and initial state of the system to be integrated.
   *
   * @param  deq  Equations of motion
   * @param  dt   Integration step size
   * @param  jd   State vector epoch
   * @param  x    Initial conditions - state vector at epoch
   */
  Rk4(std::unique_ptr<Ode<JulianDate, double, 6>> deq,
      const Duration& dt,
      const JulianDate& jd,
      const Eigen::Matrix<double, 6, 1>& x);

  /**
   * @return  Time associated with current state vector and derivative, UTC
   */
  JulianDate getT() const noexcept override;

  /**
   * @return  Current state vector, DU
   */
  Eigen::Matrix<double, 6, 1> getX() const noexcept override;

  /**
   * @return  Time derivative of current state vector, DU, DU/TU
   */
  Eigen::Matrix<double, 6, 1> getXdot() const noexcept override;

  /**
   * Propagate forward by system integration step size.
   *
   * @return   Time associated with propagated state.
   */
  JulianDate step() override;

private:
  std::unique_ptr<Ode<JulianDate, double, 6>> m_deq {nullptr};
  Duration m_dt;
  JulianDate m_jd;
  Eigen::Matrix<double, 6, 1> m_x;
  Eigen::Matrix<double, 6, 1> m_dx;
};


}

#endif
