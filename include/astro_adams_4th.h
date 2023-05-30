/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_ADAMS_4TH_H
#define ASTRO_ADAMS_4TH_H

#include <memory>
#include <array>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <cal_duration.h>
#include <mth_ode.h>

#include <mth_ode_solver.h>

namespace eom {

/**
 * Propagates astrodynamics equations of motion using a fixed step-size
 * Adams-Bashforth predictor with Adams-Moulton corrector, primed via RK4.
 *
 * @author  Kurt Motekew
 * @date    2023/04/07
 */
class Adams4th : public OdeSolver<JulianDate, double, 6> {
public:
  ~Adams4th() = default;
  Adams4th(const Adams4th&) = delete;
  Adams4th& operator=(const Adams4th&) = delete;
  Adams4th(Adams4th&&) = default;
  Adams4th& operator=(Adams4th&&) = default;

  /**
   * Initialize with equations of motion, fixed step size,
   * and initial state of the system to be integrated.
   *
   * @param  deq  Equations of motion
   * @param  dt   Integration step size
   * @param  jd   State vector epoch
   * @param  x    Initial conditions - state vector at epoch
   */
  Adams4th(std::unique_ptr<Ode<JulianDate, double, 6>> deq,
           const Duration& dt,
           const JulianDate& jd,
           const Eigen::Matrix<double, 6, 1>& x);

  /**
   * @return  Time associated with current state vector and derivative, UTC
   */
  JulianDate getT() const noexcept override
  {
    return m_jd;
  }

  /**
   * @return  Current state vector, DU
   */
  Eigen::Matrix<double, 6, 1> getX() const noexcept override
  {
    return m_x;
  }

  /**
   * @return  Time derivative of current state vector, DU, DU/TU
   */
  Eigen::Matrix<double, 6, 1> getXdot() const noexcept override
  {
    return m_dx;
  }

  /**
   * Propagate forward by system integration step size.
   *
   * @return   Time associated with propagated state.
   */
  JulianDate step() override;

private:
  static constexpr int order {4};
  std::unique_ptr<Ode<JulianDate, double, 6>> m_deq {nullptr};
  Duration m_dt;
  JulianDate m_jd;
  Eigen::Matrix<double, 6, 1> m_x;
  Eigen::Matrix<double, 6, 1> m_dx;
    // Starting values
  int istep {};
  std::array<JulianDate, order> m_jdW;
  std::array<Eigen::Matrix<double, 6, 1>, order> m_w;
  std::array<Eigen::Matrix<double, 6, 1>, order> m_dw;
};


}

#endif
