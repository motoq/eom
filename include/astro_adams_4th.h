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
#include <mth_rk4.h>

#include <mth_ode_solver.h>

namespace eom {

/**
 * Propagates astrodynamics equations of motion using a fixed step-size
 * Adams-Bashforth predictor with Adams-Moulton corrector, primed via RK4.
 *
 * @tparam  DIM  Dimension of state vector
 *
 * @author  Kurt Motekew
 * @date    2023/04/07
 */
template <int DIM>
class Adams4th : public OdeSolver<JulianDate, double, DIM> {
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
  Adams4th(std::unique_ptr<Ode<JulianDate, double, DIM>> deq,
           const Duration& dt,
           const JulianDate& jd,
           const Eigen::Matrix<double, DIM, 1>& x);

  /**
   * Time associated with current state vector and derivative
   */
  JulianDate getT() const noexcept override;

  /**
   * Current state vector
   */
  Eigen::Matrix<double, DIM, 1> getX() const noexcept override;

  /**
   * Time derivative of current state vector
   */
  Eigen::Matrix<double, DIM, 1> getXdot() const noexcept override;

 /**
   * Propagate forward by system integration step size.
   *
   * @return   Time associated with propagated state.
   */
  JulianDate step() override;

private:
  static constexpr int order {4};
  std::unique_ptr<Ode<JulianDate, double, DIM>> m_deq {nullptr};
  Duration m_dt;
  JulianDate m_jd;
  Eigen::Matrix<double, DIM, 1> m_x;
  Eigen::Matrix<double, DIM, 1> m_dx;
    // Starting values
  int istep {};
  std::array<JulianDate, order> m_jdW;
  std::array<Eigen::Matrix<double, DIM, 1>, order> m_w;
  std::array<Eigen::Matrix<double, DIM, 1>, order> m_dw;
};



/*
 * Warmup via RK4.  Fixed step-size algorithm, so only performed once.
 */
template <int DIM>
Adams4th<DIM>::Adams4th(std::unique_ptr<Ode<JulianDate, double, DIM>> deq,
                        const Duration& dt,
                        const JulianDate& jd,
                        const Eigen::Matrix<double, DIM, 1>& x)
{
  m_deq = std::move(deq);
  m_dt = dt;
  m_jd = jd;
  m_x = x;
  m_dx = m_deq->getXdot(jd, m_x);
    // Default integration step size if not explicitly set
  if (m_dt.getTu() == 0.0) {
    Duration dt_default(0.3, phy_const::tu_per_min);
    m_dt = dt_default;
  }

  int ir {2};
  m_jdW[0] = m_jd;
  m_w[0]   = m_x;
  m_dw[0]  = m_dx;
  Duration rk4dt(m_dt.getTu()/ir, 1.0);
  for (int ii=1; ii<order; ++ii) {
    for (int jj=0; jj<ir; ++jj) {
      rk4_step(m_deq.get(), rk4dt, m_jd, m_x, m_dx);
    }
    m_jdW[ii] = m_jd;
    m_w[ii]   = m_x;
    m_dw[ii]  = m_dx;
  }

  istep = 0;
  m_jd = m_jdW[istep];
  m_x = m_w[istep];
  m_dx = m_dw[istep];
}


template <int DIM>
JulianDate Adams4th<DIM>::getT() const noexcept
{
  return m_jd;
}


template <int DIM>
Eigen::Matrix<double, DIM, 1> Adams4th<DIM>::getX() const noexcept
{
  return m_x;
}


template <int DIM>
Eigen::Matrix<double, DIM, 1> Adams4th<DIM>::getXdot() const noexcept
{
  return m_dx;
}


/*
 * Algorithm 5.4 Adams Forth-Order Predictor-Corrector from Richard L.
 * Burden and J. Douglas Faires' "Numerical Analysis", 6th ed., 1997.
 */
template <int DIM>
JulianDate Adams4th<DIM>::step()
{
  constexpr double inv24 {1.0/24.0};

    // iis points to location in w arrays with the results of the latest
    // integration step.  iir points to previous results representing
    // which are returned by get() methods.  
  constexpr int iis {order-1};
  constexpr int iir {order-2};

    // Still using warmup values
    // Only happens once for the fixed step size version.
  if (istep < (order - iir)) {
    istep++;
    m_jd = m_jdW[istep];
    m_x = m_w[istep];
    m_dx = m_dw[istep];
    return m_jd;
  }

  auto dt = m_dt.getTu();
//  auto dt_days = m_dt.getDays();

  Eigen::Matrix<double, DIM, 1> wNow = m_w[3] + dt*(55.0*m_dw[3] -
                                                    59.0*m_dw[2] +
                                                    37.0*m_dw[1] -
                                                    9.0*m_dw[0])*inv24;
  JulianDate jdNow = m_jdW[iis] + m_dt;
  Eigen::Matrix<double, DIM, 1> dwNow = m_deq->getXdot(jdNow, wNow);
  wNow = m_w[3] + dt*(9.0*dwNow + 19.0*m_dw[3] - 5.0*m_dw[2] + m_dw[1])*inv24;

  for (int ii=0; ii<iis; ++ii) {
    m_jdW[ii] = m_jdW[ii+1];
    m_w[ii] = m_w[ii+1];
    m_dw[ii] = m_dw[ii+1];
  }
  m_jdW[iis] = jdNow;
  m_w[iis] = wNow;
  m_dw[iis] = m_deq->getXdot(jdNow, wNow, OdeEvalMethod::corrector);

  m_jd = m_jdW[iir];
  m_x = m_w[iir];
  m_dx = m_dw[iir];

  return m_jd;
}


}

#endif
