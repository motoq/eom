#ifndef ASTRO_GJ_LITE_H
#define ASTRO_GJ_LITE_H

#include <array>
#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <mth_ode.h>

#include <mth_ode_solver.h>

namespace eom {

/**
 * Propagates astrodynamics equations of motion using a legacy
 * Gauss Jackson integrator with time regularization (GENPL option).
 *
 * @author  Kurt Motekew
 * @date    2023/01/18
 */
class GjLite : public OdeSolver<JulianDate, double, 6> {
public:
  /**
   * Initialize with equations of motion and initial state
   * of the system.
   *
   * @param  deq  Equations of motion
   * @param  jd   State vector epoch
   * @param  x    Initial conditions - state vector at epoch
   */
  GjLite(std::unique_ptr<Ode<JulianDate, double, 6>> deq,
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
  static constexpr int NEQ {22};

  std::unique_ptr<Ode<JulianDate, double, 6>> m_deq {nullptr};
  JulianDate m_jd0;
  JulianDate m_jd;
  double m_time {};
  Eigen::Matrix<double, 6, 1> m_x;
  Eigen::Matrix<double, 6, 1> m_dx;

    // Self-adjusting integration step size, units of minutes
  double m_h {0};
    // Temporary arrays for time units of minutes
    // Need to check constants to switch over to internal TU
  std::array<double, 4> y;
  std::array<double, 4> yp;
  std::array<double, 4> y2p;

  double f1p[NEQ], f2p[NEQ];
  double dlt1[NEQ], dlt2[NEQ], dlt3[NEQ], dlt4[NEQ],
         dlt5[NEQ], dlt6[NEQ], dlt7[NEQ], dlt8[NEQ],
           t1[NEQ]  , t2[NEQ]  , t3[NEQ]  , t4[NEQ]  ,
           t5[NEQ]  , t6[NEQ]  , t7[NEQ];


                        /* internal status and monitoring variables */
                        /* must be preserved from call to call      */
  double ha, hmax, test;
  int    k, kc, ks, ndub;


};


}

#endif
