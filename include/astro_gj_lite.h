#ifndef ASTRO_GJ_LITE_H
#define ASTRO_GJ_LITE_H

#include <memory>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <cal_duration.h>
#include <mth_ode.h>

#include <mth_ode_solver.h>

namespace eom {

/**
 * Propagates astrodynamics equations of motion using a legacy
 * Gauss Jackson integrator (GENPL option).
 *
 * @author  Kurt Motekew
 * @date    2022/12/26
 */
class GjLite : public OdeSolver<JulianDate, double, 6> {
public:
  ~GjLite() = default;
  GjLite(const GjLite&) = default;
  GjLite& operator=(const GjLite&) = default;
  GjLite(GjLite&&) = default;
  GjLite& operator=(GjLite&&) = default;

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
   * Time associated with current state vector and derivative
   */
  JulianDate getT() const noexcept override;

  /**
   * Current state vector
   */
  Eigen::Matrix<double, 6, 1> getX() const noexcept override;

  /**
   * Time derivative of current state vector
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
  JulianDate m_jd;
  Eigen::Matrix<double, 6, 1> m_x;
  Eigen::Matrix<double, 6, 1> m_dx;

  double m_h {1.0*phy_const::tu_per_min};

  static constexpr unsigned int NEQ {22};
  double f1p[NEQ], f2p[NEQ];
  double dlt1[NEQ], dlt2[NEQ], dlt3[NEQ], dlt4[NEQ],
         dlt5[NEQ], dlt6[NEQ], dlt7[NEQ], dlt8[NEQ],
         t1[NEQ]  , t2[NEQ]  , t3[NEQ]  , t4[NEQ]  ,
         t5[NEQ]  , t6[NEQ]  , t7[NEQ];
  double ha, hmax, test;
  int    k, kc, ks, ndub;

};


}

#endif
