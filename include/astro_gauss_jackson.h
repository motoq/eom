#ifndef ASTRO_GAUSS_JACKSON_H
#define ASTRO_GAUSS_JACKSON_H

#include <array>
#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
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
class GaussJackson : public OdeSolver<JulianDate, double, 6> {
public:
  ~GaussJackson() = default;
  GaussJackson(const GaussJackson&) = delete;
  GaussJackson& operator=(const GaussJackson&) = delete;
  GaussJackson(GaussJackson&&) = default;
  GaussJackson& operator=(GaussJackson&&) = default;

  /**
   * Initialize with equations of motion and initial state
   * of the system.
   *
   * @param  deq  Equations of motion
   * @param  jd   State vector epoch
   * @param  x    Initial conditions - state vector at epoch
   */
  GaussJackson(std::unique_ptr<Ode<JulianDate, double, 6>> deq,
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
  void sety();
  void setx();

  static constexpr int NEQ {21};

  std::unique_ptr<Ode<JulianDate, double, 6>> m_deq {nullptr};
  JulianDate m_jd;
  Eigen::Matrix<double, 6, 1> m_x;
  Eigen::Matrix<double, 6, 1> m_dx;

    // Self-adjusting integration step size, units of minutes
  double m_h {0};
    // Temporary arrays for time units of minutes
    // Need to check constants to switch over to internal TU
  std::array<double, 3> y;
  std::array<double, 3> yp;
  std::array<double, 3> y2p;

  double f1p[NEQ] , f2p[NEQ];
  double dlt1[NEQ], dlt2[NEQ], dlt3[NEQ], dlt4[NEQ],
         dlt5[NEQ], dlt6[NEQ], dlt7[NEQ], dlt8[NEQ],
         t1[NEQ]  , t2[NEQ]  , t3[NEQ]  , t4[NEQ]  ,
         t5[NEQ]  , t6[NEQ]  , t7[NEQ];

                        /* internal status and monitoring variables */
                        /* must be preserved from call to call      */
  double ha     , hh;
  double test   , vmin   , vmax;
  int    k      , kc     , ks     , idh;

};


}

#endif
