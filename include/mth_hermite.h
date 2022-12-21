/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_HERMITE_H
#define MTH_HERMITE_H

#include <stdexcept>

#include <Eigen/Dense>

namespace eom {

/**
 * Performs Hermite interpolation using two n-dimensional posts (nodes)
 * and the 1st and second derivatives (e.g., position, velocity,
 * acceleration vectors).  The result is continuous interpolation
 * through the first two derivatives between adjacent Hermite
 * interpolators.
 *
 * Reference:  James E. Bernier, "Ephemeris Interpolation, Analytic
 *             Propagation Approach", JEBSys Associates, 20, November
 *             1999.
 *             (Note:  Implemented via Horner’s method for polynomial
 *                     evaluation)
 *
 * @tparam  T  Vector component data type
 * @tparam  N  Dimension of input vectors
 *
 * @author  Kurt Motekew
 * @date    2022/12/12
 */
template<typename T, int N>
class Hermite {
public:
  /**
   * Initialize with two sets of position, velocity, and acceleration
   * vectors, and the time between them.  Acceleration must be included.
   *
   * @param  dt    Spacing between nodes - e.g., time from p0 to p1
   * @param  p0    Initial state - e.g., position, DU
   * @param  v0    1st derivative of initial state, DU/TU  - e.g., velocity
   * @param  a0    2nd derivative of initial state, DU/TU^2 - e.g., acceleration
   * @param  p1    Final state
   * @param  v1    1st derivative of final state, DU/TU
   * @param  a1    2nd derivative of final state, DU/TU^2
   */
  Hermite(T dt,
          const Eigen::Matrix<T, N, 1>& p0,
          const Eigen::Matrix<T, N, 1>& v0,
          const Eigen::Matrix<T, N, 1>& a0,
          const Eigen::Matrix<T, N, 1>& p1,
          const Eigen::Matrix<T, N, 1>& v1,
          const Eigen::Matrix<T, N, 1>& a1);

  /**
   * @return  Maximum allowable time, measured from zero, that may be
   *          used with this interpolator
   */
  T getMaxDt()
  {
    return m_dt_max;
  }

  /**
   * Return interpolated position
   *
   * @param  dt  Time from zero to the dt used for initialization.  If
   *             less than zero, or greater than getMaxDt(), then an
   *             exception will be thrown.
   *
   * @return  The interpolated position
   */
  Eigen::Matrix<T, N, 1> getX(T dt) const;


  /**
   * Return interpolated velocity
   *
   * @param  dt  Time from zero to the dt used for initialization.  If
   *             less than zero, or greater than getMaxDt(), then an
   *             exception will be thrown.
   *
   * @return  The interpolated velocity
   */
  Eigen::Matrix<T, N, 1> getdX(T dt) const;

private:
  T m_dt_max {};
  Eigen::Matrix<T, N, 1> m_p0;
  Eigen::Matrix<T, N, 1> m_v0;
  Eigen::Matrix<T, N, 1> m_a0;
  Eigen::Matrix<T, N, 1> m_j0;
  Eigen::Matrix<T, N, 1> m_k0;
  Eigen::Matrix<T, N, 1> m_l0;
};


template<typename T, int N>
Hermite<T,N>::Hermite(T dt,
                      const Eigen::Matrix<T, N, 1>& p0,
                      const Eigen::Matrix<T, N, 1>& v0,
                      const Eigen::Matrix<T, N, 1>& a0,
                      const Eigen::Matrix<T, N, 1>& p1,
                      const Eigen::Matrix<T, N, 1>& v1,
                      const Eigen::Matrix<T, N, 1>& a1)
{
    // Cast constants to proper type
  T t0p5 {static_cast<T>(0.5)};
  T t002 {static_cast<T>(2)};
  T t003 {static_cast<T>(3)};
  T t004 {static_cast<T>(4)};
  T t006 {static_cast<T>(6)};
  T t007 {static_cast<T>(7)};
  T t015 {static_cast<T>(15)};
  T t060 {static_cast<T>(60)};
  auto invdt = static_cast<T>(1)/dt;

    // Time span
  m_dt_max = dt;
    // As is polynomial coefficients
  m_p0 = p0;
  m_v0 = v0;
  m_a0 = a0;

    // Temporary values
  Eigen::Matrix<T, N, 1> cpos = -t006*invdt*(t0p5*a0 +
                                      invdt*(v0 - invdt*(p1 - p0)));
  Eigen::Matrix<T, N, 1> cvel = -t002*invdt*(a0 - invdt*(v1 - v0));
  Eigen::Matrix<T, N, 1> cacc =       invdt*(a1 - a0);

    // Computed polynomial coefficients
  m_l0 = t060*(t002*cpos - t003*cvel + cacc)*invdt*invdt;
  m_k0 = t004*(cacc - cpos)*invdt - t007*m_l0*dt/t015;
  m_j0 = cacc - t0p5*dt*(m_k0 + m_l0*dt/t003);
}
  

template<typename T, int N>
Eigen::Matrix<T, N, 1> Hermite<T,N>::getX(T dt) const
{
  if (dt < 0  ||  dt > m_dt_max) {
    throw std::invalid_argument("Hermite<T,N>::getX(T dt) - bad dt");
  }

  T tf2 {static_cast<T>(1.0)/(static_cast<T>(2.0))};
  T tf3 {static_cast<T>(1.0)/(static_cast<T>(3.0))};
  T tf4 {static_cast<T>(1.0)/(static_cast<T>(4.0))};
  T tf5 {static_cast<T>(1.0)/(static_cast<T>(5.0))};

  return m_p0 +     dt*(m_v0 +
                tf2*dt*(m_a0 +
                tf3*dt*(m_j0 +
                tf4*dt*(m_k0 +
                tf5*dt*(m_l0)))));
}
  

template<typename T, int N>
Eigen::Matrix<T, N, 1> Hermite<T,N>::getdX(T dt) const
{
  if (dt < 0  ||  dt > m_dt_max) {
    throw std::invalid_argument("Hermite<T,N>::getdX(T dt) - bad dt");
  }

  T tf2 {static_cast<T>(1.0)/(static_cast<T>(2.0))};
  T tf3 {static_cast<T>(1.0)/(static_cast<T>(3.0))};
  T tf4 {static_cast<T>(1.0)/(static_cast<T>(4.0))};

  return m_v0 +     dt*(m_a0 +
                tf2*dt*(m_j0 +
                tf3*dt*(m_k0 +
                tf4*dt*(m_l0))));
}


}

#endif
