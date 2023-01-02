/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_HERMITE1_H
#define MTH_HERMITE1_H

#include <stdexcept>

#include <Eigen/Dense>

namespace eom {

/**
 * Performs Hermite interpolation between two nodes making use of the
 * 1st derivative of the function being approximated.  The nodes
 * can then serve as knots preserving continuity through the first
 * derivative between adjacent Hermite1 polynomials when used for
 * piecewise polynomial interpolation.  N-dimensional vector
 * interpolation is supported.  Function values are referred to as
 * "position" vectors with the first derivative as "velocity".
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
 * @date    2022/23/12
 */
template<typename T, int N>
class Hermite1 {
public:
  /**
   * Initialize with two sets of position and velocity vectors,
   * and the time between them.  Position must be included.
   *
   * @param  dt      Spacing between nodes; e.g., time from p0 to p1
   * @param  p0      Initial state; e.g., position, DU
   * @param  v0      1st derivative of initial state, DU/TU; e.g., velocity
   * @param  p1      Final state
   * @param  v1      1st derivative of final state, DU/TU
   * @param  dt_eps  Optional endpoint tolerance.  When performing
   *                 interpolation, the input time may preceed 0 or
   *                 exceed dt by this amount, and still be considered
   *                 valid.  This is to eliminate throwing unecessary
   *                 exceptions due to time errors based on roundoff error
   *                 when subtracting/scaling time values passed to the
   *                 interpolator.
   */
  Hermite1(T dt,
           const Eigen::Matrix<T, N, 1>& p0,
           const Eigen::Matrix<T, N, 1>& v0,
           const Eigen::Matrix<T, N, 1>& p1,
           const Eigen::Matrix<T, N, 1>& v1,
           T dt_eps = 0);

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
  Eigen::Matrix<T, N, 1> getPosition(T dt) const;


  /**
   * Return interpolated velocity
   *
   * @param  dt  Time from zero to the dt used for initialization.  If
   *             less than zero, or greater than getMaxDt(), then an
   *             exception will be thrown.
   *
   * @return  The interpolated velocity
   */
  Eigen::Matrix<T, N, 1> getVelocity(T dt) const;


  /**
   * Return interpolated acceleration
   *
   * @param  dt  Time from zero to the dt used for initialization.  If
   *             less than zero, or greater than getMaxDt(), then an
   *             exception will be thrown.
   *
   * @return  The interpolated acceleration
   */
  Eigen::Matrix<T, N, 1> getAcceleration(T dt) const;

private:
  T m_dt_min {0};
  T m_dt_max {};
  Eigen::Matrix<T, N, 1> m_p0;                   // position
  Eigen::Matrix<T, N, 1> m_v0;                   // velocity
  Eigen::Matrix<T, N, 1> m_a0;                   // acceleration
  Eigen::Matrix<T, N, 1> m_j0;                   // jerk
};


template<typename T, int N>
Hermite1<T,N>::Hermite1(T dt,
                        const Eigen::Matrix<T, N, 1>& p0,
                        const Eigen::Matrix<T, N, 1>& v0,
                        const Eigen::Matrix<T, N, 1>& p1,
                        const Eigen::Matrix<T, N, 1>& v1) :
                        m_dt_max(dt), m_p0(p0), m_v0(v0);
{
    // Cast constants to proper type
  constexpr T t0p5 {static_cast<T>(0.5)};
  constexpr T t002 {static_cast<T>(2)};
  constexpr T t006 {static_cast<T>(6)};

    // Computed polynomial coefficients
  auto invdt = static_cast<T>(1)/dt;
  m_j0 = t006*invdt*invdt*(v0 + v1 - 2.0*invdt*(p1 - p0));
  m_a0 = invdt*(v1 - v0) - t0p5*m_j0*dt;

  m_dt_min -= dt_eps;
  m_dt_max += dt_eps;
}
  

template<typename T, int N>
Eigen::Matrix<T, N, 1> Hermite1<T,N>::getPosition(T dt) const
{
  if (dt < m_dt_min  ||  dt > m_dt_max) {
    throw std::invalid_argument("Hermite1<T,N>::getPosition(T dt) - bad dt");
  }

  constexpr T tf2 {static_cast<T>(1.0)/(static_cast<T>(2.0))};
  constexpr T tf3 {static_cast<T>(1.0)/(static_cast<T>(3.0))};

  return m_p0 +     dt*(m_v0 +
                tf2*dt*(m_a0 +
                tf3*dt*(m_j0)));
}
  

template<typename T, int N>
Eigen::Matrix<T, N, 1> Hermite1<T,N>::getVelocity(T dt) const
{
  if (dt < m_dt_min  ||  dt > m_dt_max) {
    throw std::invalid_argument("Hermite1<T,N>::getVelocity(T dt) - bad dt");
  }

  constexpr T tf2 {static_cast<T>(1.0)/(static_cast<T>(2.0))};

  return m_v0 +     dt*(m_a0 +
                tf2*dt*(m_j0));
}
  

template<typename T, int N>
Eigen::Matrix<T, N, 1> Hermite1<T,N>::getAcceleration(T dt) const
{
  if (dt < m_dt_min  ||  dt > m_dt_max) {
    throw std::invalid_argument("Hermite1<T,N>::getAcceleration() - bad dt");
  }

  return m_a0 + dt*m_j0;
}


}

#endif
