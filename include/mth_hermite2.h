/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_HERMITE2_H
#define MTH_HERMITE2_H

#include <stdexcept>

#include <Eigen/Dense>

namespace eom {

/**
 * Performs Hermite interpolation between two nodes making use of the
 * 1st and 2nd derivatives of the function being approximated.  The
 * nodes can then serve as knots preserving continuity through the
 * second derivative between adjacent Hermite2 polynomials when used
 * for piecewise polynomial interpolation.  N-dimensional vector
 * interpolation is supported.  Function values are referred to as
 * "position" vectors, with the first and second derivatives "velocity"
 * and "acceleration".
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
class Hermite2 {
public:
  /**
   * Initialize with two sets of position, velocity, and acceleration
   * vectors, and the time between them.  Acceleration must be included.
   *
   * @param  dt      Spacing between nodes; e.g., time from p0 to p1
   * @param  p0      Initial state; e.g., position, DU
   * @param  v0      1st derivative of initial state, DU/TU; e.g., velocity
   * @param  a0      2nd derivative of initial state, DU/TU^2;
   *                 e.g., acceleration
   * @param  p1      Final state
   * @param  v1      1st derivative of final state, DU/TU
   * @param  a1      2nd derivative of final state, DU/TU^2
   * @param  dt_eps  Optional endpoint tolerance.  When performing
   *                 interpolation, the input time may preceed 0 or
   *                 exceed dt by this amount, and still be considered
   *                 valid.  This is to eliminate throwing unecessary
   *                 exceptions due to time errors based on roundoff error
   *                 when subtracting/scaling time values passed to the
   *                 interpolator.
   */
  Hermite2(T dt,
           const Eigen::Matrix<T, N, 1>& p0,
           const Eigen::Matrix<T, N, 1>& v0,
           const Eigen::Matrix<T, N, 1>& a0,
           const Eigen::Matrix<T, N, 1>& p1,
           const Eigen::Matrix<T, N, 1>& v1,
           const Eigen::Matrix<T, N, 1>& a1,
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
   *
   * @throws  invalid_argument if the requested time is out of the
   *          polynomial range.
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
   *
   * @throws  invalid_argument if the requested time is out of the
   *          polynomial range.
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
   *
   * @throws  invalid_argument if the requested time is out of the
   *          polynomial range.
   */
  Eigen::Matrix<T, N, 1> getAcceleration(T dt) const;

private:
  T m_dt_min {0};
  T m_dt_max {};
  Eigen::Matrix<T, N, 1> m_p0;                   // position
  Eigen::Matrix<T, N, 1> m_v0;                   // velocity
  Eigen::Matrix<T, N, 1> m_a0;                   // acceleration
  Eigen::Matrix<T, N, 1> m_j0;                   // jerk
  Eigen::Matrix<T, N, 1> m_k0;                   // jitter
  Eigen::Matrix<T, N, 1> m_l0;                   // dither ("Thanks Ed")
};


template<typename T, int N>
Hermite2<T,N>::Hermite2(T dt,
                        const Eigen::Matrix<T, N, 1>& p0,
                        const Eigen::Matrix<T, N, 1>& v0,
                        const Eigen::Matrix<T, N, 1>& a0,
                        const Eigen::Matrix<T, N, 1>& p1,
                        const Eigen::Matrix<T, N, 1>& v1,
                        const Eigen::Matrix<T, N, 1>& a1,
                        T dt_eps) :
                        m_dt_max(dt), m_p0(p0), m_v0(v0), m_a0(a0)
{
    // Cast constants to proper type
  constexpr T t0p5 {static_cast<T>(0.5)};
  constexpr T t002 {static_cast<T>(2)};
  constexpr T t003 {static_cast<T>(3)};
  constexpr T t004 {static_cast<T>(4)};
  constexpr T t006 {static_cast<T>(6)};
  constexpr T t007 {static_cast<T>(7)};
  constexpr T t015 {static_cast<T>(15)};
  constexpr T t060 {static_cast<T>(60)};
    // Temporary values
  auto invdt = static_cast<T>(1)/dt;
  Eigen::Matrix<T, N, 1> cpos = -t006*invdt*(t0p5*a0 +
                                      invdt*(v0 - invdt*(p1 - p0)));
  Eigen::Matrix<T, N, 1> cvel = -t002*invdt*(a0 - invdt*(v1 - v0));
  Eigen::Matrix<T, N, 1> cacc =       invdt*(a1 - a0);

    // Computed polynomial coefficients
  m_l0 = t060*(t002*cpos - t003*cvel + cacc)*invdt*invdt;
  m_k0 = t004*(cacc - cpos)*invdt - t007*m_l0*dt/t015;
  m_j0 = cacc - t0p5*dt*(m_k0 + m_l0*dt/t003);

  m_dt_min -= dt_eps;
  m_dt_max += dt_eps;
}
  

template<typename T, int N>
Eigen::Matrix<T, N, 1> Hermite2<T,N>::getPosition(T dt) const
{
  if (dt < m_dt_min  ||  dt > m_dt_max) {
    throw std::invalid_argument("Hermite2<T,N>::getPosition(T dt) - bad dt");
  }

  constexpr T tf2 {static_cast<T>(1.0)/(static_cast<T>(2.0))};
  constexpr T tf3 {static_cast<T>(1.0)/(static_cast<T>(3.0))};
  constexpr T tf4 {static_cast<T>(1.0)/(static_cast<T>(4.0))};
  constexpr T tf5 {static_cast<T>(1.0)/(static_cast<T>(5.0))};

  return m_p0 +     dt*(m_v0 +
                tf2*dt*(m_a0 +
                tf3*dt*(m_j0 +
                tf4*dt*(m_k0 +
                tf5*dt*(m_l0)))));
}
  

template<typename T, int N>
Eigen::Matrix<T, N, 1> Hermite2<T,N>::getVelocity(T dt) const
{
  if (dt < m_dt_min  ||  dt > m_dt_max) {
    throw std::invalid_argument("Hermite2<T,N>::getVelocity(T dt) - bad dt");
  }

  constexpr T tf2 {static_cast<T>(1.0)/(static_cast<T>(2.0))};
  constexpr T tf3 {static_cast<T>(1.0)/(static_cast<T>(3.0))};
  constexpr T tf4 {static_cast<T>(1.0)/(static_cast<T>(4.0))};

  return m_v0 +     dt*(m_a0 +
                tf2*dt*(m_j0 +
                tf3*dt*(m_k0 +
                tf4*dt*(m_l0))));
}
  

template<typename T, int N>
Eigen::Matrix<T, N, 1> Hermite2<T,N>::getAcceleration(T dt) const
{
  if (dt < m_dt_min  ||  dt > m_dt_max) {
    throw std::invalid_argument("Hermite2<T,N>::getAcceleration() - bad dt");
  }

  constexpr T tf2 {static_cast<T>(1.0)/(static_cast<T>(2.0))};
  constexpr T tf3 {static_cast<T>(1.0)/(static_cast<T>(3.0))};

  return m_a0 +     dt*(m_j0 +
                tf2*dt*(m_k0 +
                tf3*dt*(m_l0)));
}


}

#endif
