/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_GRANULE_H
#define ASTRO_GRANULE_H

#include <array>
#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>
#include <mth_chebyshev.m>

namespace eom {

/**
 * Creates an ephemeris "granule" composed of 
 *
 * Reference:  X X Newhall, "Numerical Representation of Planetary
 *             Ephemerides", Jet Propulsion Laboratory, 1989.
 *
 * @tparam  ORDER  Order of the polynomial (highest exponent)
 * @tparam  N      Number of fit points
 *
 * @author  Kurt Motekew
 * @date    2023/02/22
 */
template<int ORDER, int N>
class Granule {
public:
  /**
   * Initialize with a set of position and velocity vectors
   *
   * @param  ts  Times associated with each position and velocity
   *             vector
   * @param  ps  3xN matrix of position vectors
   * @param  vs  3xN matrix of velocity vectors
   */
  Granule(const std::array<JulianDate, N>& ts,
          const Eigen::Matrix<double, 3, N>& ps,
          const Eigen::Matrix<double, 3, N>& vs);

  /**
   * @return  Earliest time for which state can be retrieved
   */
  JulianDate getBeginTime() const
  {
    return m_jdStart;
  }

  /**
   * @return  Latest time for which state can be retrieved
   */
  JulianDate getEndTime() const
  {
    return m_jdStop;
  }

  /**
   * Return interpolated position
   *
   * @param  jd  Time for which to retrieve position from this granule
   *
   * @return  The interpolated position, DU
   *
   * @throws  invalid_argument if the requested time is out of range
   */
  Eigen::Matrix<double, 3, 1> getPosition(const JulianDate& jd) const;


  /**
   * Return interpolated velocity
   *
   * @param  jd  Time for which to retrieve velocity from this granule
   *
   * @return  The interpolated velocity, DU/TU
   *
   * @throws  invalid_argument if the requested time is out of range
   */
  Eigen::Matrix<double, 3, 1> getVelocity(JulianDate jd) const;

private:
  JulianDate m_jdStart;
  JulianDate m_jdStop;
  double m_days {};
  double m_dt_norm {1.0};
  double m_dt_shift {};
  Eigen::Matrix<double, 3, 1> m_ax;
  Eigen::Matrix<double, 3, 1> m_ay;
  Eigen::Matrix<double, 3, 1> m_az;
};


template<int ORDER, int N>
Granule<ORDER,N>::Granule(const std::array<JulianDate, N>& ts,
                          const Eigen::Matrix<double, 3, N>& ps,
                          const Eigen::Matrix<double, 3, N>& vs)
{
  // N > O assert

  m_jdStart = ts[0];
  m_jdStop = ts[N-1];
  
  m_days = m_jdStop - m_jdStart;
  m_dt_norm = 0.5*phy_const::tu_per_day*m_days;
  m_dt_shift = 0.0 + m_dt_norm;

    // N is the number of position or velocity observations
    // 2N is the total number of observations - alternate pos then vel
  Eigen::Matrix<double, 2*N, ORDER+1> tmat;
  Eigen::Matrix<double, 2*N, 1> fx;
  Eigen::Matrix<double, 2*N, 1> fy;
  Eigen::Matrix<double, 2*N, 1> fz;

  for (int ii=0; ii<N; ++ii) {
    double tu {phy_const::tu_per_day*(ts[ii] - m_jdStart)};
    double dt {(tu - m_dt_shift)/m_dt_norm};
    int prow {2*ii};
    int vrow {prow+1};
    tmat.block(prow,0,0,ORDER) = chebyshev::poly<double, ORDER>(dt);
    tmat.block(vrow,0,0,ORDER) = chebyshev::poly_dot<double, ORDER>(dt);
    fx(prow) = ps(0, ii);
    fx(vrow) = vs(0, ii);
    fy(prow) = ps(1, ii);
    fy(vrow) = vs(1, ii);
    fz(prow) = ps(2, ii);
    fz(vrow) = vs(2, ii);
  }
  //ColPivHouseholderQR<Matrix3f> dec(A);
  //Vector3f x = dec.solve(b);
  m_ax = tmat.colPivHouseholderQr().solve(fx);
  m_ay = tmat.colPivHouseholderQr().solve(fy);
  m_az = tmat.colPivHouseholderQr().solve(fz);
}
  

template<int ORDER, int N>
Eigen::Matrix<double, 3, 1>
Granule<ORDER,N>::getPosition(const JulianDate& jd) const
{
  double tu {phy_const::tu_per_day*(jd - m_jdStart)};
  if (tu < -phy_const::epsdt  ||  tu > phy_const::epsdt) {
    throw std::invalid_argument("Granule<T,N>::getPosition() - bad jd");
  }
  double dt {(tu - m_dt_shift)/m_dt_norm};
  Eigen::Matrix<double, 1, ORDER+1> tpoly = chebyshev::poly<double, ORDER>(dt);
  Eigen::Matrix<double, 3, 1> pos = { m_ax.transpose()*tpoly,
                                      m_ay.transpose()*tpoly,
                                      m_az.transpose()*tpoly };
  return pos;
}
  

template<int ORDER, int N>
Eigen::Matrix<double, 3, 1>
Granule<ORDER,N>::getVelocity(const JulianDate& jd) const
{
  double tu {phy_const::tu_per_day*(jd - m_jdStart)};
  if (tu < -phy_const::epsdt  ||  tu > phy_const::epsdt) {
    throw std::invalid_argument("Granule<T,N>::getVelocity() - bad jd");
  }
  double dt {(tu - m_dt_shift)/m_dt_norm};
  Eigen::Matrix<double, 1, ORDER+1> tpoly = chebyshev::poly_dot<double,
                                                                ORDER>(dt);
  Eigen::Matrix<double, 3, 1> vel = { m_ax.transpose()*tpoly,
                                      m_ay.transpose()*tpoly,
                                      m_az.transpose()*tpoly };
  return vel;
}


}

#endif
