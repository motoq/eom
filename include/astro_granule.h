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
#include <mth_chebyshev.h>

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
  Eigen::Matrix<double, 3, 1> getVelocity(const JulianDate& jd) const;

private:
  JulianDate m_jdStart;
  JulianDate m_jdStop;
  double m_days {};
  double m_dt_norm {1.0};
  double m_dt_shift {};
  Eigen::Matrix<double, ORDER+1, 1> m_ax;
  Eigen::Matrix<double, ORDER+1, 1> m_ay;
  Eigen::Matrix<double, ORDER+1, 1> m_az;
  Eigen::Matrix<double, ORDER+1, 1> m_adx;
  Eigen::Matrix<double, ORDER+1, 1> m_ady;
  Eigen::Matrix<double, ORDER+1, 1> m_adz;
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
  Eigen::Matrix<double, N, ORDER+1> tmat;

  for (int ii=0; ii<N; ++ii) {
    double tu {phy_const::tu_per_day*(ts[ii] - m_jdStart)};
    double dt {(tu - m_dt_shift)/m_dt_norm};
    tmat.block(ii,0,1,ORDER+1) = chebyshev::poly<double, ORDER>(dt);
  }
  //ColPivHouseholderQR<Matrix3f> dec(A);
  //Vector3f x = dec.solve(b);
  m_ax = tmat.colPivHouseholderQr().solve(ps.block(0,0,1,N).transpose());
  m_ay = tmat.colPivHouseholderQr().solve(ps.block(1,0,1,N).transpose());
  m_az = tmat.colPivHouseholderQr().solve(ps.block(2,0,1,N).transpose());
  m_adx = tmat.colPivHouseholderQr().solve(vs.block(0,0,1,N).transpose());
  m_ady = tmat.colPivHouseholderQr().solve(vs.block(1,0,1,N).transpose());
  m_adz = tmat.colPivHouseholderQr().solve(vs.block(2,0,1,N).transpose());
}
  

template<int ORDER, int N>
Eigen::Matrix<double, 3, 1>
Granule<ORDER,N>::getPosition(const JulianDate& jd) const
{
  double dtlim {1.0 + phy_const::epsdt/m_dt_norm};
  double tu {phy_const::tu_per_day*(jd - m_jdStart)};
  double dt {(tu - m_dt_shift)/m_dt_norm};
  if (dt < -dtlim  ||  dt > dtlim) {
    throw std::invalid_argument("Granule<T,N>::getPosition() - bad jd");
  }
  Eigen::Matrix<double, 1, ORDER+1> tpoly = chebyshev::poly<double, ORDER>(dt);
  Eigen::Matrix<double, 3, 1> pos = {tpoly*m_ax, tpoly*m_ay, tpoly*m_az};

  return pos;
}
  

template<int ORDER, int N>
Eigen::Matrix<double, 3, 1>
Granule<ORDER,N>::getVelocity(const JulianDate& jd) const
{
  double dtlim {1.0 + phy_const::epsdt/m_dt_norm};
  double tu {phy_const::tu_per_day*(jd - m_jdStart)};
  double dt {(tu - m_dt_shift)/m_dt_norm};
  if (dt < -dtlim  ||  dt > dtlim) {
    throw std::invalid_argument("Granule<T,N>::getVelocity() - bad jd");
  }
  Eigen::Matrix<double, 1, ORDER+1> tpoly = chebyshev::poly<double, ORDER>(dt);
  Eigen::Matrix<double, 3, 1> vel = {tpoly*m_adx, tpoly*m_ady, tpoly*m_adz};

  return vel;
}


}

#endif
