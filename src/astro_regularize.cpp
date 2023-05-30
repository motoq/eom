#include <astro_regularize.h>

#include <cmath>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_duration.h>
#include <astro_keplerian.h>

namespace {
  constexpr int nint {18};
  constexpr double dtheta {utl_const::pi/nint};
    // Regularization order and integration steps per revolution
  constexpr double ord {1.5};
  constexpr double sdiv {120.0};
}

namespace eom {

/*
 * Consult Berry & Healy, "The generalized Sundman transformation for
 * propagation of high-eccentricity elliptical orbits", for background
 * info.
 *
 * Still a little bit experimental...
 */
Regularize::Regularize(const Eigen::Matrix<double, 6, 1>& x,
                       const Eigen::Matrix<double, 6, 1>& dx)
{
    // Determine maximum step size - get orbit parameters
  Keplerian oe(x);
  double ecc {oe.getEccentricity()};
    // Integrate over intermediate anomaly with order 1.5
  double theta  {0.0};
    // (1 + e)^(1.5 - 2) + (1 - e)^(1.5 - 2)
  double s2pi {0.5*(1.0/std::sqrt(1.0 + ecc) + 1.0/std::sqrt(1.0 - ecc))};
  for (int ii=0; ii<(nint-1); ++ii) {
    theta += dtheta;
    s2pi += 1.0/std::sqrt(1.0 + ecc*std::cos(theta));
  }
  s2pi *= 2.0*dtheta;
  m_ds = s2pi/sdiv;

  Duration time(0.0, 0.0);
  setTimeState(time, x, dx);
}


void Regularize::setTimeState(const Duration& time,
                              const Eigen::Matrix<double, 6, 1>& x,
                              const Eigen::Matrix<double, 6, 1>& dx)
{
  m_time = time;
  m_x = x;
  m_dx = dx;

  Eigen::Matrix<double, 3, 1> rvec = m_x.block<3,1>(0,0);
  Eigen::Matrix<double, 3, 1> vvec = m_x.block<3,1>(3,0);
  double r2 {rvec.dot(rvec)};
  double rmag {std::sqrt(r2)};
  double rdotv {rvec.dot(vvec)};
  double orvu {ord*rmag*rdotv/phy_const::gm};
  double dxsdyt2inv {rmag*r2/phy_const::gm};
  double dxsdytinv {std::sqrt(dxsdyt2inv)};

  m_y(0) = m_time.getTu();
  m_y(4) = dxsdytinv;
  m_dy(0) = m_y(4);
  m_dy(4) = orvu;
  m_y.block<3,1>(1,0) = m_x.block<3,1>(0,0);
  m_y.block<3,1>(5,0) = dxsdytinv*m_dx.block<3,1>(0,0);
  m_dy.block<3,1>(1,0) = m_y.block<3,1>(5,0);
  m_dy.block<3,1>(5,0) = orvu*m_dx.block<3,1>(0,0) +
                         dxsdyt2inv*m_dx.block<3,1>(3,0);
}


void Regularize::setRegularizedState(const Eigen::Matrix<double, 8, 1>& y,
                                     const Eigen::Matrix<double, 8, 1>& dy)
{
  m_y = y;
  m_dy = dy;
  m_time.set(m_y(0), 1.0);

  double r2 {m_y(1)*m_y(1) + m_y(2)*m_y(2) + m_y(3)*m_y(3)};
  double rmag {std::sqrt(r2)};
  double dxsdyt2 {phy_const::gm/(rmag*r2)};
  double dxsdyt {std::sqrt(dxsdyt2)};
  Eigen::Matrix<double, 3, 1> rvec = m_x.block<3,1>(0,0);
  Eigen::Matrix<double, 3, 1> vvec = m_x.block<3,1>(3,0);
  double rdotv {rvec.dot(vvec)};

  m_x.block<3,1>(0,0) = m_y.block<3,1>(1,0);
  m_dx.block<3,1>(0,0) = dxsdyt*m_dy.block<3,1>(1,0);
  m_x.block<3,1>(3,0) = m_dx.block<3,1>(0,0);
  m_dx.block<3,1>(3,0) = dxsdyt2*m_dy.block<3,1>(5,0) -
                         (ord*rdotv/r2)*m_dx.block<3,1>(0,0);
}


}
