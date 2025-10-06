#include <astro_secular_j2.h>

#include <array>
#include <cmath>
#include <memory>
#include <utility>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_keplerian.h>

namespace {
  constexpr int io {3};           // RAAN
  constexpr int iw {4};           // Argument of perigee
}

namespace eom {

SecularJ2::SecularJ2(const std::string& orbit_name,
                     const JulianDate& epoch,
                     const Eigen::Matrix<double, 6, 1>& xeci,
                     std::shared_ptr<const EcfEciSys> ecfeciSys)
{
  m_name = orbit_name;
  m_jd0 = epoch;
  m_ecfeci = std::move(ecfeciSys);
    // A true equator ECI frame is required for propagation
  Eigen::Matrix<double, 6, 1> xecf = m_ecfeci->eci2ecf(m_jd0,
                                                       xeci.block<3,1>(0,0),
                                                       xeci.block<3,1>(3,0));
  Eigen::Matrix<double, 6, 1> xteme = m_ecfeci->ecf2teme(m_jd0,
                                                         xecf.block<3,1>(0,0),
                                                         xecf.block<3,1>(3,0));
  Keplerian kep {xteme};
    // Epoch state
  m_oe0 = kep.getOrbitalElements();
  m_m0 = kep.getMeanAnomaly();

  const double ecc {kep.getEccentricity()};
  const double slr {kep.getSemilatusRectum()};
  const double p2 {slr*slr};
  const double re2j2_p2 {phy_const::re2j2/p2};
  const double inc {kep.getInclination()};
  const double si {std::sin(inc)};
  const double si2 {si*si};

  const double m0dot_n {1.5*re2j2_p2*std::sqrt(1.0 - ecc*ecc)*(1.0 - 1.5*si2)};
  m_nbar = (1.0 + m0dot_n)*kep.getMeanMotion();

  const double xx_x {m_nbar*re2j2_p2};
  m_odot = -1.5*xx_x*std::cos(inc);
  m_wdot =  1.5*xx_x*(2.0 - 2.5*si2);
}


Eigen::Matrix<double, 6, 1> SecularJ2::getStateVector(const JulianDate& jd,
                                                      EphemFrame frame) const
{
  double dt {phy_const::tu_per_day*(jd - m_jd0)};
  std::array<double, 6> oe = m_oe0;
    // Update slow moving OE subject to 1st order J2 secular effects
  oe[io] += dt*m_odot;
  oe[iw] += dt*m_wdot;
  oe[io] = std::fmod(oe[io], utl_const::tpi);
  oe[iw] = std::fmod(oe[iw], utl_const::tpi);

  Keplerian kep {oe};
    // Update fast parameter
  kep.setWithMeanAnomaly(m_m0 + m_nbar*dt);
    // Cartesian - TEME to earth fixed
  Eigen::Matrix<double, 6, 1> xteme {kep.getCartesian()};
  Eigen::Matrix<double, 6, 1> xecf = m_ecfeci->teme2ecf(jd,
                                                        xteme.block<3,1>(0,0),
                                                        xteme.block<3,1>(3,0));
  if (frame == EphemFrame::eci) {
    return m_ecfeci->ecf2eci(jd, xecf.block<3,1>(0,0), xecf.block<3,1>(3,0));
  }
  return xecf;
}


Eigen::Matrix<double, 3, 1> SecularJ2::getPosition(const JulianDate& jd,
                                                   EphemFrame frame) const
{
  Eigen::Matrix<double, 6, 1> xdx_ecf = getStateVector(jd,
                                                       EphemFrame::ecf);
  Eigen::Matrix<double, 3, 1> xecf = xdx_ecf.block<3,1>(0,0);

  if (frame == EphemFrame::eci) {
    return m_ecfeci->ecf2eci(jd, xecf);
  }
  return xecf;
}


}
