/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_gpx.h>

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

/*
 * Currently J2 secular with 2nd order components.
 */
GpX::GpX(const std::string& orbit_name,
         const JulianDate& epoch,
         const Eigen::Matrix<double, 6, 1>& xeci,
         std::shared_ptr<const EcfEciSys> ecfeciSys) :
         m_name {orbit_name},
         m_jd0 {epoch},
         m_ecfeci {std::move(ecfeciSys)}
{
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

    // Compute effective mean motion
  const double ecc {kep.getEccentricity()};
  const double ecc2 {ecc*ecc};
  const double ecc4 {ecc2*ecc2};
  const double slr {kep.getSemilatusRectum()};
  const double p2 {slr*slr};
  const double re2j2_p2 {phy_const::re2j2/p2};
  const double inc {kep.getInclination()};
  const double si {std::sin(inc)};
  const double si2 {si*si};

  m_n = kep.getMeanMotion();
  m_mdot = 1.5*m_n*re2j2_p2*std::sqrt(1.0 - ecc2)*(1.0 - 1.5*si2) +
           3.0*m_n*re2j2_p2*re2j2_p2*(320.0*ecc2 - 280.0*ecc4 +
                                      (1600.0 - 1568.0*ecc2 + 328.0*ecc4)*si2 +
                                      (-2096.0 + 1072*ecc2 + 79.0*ecc4)*si2*si2
                                     )/(512.0*std::sqrt(1.0 - ecc2));
    // RAAN and argument of perigee rates
  const double xx_x {m_n*re2j2_p2};
  m_odot = -1.5*xx_x*std::cos(inc) +
            3.0*m_n*re2j2_p2*re2j2_p2*(12.0 - 4.0*ecc2 -
                                       (80.0 + 5.0*ecc2)*si2
                                      )/32.0;
  m_wdot =  1.5*xx_x*(2.0 - 2.5*si2) +
            9.0*m_n*re2j2_p2*re2j2_p2*(56.0*ecc2 +
                                       (760.0 - 36.0*ecc2)*si2 -
                                       (890.0 + 45.0*ecc2)*si2*si2
                                      )/384.0;
}


std::string GpX::getName() const
{
  return m_name;
}


JulianDate GpX::getEpoch() const
{
  return m_jd0;
}


JulianDate GpX::getBeginTime() const
{
  return m_ecfeci->getBeginTime();
}


JulianDate GpX::getEndTime() const
{
  return m_ecfeci->getEndTime();
}


Eigen::Matrix<double, 6, 1> GpX::getStateVector(const JulianDate& jd,
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
  kep.setWithMeanAnomaly(std::fmod(m_m0 + dt*(m_n + m_mdot), utl_const::tpi));
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


Eigen::Matrix<double, 3, 1> GpX::getPosition(const JulianDate& jd,
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
