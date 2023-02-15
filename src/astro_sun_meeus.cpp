/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_sun_meeus.h>

#include <memory>
#include <cmath>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_julian_date.h>
#include <cal_leap_seconds.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>

namespace eom {


SunMeeus::SunMeeus(const std::shared_ptr<const EcfEciSys>& ecfeciSys,
                   const std::string& name)
{
  m_ecfeci = ecfeciSys;
  m_name = name;
  m_jdStart = m_ecfeci->getBeginTime();
  m_jdStop = m_ecfeci->getEndTime();
}


Eigen::Matrix<double, 6, 1> SunMeeus::getStateVector(const JulianDate& jd,
                                                     EphemFrame frame) const 
{
  Eigen::Matrix<double, 6, 1> xeci;
  xeci.block<3,1>(0,0) = this->getPosition(jd, EphemFrame::eci);

    // One minute delta_t
  xeci.block<3,1>(3,0) =
      (this->getPosition(jd +  1.0*utl_const::day_per_min, EphemFrame::eci) -
       this->getPosition(jd + -1.0*utl_const::day_per_min, EphemFrame::eci))/
       (2.0*phy_const::tu_per_min);

  if (frame == EphemFrame::ecf) {
    return m_ecfeci->eci2ecf(jd, xeci.block<3,1>(0,0), xeci.block<3,1>(3,0));
  }
  return xeci;
}


Eigen::Matrix<double, 3, 1> SunMeeus::getPosition(const JulianDate& jd,
                                                EphemFrame frame) const
{
  if (jd < m_jdStart  ||  m_jdStop < jd) {
    throw std::out_of_range("SunMeeus::getPosition() - bad time");
  }

  auto jdTT = eom::LeapSeconds::getInstance().utc2tt(jd);
  double jdCent {jdTT.getJulianCenturies()};

    // Eccentricity of Earth orbit
  auto ecc = 0.016708634 - jdCent*(0.000042037 + jdCent*0.0000001267);
    // Mean longitude of the Sun w.r.t. the mean equinox of date
    // and mean anomaly (working in degrees)
  auto el0 = 280.46646 + jdCent*(36000.76983 + jdCent*0.0003032);
  auto em  = 357.52911 + jdCent*(35999.05029 - jdCent*0.0001537);
    // Sun equation of center
  auto em_rad = em*utl_const::rad_per_deg;
  auto cee = (1.914602 - jdCent*(0.004817 + jdCent*0.000014))*std::sin(em_rad) +
             (0.019993 - jdCent*0.000101)*std::sin(2.0*em_rad) +
             0.000289*std::sin(3.0*em_rad);

    // Sun true longitude w.r.t. the mean equinox of the date,
  auto lon_sun = el0 + cee;
    // true anomaly
  auto nu_sun = em + cee;
    // radial distance (AU)
  auto r_sun = 1.000001018*(1.0 - ecc*ecc)/
              (1.0 + ecc*std::cos(utl_const::rad_per_deg*nu_sun));
    // Convert to DU using Meeus value
  r_sun *= 149597870.0*phy_const::du_per_km;

    // Mean obliquity of the ecliptic, seconds
  auto e0 = 21.448 + 60.0*(26.0 + 60.0*23) - jdCent*(46.8150 +
                                             jdCent*(0.00059 -
                                             jdCent*(0.001813)));
    // degrees
  e0 /= 3600.0;

    // Right ascension and declination to Cartesian MOD
  auto e0_rad = utl_const::rad_per_deg*e0;
  auto se0 = std::sin(e0_rad);
  auto ce0 = std::cos(e0_rad);
    //
  auto lon_sun_rad = utl_const::rad_per_deg*lon_sun;
  auto slon = std::sin(lon_sun_rad);
  auto clon = std::cos(lon_sun_rad);
    //
  auto ra = std::atan2(ce0*slon, clon);
  auto sde = se0*slon;
  auto de = std::asin(sde);
  auto cde = cos(de);
  Eigen::Matrix<double, 3, 1> xeci {r_sun*cde*std::cos(ra),
                                    r_sun*cde*std::sin(ra),
                                    r_sun*sde};
    // MOD to J2000
  xeci = m_ecfeci->mod2eci(jd, xeci);

  if (frame == EphemFrame::ecf) {
    return m_ecfeci->eci2ecf(jd, xeci);
  }
  return xeci;
}


}
