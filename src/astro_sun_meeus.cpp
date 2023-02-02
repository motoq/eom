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
}


Eigen::Matrix<double, 6, 1> SunMeeus::getStateVector(const JulianDate& jd,
                                                     EphemFrame frame) const 
{
  Eigen::Matrix<double, 6, 1> xeci;
  xeci.block<3,1>(0,0) = this->getPosition(jd, EphemFrame::eci);

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
  LeapSeconds& ls = eom::LeapSeconds::getInstance();
  auto jdTT = ls.utc2tt(jd);
  double jdCent {jdTT.getJulianCenturies()};

    // Eccentricity of Earth orbit
  auto ecc = 0.016708634 - jdCent*(0.000042037 + jdCent*0.0000001267);
    // Mean longitude of the Sun w.r.t. the mean equinox of date
    // and mean anomaly (degrees to radians)
  auto el0 = 280.46646 + jdCent*(36000.76983 + jdCent*0.0003032);
  auto em  = 357.52911 + jdCent*(35999.05029 - jdCent*0.0001537);
  el0 *= utl_const::rad_per_deg;
  em *= utl_const::rad_per_deg;
    // Sun equation of center
  auto cee = (1.914602 - jdCent*(0.004817 + jdCent*0.000014))*std::sin(em) +
             (0.019993 - jdCent*0.000101)*std::sin(2.0*em) +
             0.000289*std::sin(3.0*em);
  cee *= utl_const::rad_per_deg;
    // Sun true longitude w.r.t. the mean equinox of the date,
    // true anomaly, and radial distance (AU)
  auto lon_sun = el0 + cee;
  auto nu_sun = em + cee;
  auto r_sun = 1.000001018*(1.0 - ecc*ecc)/(1.0 + ecc*std::cos(nu_sun));
    // Convert to DU using Meeus value
  r_sun *= 149597870.0*phy_const::du_per_km;
    // J2000 equinox
  lon_sun -= utl_const::rad_per_deg*(0.01397*jdTT.getMjd2000()/364.25);

    // Mean obliquity of the ecliptic
  auto e0 = (23.0 + 26.0/60.0 + 21.448/3600.0) - jdCent*(46.8150/3600.0 +
                                                 jdCent*(0.00059/3600.0 -
                                                 jdCent*(0.001813/3600.0)));
  e0 *= utl_const::rad_per_deg;

  auto ra = std::atan2(std::cos(e0)*std::sin(lon_sun), std::cos(lon_sun));
  auto de = std::asin(std::sin(e0)*std::sin(lon_sun));
  Eigen::Matrix<double, 3, 1> xeci;
  xeci(0) = r_sun*std::cos(de)*std::cos(ra);
  xeci(1) = r_sun*std::cos(de)*std::sin(ra);
  xeci(2) = r_sun*std::sin(de);

  if (frame == EphemFrame::ecf) {
    return m_ecfeci->eci2ecf(jd, xeci);
  }
  return xeci;
}


}
