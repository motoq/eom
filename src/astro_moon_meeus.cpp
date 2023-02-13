/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_moon_meeus.h>

#include <memory>
#include <cmath>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_julian_date.h>
#include <cal_leap_seconds.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>

#include "astro_meeus_t47.h"

namespace {

/**
 * Reduce the input angle to less than 360 deg
 *
 * @param  x  Input angle, degrees
 *
 * @return  Output angle reduced to less than 360 degrees
 */
constexpr double fnred(double x)
{
  return x - 360.0*static_cast<int>(x/360.0);
}

}

namespace eom {


MoonMeeus::MoonMeeus(const std::shared_ptr<const EcfEciSys>& ecfeciSys,
                     const std::string& name)
{
  m_ecfeci = ecfeciSys;
  m_name = name;
  m_jdStart = m_ecfeci->getBeginTime();
  m_jdStop = m_ecfeci->getEndTime();
}


Eigen::Matrix<double, 6, 1> MoonMeeus::getStateVector(const JulianDate& jd,
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


Eigen::Matrix<double, 3, 1> MoonMeeus::getPosition(const JulianDate& jd,
                                                   EphemFrame frame) const
{
  using namespace meeus_t47;

  if (jd < m_jdStart  ||  m_jdStop < jd) {
    throw std::out_of_range("MoonMeeus::getPosition() - bad time");
  }

  auto jdTT = eom::LeapSeconds::getInstance().utc2tt(jd);
  double jdCent {jdTT.getJulianCenturies()};

    // Moon's mean longitude, w.r.t. mean equinox of date
  auto el_prime = 218.3164477 + jdCent*(481267.88123421 -
                                jdCent*(0.0015786 -
                                jdCent*(1.0/538841.0 -
                                jdCent*(1.0/65194000.0))));
  el_prime = fnred(el_prime);
    // Mean elongation of the Moon
  auto dee = 297.8501921 + jdCent*(445267.1114034 -
                           jdCent*(0.0018819 -
                           jdCent*(1.0/545868.0 -
                           jdCent*(1.0/113065000.0))));
  dee = fnred(dee);
    // Sun's mean anomaly
  auto em  = 357.5291092 + jdCent*(35999.0502909 -
                           jdCent*(0.0001536 -
                           jdCent*(1.0/24490000.0)));
  em = fnred(em);
    // Moon's mean anomaly
  auto em_prime = 134.9633964 + jdCent*(477198.8675055 +
                                jdCent*(0.0087414 +
                                jdCent*(1.0/69699.0 -
                                jdCent*(1.0/14712000.0))));
  em_prime = fnred(em_prime);
    // Moon's argument of latitude
  auto eff = 93.2720950 + jdCent*(483202.0175233 -
                          jdCent*(0.0036539 +
                          jdCent*(1.0/3526000.0 -
                          jdCent*(1.0/863310000.0))));
  eff = fnred(eff);
    // Three more arguments
  auto a1 = 119.75 + jdCent*131.849;
  a1 = fnred(a1);
  auto a2 = 53.09 + jdCent*479264.290;
  a2 = fnred(a2);
  auto a3 = 313.45 + jdCent*481266.484;
  a3 = fnred(a3);

  auto ecc = 1.0 - jdCent*(0.002516 + jdCent*0.0000074);

  double sum_lon {0.0};
  double sum_rng {0.0};
  double sum_lat {0.0};
  for (int ii=(nt47-1); ii>=0; --ii) {
    auto aoff = ii*atcols;
    auto boff = ii*btcols;

      // Periodic term for longitude and distance
    double lonrt {terms47a[aoff]*dee +
                  terms47a[aoff+1]*em +
                  terms47a[aoff+2]*em_prime +
                  terms47a[aoff+3]*eff};
    auto lon = std::sin(lonrt)*terms47a[aoff+4];
    auto rng = std::cos(lonrt)*terms47a[aoff+5];
      // Periodic term for latitude
    double latt {terms47b[boff]*dee +
                 terms47b[boff+1]*em +
                 terms47b[boff+2]*em_prime +
                 terms47b[boff+3]*eff};
    auto lat = std::sin(latt)*terms47b[boff+4];
      // Adjust for earth orbit eccentricity
    for (int jj=0; jj<std::abs(terms47a[aoff+1]); ++jj) {
      lon *= ecc;
      rng *= ecc;
    }
    for (int jj=0; jj<std::abs(terms47b[boff+1]); ++jj) {
      lat *= ecc;
    }
    sum_lon += lon;
    sum_rng += rng;
    sum_lat += lat;
  }
  auto a1_rad = utl_const::rad_per_deg*a1;
  auto a2_rad = utl_const::rad_per_deg*a2;
  auto a3_rad = utl_const::rad_per_deg*a3;
  auto el_prime_rad = utl_const::rad_per_deg*el_prime;
  auto em_prime_rad = utl_const::rad_per_deg*em_prime;
  auto eff_rad = utl_const::rad_per_deg*eff;
  sum_lon += 3958.0*std::sin(a1_rad) +
             1962.0*std::sin(el_prime_rad - eff_rad) +
              318.0*std::sin(a2_rad);
  sum_lat += -2235.0*std::sin(el_prime_rad) +
               382.0*std::sin(a3_rad) +
               175.0*std::sin(a1_rad - eff_rad) +
               175.0*std::sin(a1_rad + eff_rad) +
               127.0*std::sin(el_prime_rad - em_prime_rad) -
               115.0*std::sin(el_prime_rad + em_prime_rad);

  auto lon_moon = el_prime + sum_lon/1000000.0;
  auto lat_moon = sum_lat/1000000.0;
  auto rng_moon = 385000.56 + sum_rng/1000;

    // Mean obliquity of the ecliptic, seconds
  auto e0 = 21.448 + 60.0*(26.0 + 60.0*23) - jdCent*(46.8150 +
                                             jdCent*(0.00059 -
                                             jdCent*(0.001813)));
    // degrees
  e0 /= 3600.0;

    // Right ascension (J2000) and declination to Cartesian
  auto e0_rad = utl_const::rad_per_deg*e0;
  auto se0 = std::sin(e0_rad);
  auto ce0 = std::cos(e0_rad);
    //
  auto lon_moon_rad = utl_const::rad_per_deg*lon_moon;
  auto slon = std::sin(lon_moon_rad);
  auto clon = std::cos(lon_moon_rad);
    //
  auto rng_moon_du = phy_const::du_per_km*rng_moon;
    //
  auto ra = std::atan2(ce0*slon, clon);
  auto sde = se0*slon;
  auto de = std::asin(sde);
  auto cde = cos(de);
  Eigen::Matrix<double, 3, 1> xeci {rng_moon_du*cde*std::cos(ra),
                                    rng_moon_du*cde*std::sin(ra),
                                    rng_moon_du*sde};
  xeci = m_ecfeci->mod2eci(jd, xeci);

  if (frame == EphemFrame::ecf) {
    return m_ecfeci->eci2ecf(jd, xeci);
  }
  return xeci;
}


}
