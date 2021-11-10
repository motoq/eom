/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_ecfeci_sys.h>

#include <stdexcept>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <sofa.h>

#include <cal_const.h>
#include <cal_julian_date.h>
#include <cal_greg_date.h>
#include <cal_duration.h>
#include <cal_leap_seconds.h>

/*
 * Local utility for converting a C double[3][3] to an Eigen Matrix3d.
 * Should hopefully be the only need for this utility in the codebase
 * unless other additional libraries, similar to SOFA, rely on C array
 * matrix representations.
 */
static Eigen::Matrix3d from3x3(double mtx[3][3]);

/*
 * Nomenclature (IAU typically uses "System" vs. "Frame", e.g., GCRS
 * vs.. GCRF.
 *   ITRF  International terrestrial reference frame (ECEF)
 *   TIRF  Terrestrial intermediate reference frame (PEF)
 *   CIRF  Celestial intermediate reference frame (non-equinox based TETE)
 *   GCRF  Geocentric  celestial reference frame (ECI)
 */

namespace {
    // Default buffer from jdStart and jdStop with a zero rate_days
  constexpr double offset_0 {0.25};
    // Percentage of rate_days to use as buffer for jdStart and jdStop
  constexpr double p_offset {0.25};
}

namespace eom {

EcfEciSys::EcfEciSys(const JulianDate& startTime, const JulianDate& stopTime,
                     const Duration& dt, bool interpolate) :
                     jdStart {startTime}, jdStop {stopTime},
                     rate_days {dt.getDays()}, interpolate_bpnpm {interpolate}
{
    // Put a small buffer around start and stop to minimize logic
    // locating ECFECI data
  if (rate_days == 0.0) {
    double offset {offset_0};
    jdStart += -offset;
    jdStop  +=  offset;
  } else {
    double offset {p_offset*rate_days};
    jdStart += -offset;
    jdStop  +=  offset;
  }

    // Replace with EOP data source when parsing implemented
  double xp {0.0};
  double yp {0.0};
    // Variables used for intermediate calculations within loop.
  double x;                       // Celestial Intermediate Pole x coordinate
  double y;                       // Celestial Intermediate Pole y coordinate
  double s;                       // CIO locator, radians
  double cirf2gcrf[3][3];         // CIRF to GCRF
  double itrf2tirf[3][3];         // ITRF to TIRF
    //
  eom::LeapSeconds& ls = eom::LeapSeconds::getInstance();
  auto jd = jdStart;
    // If only one message, compute for middle time and set rate
    // to be the full duration of the time period
  if (rate_days == 0.0) {
    rate_days = jdStop - jdStart;
    jd += rate_days/2.0;
  }
  while (jd <= jdStop) {
    auto jdTT = ls.utc2tt(jd);
      // Pole locations for BPN
    iauXys06a(jdTT.getJdHigh(), jdTT.getJdLow(), &x, &y, &s);
      // Transpose of BPN, to BPN
    iauC2ixys(x, y, s, cirf2gcrf); 
    iauTr(cirf2gcrf, cirf2gcrf);
      // CIO locator used for polar motion transformation
    double sp {iauSp00(jdTT.getJdHigh(), jdTT.getJdLow())};
      // Polar motion, ECF to ECI after transpose
    iauPom00(xp, yp, sp, itrf2tirf);
    iauTr(itrf2tirf, itrf2tirf);
      // Convert to final form and insert
    Eigen::Matrix3d mbpn = from3x3(cirf2gcrf);
    Eigen::Matrix3d mpm = from3x3(itrf2tirf);
    Eigen::Quaterniond qbpn(mbpn);
    Eigen::Quaterniond qpm(mpm);
    ecf_eci f2i {jd.getMjd2000(), 0.0, 0.0, qpm, qbpn};
    f2iData.push_back(f2i);

    jd += rate_days;
  }
    // Size fixed after initialization
  nfi = static_cast<unsigned long>(f2iData.size());
}


ecf_eci EcfEciSys::getEcfEciData(const JulianDate& utc) const
{
    // Check for valid date
  double days {utc - jdStart};
  if (days < 0.0  ||  jdStop - utc < 0.0) {
    throw std::out_of_range ("EcfEciSys::getEcfEciData Time out of range");
  }

    // Always need first index
  unsigned long int ndx1 {static_cast<unsigned long int>(days/rate_days)};
  const ecf_eci& f2i1 = f2iData[ndx1];

    // Get second data set if interpolating - otherwise,
    // return data less than or equal to requested time
  if (nfi == 1UL) {
    return f2iData[0UL];
  } else if (interpolate_bpnpm) {
    double mjd2000 {utc.getMjd2000()};
    double dt_days {mjd2000 - f2i1.mjd2000};
    double dt {dt_days/rate_days};
    unsigned long int ndx2 {ndx1 + 1UL};
    const ecf_eci& f2i2 = f2iData[ndx2];
    Eigen::Quaterniond bpn {f2i1.bpn.slerp(dt, f2i2.bpn)};
    Eigen::Quaterniond pm {f2i1.pm.slerp(dt, f2i2.pm)};
    ecf_eci f2i {mjd2000, f2i1.ut1mutc, f2i1.lod, pm, bpn};
    return f2i;
  } else {
    return f2i1;
  }
}


Eigen::Matrix<double, 3, 1>
EcfEciSys::ecf2eci(const JulianDate& utc,
                   const Eigen::Matrix<double, 3, 1>& posf) const
{

  ecf_eci f2i {this->getEcfEciData(utc)};
  auto ut1 {utc + cal_const::day_per_sec*f2i.ut1mutc};
  double era {iauEra00(ut1.getJdHigh(), ut1.getJdLow())};
  Eigen::Quaterniond qera{Eigen::AngleAxisd(era, Eigen::Vector3d::UnitZ())};
  Eigen::Matrix<double, 3, 1> posi = f2i.bpn*qera*f2i.pm*posf;
  return posi;
}


Eigen::Matrix<double, 3, 1>
EcfEciSys::eci2ecf(const JulianDate& utc,
                  const Eigen::Matrix<double, 3, 1>& posi) const
{
  ecf_eci f2i {this->getEcfEciData(utc)};
  auto ut1 {utc + cal_const::day_per_sec*f2i.ut1mutc};
  double era {iauEra00(ut1.getJdHigh(), ut1.getJdLow())};
  Eigen::Quaterniond qera{Eigen::AngleAxisd(-era, Eigen::Vector3d::UnitZ())};
  Eigen::Matrix<double, 3, 1> posf =
                              f2i.pm.conjugate()*qera*f2i.bpn.conjugate()*posi;

  return posf;
}


}

static Eigen::Matrix3d from3x3(double mtx[3][3])
{
  Eigen::Matrix3d emtx;
  emtx(0,0) = mtx[0][0];
  emtx(1,0) = mtx[1][0];
  emtx(2,0) = mtx[2][0];
  emtx(0,1) = mtx[0][1];
  emtx(1,1) = mtx[1][1];
  emtx(2,1) = mtx[2][1];
  emtx(0,2) = mtx[0][2];
  emtx(1,2) = mtx[1][2];
  emtx(2,2) = mtx[2][2];

  return emtx;
}
