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

namespace eom {

EcfEciSys::EcfEciSys(JulianDate startTime, JulianDate stopTime, Duration dt) :
                     jdStart {startTime}, jdStop {stopTime},
                                          dt_days {dt.getDays()}
{
    // Put a small buffer around start and stop to minimize logic
    // locating ECFECI data
  double offset {0.25*dt_days};
  jdStart += -offset;
  jdStop  +=  offset;

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
    ecf_eci f2i {jd.getJd2000(), 0.0, 0.0, qpm, qbpn};
    f2iData.push_back(f2i);

    jd += dt_days;
  }
}


ecf_eci EcfEciSys::getEcfEciData(JulianDate& utc)
{
  double days {utc - jdStart};
  if (days < 0.0) {
    throw std::out_of_range ("EcfEciSys::getEcfEciData Time out of range");
  }

  unsigned long int ndx0 {static_cast<unsigned long int>(days/dt_days)};
  unsigned long int ndx1 {ndx0 + 1UL};

  ecf_eci& f2i0 = f2iData[ndx0];
  ecf_eci& f2i1 = f2iData[ndx1];

  ecf_eci f2i;
  f2i.jd2000 = (f2i0.jd2000 + f2i1.jd2000)/2.0;
  return f2i;
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
