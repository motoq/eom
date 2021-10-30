/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_ecfeci_sys.h>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <sofa.h>

#include <cal_julian_date.h>
#include <cal_greg_date.h>
#include <cal_duration.h>

namespace eom {

/*
 * Nomenclature
 *   ITRF  International terrestrial reference frame (ECEF)
 *   TIRF  Terrestrial intermediate reference frame (PEF)
 *   CIRF  Celestial intermediate reference frame (non-equinox based TETE)
 *   GCRF  Geocentric  celestial reference frame (ECI)
 */

EcfEciSys::EcfEciSys(JulianDate startTime, JulianDate stopTime, Duration dt) :
                     jdStart {startTime}, jdStop {stopTime},
                                          dt_days {dt.getDays()}
{
  auto jd = jdStart;
    // Convert to TT
  double x;                       // Celestial Intermediate Pole x coordinate
  double y;                       // Celestial Intermediate Pole y coordinate
  double s;                       // CIO locator, radians
  double cirf2gcrf[3][3];              // intermediate-to-Celestial matrix
  double itrf2tirf[3][3];         // ITRF to TIRF
  while (jd < jdStop) {
    auto jdTT = jd;
      // Pole locations for BPN
    iauXys06a(jdTT.getJdHigh(), jdTT.getJdLow(), &x, &y, &s);
      // Transpose of BPN, to BPN
    iauC2ixys(x, y, s, cirf2gcrf); 
    iauTr(cirf2gcrf, cirf2gcrf);
      // CIO locator used for polar motion transformation
    double sp {iauSp00(jdTT.getJdHigh(), jdTT.getJdLow())};
      // Polar motion, ECF to ECI after transpose
    double xp {0.0};
    double yp {0.0};
    iauPom00(xp, yp, sp, itrf2tirf);
    iauTr(itrf2tirf, itrf2tirf);

    Eigen::Matrix3d mtmp;
    Eigen::Quaterniond q(mtmp);
    ecf_eci f2i {jd.getJd2000(), 0.0, 0.0, q, q};
    f2iData.push_back(f2i);

    jd += dt_days;
  }
    
}

}
