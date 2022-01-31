/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
//#include <algorithm>
//#include <cmath>
//#include <string>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <astro_earth_surf.h>
#include <astro_earth_xt.h>
#include <astro_ground_point.h>

#include <eom_test.h>

namespace eom_app {

void eom_test_earth_x()
{

  std::cout << "\n\n  === Test:  EarthXt ===";


  eom::GroundPoint posGd (utl_const::pi/6.0, utl_const::pi/3, 1.0);
  eom::GroundPoint locGd (utl_const::pi/3.0, utl_const::pi/5, 0.0);
  Eigen::Matrix<double, 3, 1> pos {posGd.getCartesian()};
  Eigen::Matrix<double, 3, 1> loc {locGd.getCartesian()};
  Eigen::Matrix<double, 3, 1> pnt {loc - pos};
  eom::EarthXt ext {eom::HorizonMode::horizon_never};
  ext.setEarthXt(pos, pnt);
  eom::EarthSurf es {eom::HorizonMode::horizon_never};
  es.setEarthSurf(pos, pnt);

  if (ext.missed()) {
    std::cout << "\nMissed!";
  } else {
    std::cout << "\n  " << (loc - ext.getEarthXt()).norm();
  }

  if (es.missed()) {
    std::cout << "\nMissed!";
  } else {
    std::cout << "\n  " << (loc - es.getEarthSurf()).norm();
  }

  if (!ext.missed() && !es.missed()) {
    std::cout << "\n  " << (ext.getEarthXt() - es.getEarthSurf()).norm();
  }

/*
    // Test grid of lat/lon with random alt
 
    // Track maximum error in geodetic by transforming to Cartesian and
    // back.  Check for convergence issues and max iterations.
  double max_lat_err {0.0};
  double max_lon_err {0.0};
  double max_alt_err {0.0};
  int maxitr {0};
  int non_convergence {0};
  int npts {0};
  int nc1 {0};
  int nc2 {0};
  int nc3a {0};
  int nc3b {0};
  double lat {utl_const::pio2};
  double dlat {-utl_const::pi/36.0};
  while (lat >= -utl_const::pio2) {
    double lon {utl_const::pi};
    double dlon {-utl_const::pi/36.0};
    while (lon >= -utl_const::pi) {
      double alt {7.0};
      double dalt {-0.1};
      while (alt > -0.1) {
        eom::GroundPoint gpLla {lat, lon, alt};
        Eigen::Matrix<double, 3, 1> xyz = gpLla.getCartesian();
        eom::GroundPoint gpXyz {xyz};
        switch (gpXyz.getFukStarter()) {
          case eom::FukStarter::case1:
            nc1++;
            break;
          case eom::FukStarter::case2:
            nc2++;
            break;
          case eom::FukStarter::case3a:
            nc3a++;
            break;
          case eom::FukStarter::case3b:
            nc3b++;
            break;
          default:
            break;
        }
        maxitr = std::max(maxitr, gpXyz.getItr());
        if (gpXyz.getItr() < 0) {
          non_convergence++;
        }
        max_lat_err = std::max(max_lat_err,
                               std::fabs(lat - gpXyz.getLatitude()));
        max_lon_err = std::max(max_lon_err,
                               std::fabs(lon - gpXyz.getLongitude()));
        max_alt_err = std::max(max_alt_err,
                               std::fabs(alt - gpXyz.getAltitude()));
        //std::cout << '\n' << utl_const::deg_per_rad*lat << "  " <<
        //                     utl_const::deg_per_rad*lon << "  " <<
        //                     alt;
        alt += dalt;
        npts++;
      }
      lon += dlon;
    }
    lat += dlat;
  }
  std::cout << "\n\nMax Error Over Grid of " << npts << " points" <<
               "\n  Lat: " << utl_const::deg_per_rad*max_lat_err << " deg" <<
               "   Lon: " << utl_const::deg_per_rad*max_lon_err << " deg" <<
               "   Alt: " << phy_const::m_per_du*max_alt_err << " m";
  std::cout << "\nStarter 1: " << nc1 << "   Starter 2: " << nc2 <<
               "   Starter 3a: " << nc3a << "   Starter 3b: " << nc3b;
  std::cout << "\n  With max iterations " << maxitr <<
               " and " << non_convergence << " non-convergent cases";
*/


  std::cout << "\n  === End Test:  EarthXt ===\n\n";
}


}

