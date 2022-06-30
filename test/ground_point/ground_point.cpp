/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <astro_ground_point.h>

#include <eom_test.h>

static void test_xyz(Eigen::Matrix<double, 3, 1>&);
static std::string fuk_start_str(eom::FukStarter);

int main()
{
  std::cout << "\n\n  === Test:  GroundPoint ===";

    // Unit tests - init with Cartesian.  Two pairs with symmetry over
    // the equatoriall plane, one exact equatorial, two near equatorial,
    // and one near polar, and two, one at each pole.

  Eigen::Matrix<double, 3, 1> pos1P {6524.834, 6862.875, 6448.296};
  pos1P *= phy_const::du_per_km;
  test_xyz(pos1P);

  Eigen::Matrix<double, 3, 1> pos1N {6524.834, 6862.875, -6448.296};
  pos1N *= phy_const::du_per_km;
  test_xyz(pos1N);

  Eigen::Matrix<double, 3, 1> pos2P {-5552.0, -2563.0, 3258.0};
  pos2P *= phy_const::du_per_km;
  test_xyz(pos2P);

  Eigen::Matrix<double, 3, 1> pos2N {-5552.0, -2563.0, -3258.0};
  pos2N *= phy_const::du_per_km;
  test_xyz(pos2N);

  Eigen::Matrix<double, 3, 1> posEq {1.0, 0.0, 0.0};
  test_xyz(posEq);
  std::cout << "    Equatorial Test";

  Eigen::Matrix<double, 3, 1> posEqP {1.0, 0.0, 1.0*phy_const::du_per_km};
  test_xyz(posEqP);
  std::cout << "    Equatorial +dx Test";

  Eigen::Matrix<double, 3, 1> posEqN {1.0, 0.0, -1.0*phy_const::du_per_km};
  test_xyz(posEqN);
  std::cout << "    Equatorial -dx Test";

  Eigen::Matrix<double, 3, 1> posNp { 1.0*phy_const::du_per_km, 
                                     -1.0*phy_const::du_per_km, 1.0};
  test_xyz(posNp);
  std::cout << "    Near Polar Test";

  Eigen::Matrix<double, 3, 1> posNorthP { 0.0, 0.0, 1.1*phy_const::du_per_km};
  test_xyz(posNorthP);
  std::cout << "    North Pole Test";

  Eigen::Matrix<double, 3, 1> posSouthP { 0.0, 0.0, -1.1*phy_const::du_per_km};
  test_xyz(posSouthP);
  std::cout << "    South Pole Test";

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


  std::cout << "\n  === End Test:  GroundPoint ===\n\n";

}


/*
 * Convenience for Cartesian to geodetic and back test
 */
static void test_xyz(Eigen::Matrix<double, 3, 1>& xyz)
{
  eom::GroundPoint gp {xyz};
  std::cout << "\nnitr: " << gp.getItr();
  std::cout << " using " << fuk_start_str(gp.getFukStarter());
  eom::GroundPoint gp2(gp.getLatitude(), gp.getLongitude(), gp.getAltitude());
  std::cout << "  Error: " <<
      phy_const::m_per_du*(xyz - gp2.getCartesian()).norm() << " m";
}


/*
 * @return string representation of Starter method
 */
static std::string fuk_start_str(eom::FukStarter starter)
{
  switch (starter) {
    case eom::FukStarter::none:
      return "None";
    case eom::FukStarter::case1:
      return "Starter 1";
    case eom::FukStarter::case2:
      return "Starter 2";
    case eom::FukStarter::case3a:
      return "Starter 3a";
    case eom::FukStarter::case3b:
      return "Starter 3b";
    default:
      return "";
  }
}

