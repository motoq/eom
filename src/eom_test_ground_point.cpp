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
static std::string fuk_start_str(eom::FukInit);

namespace eom_app {

void eom_test_ground_point()
{

  std::cout << "\n\n  === Test:  GroundPoint ===";

    // Unit tests - init with Cartesian

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

    // Test grid of lat/lon with random alt
  
  double max_lat_err {0.0};
  double max_lon_err {0.0};
  double max_alt_err {0.0};
  int npts {0};
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
        max_lat_err = std::max(max_lat_err,
                               std::fabs(lat - gpXyz.getLatitude()));
        max_lon_err = std::max(max_lon_err,
                               std::fabs(lon - gpXyz.getLongitude()));
        max_alt_err = std::max(max_alt_err,
                               std::fabs(alt - gpXyz.getAltitude()));
        alt += dalt;
        npts++;
      }
      lon += dlon;
    }
    lat += dlat;
  }
  std::cout << "\nMax Error Over Grid of " << npts << " points" <<
               "\n Lat: " << utl_const::deg_per_rad*max_lat_err << " (deg)" <<
               "   Lon: " << utl_const::deg_per_rad*max_lon_err << " (deg)" <<
               "   Alt: " << phy_const::m_per_du*max_alt_err << " m";


  std::cout << "\n  === End Test:  GroundPoint ===\n\n";
}


}


static void test_xyz(Eigen::Matrix<double, 3, 1>& xyz)
{
  eom::GroundPoint gp {xyz};
  std::cout << "\nnitr: " << gp.getItr();
  std::cout << " using " << fuk_start_str(gp.getFukInit());
  eom::GroundPoint gp2(gp.getLatitude(), gp.getLongitude(), gp.getAltitude());
  std::cout << "  Error: " <<
      phy_const::m_per_du*(xyz - gp2.getCartesian()).norm() << " m";
  
}


static std::string fuk_start_str(eom::FukInit starter)
{
  switch (starter) {
    case eom::FukInit::none:
      return "None";
    case eom::FukInit::case1:
      return "Case1";
    case eom::FukInit::case2:
      return "Case2";
    case eom::FukInit::case3a:
      return "Case3a";
    case eom::FukInit::case3b:
      return "Case3b";
    default:
      return "";
  }
}

