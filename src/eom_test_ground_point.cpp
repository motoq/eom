/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <astro_ground_point.h>

#include <eom_test.h>

static void test_xyz(Eigen::Matrix<double, 3, 1>&);

namespace eom_app {

void eom_test_ground_point()
{

  std::cout << "\n\n  === Test:  GroundPoint ===";

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

  std::cout << "\n  === End Test:  GroundPoint ===\n\n";
}


}


static void test_xyz(Eigen::Matrix<double, 3, 1>& xyz)
{
  eom::GroundPoint gp {xyz};
  std::cout << "\nnitr: " << gp.getItr();
  eom::GroundPoint gp2(gp.getLatitude(), gp.getLongitude(), gp.getAltitude());
  std::cout << "  Error: " <<
      phy_const::m_per_du*(xyz - gp2.getCartesian()).norm() << " m";
  
}
