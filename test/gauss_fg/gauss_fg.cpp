/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cmath>
#include <iomanip>
#include <iostream>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_duration.h>
#include <astro_generate.h>


/**
 * Short unit test for implementation of the universal variable approach
 * to Gauss' problem via the f & g implementation.  This is the
 * example problem from BMW.
 */
int main()
{
    // Initial and end positions, DU, ECI, and time between
  Eigen::Matrix<double, 3, 1> r1 {0.5, 0.6, 0.7};
  Eigen::Matrix<double, 3, 1> r2 = {0.0, 1.0, 0.0};
  eom::Duration dur {13.0, phy_const::tu_per_min};

  Eigen::Matrix<double, 6, 1> x1 = eom::generate_gauss_fg_xfer(r1, r2, dur);

  std::cout << "\n\nx1:\n" << x1;

    // Useful 2-body values
  std::cout.precision(8);
  constexpr double omega {phy_const::earth_angular_velocity(0.0)};
  constexpr double geosyn_alt {std::cbrt(1.0/(omega*omega))};
  std::cout << "\n\nGeosynchronous 2-body orbital altitude: " <<
               phy_const::er_per_du*geosyn_alt << " ER, or " <<
               phy_const::m_per_du*geosyn_alt << " m";
  double period = utl_const::tpi/phy_const::earth_angular_velocity(0.0);
  std::cout << "\nPeriod = " << phy_const::sec_per_tu*period << " seconds";

  std::cout << '\n';

}

