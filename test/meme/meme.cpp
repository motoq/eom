/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iomanip>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <string>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_duration.h>
#include <cal_greg_date.h>
#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>


/**
 * Short unit test for J2000 to GCRF and MEME to GCRF
 */
int main()
{
    // Frame bias matrix from Vallado (J2000 to GCRF)
  Eigen::Matrix<double, 3, 3> fb;
  fb << 0.9999999999999942, 0.0000000707827948, -0.0000000805621738,
       -0.0000000707827974, 0.9999999999999969, -0.0000000330604088,
        0.0000000805621715, 0.0000000330604145,  0.9999999999999962;

    // Create ECFECI service
  eom::GregDate gd(2016, 3, 9);
  eom::JulianDate jd1(gd);
  eom::Duration dur(1.0, phy_const::tu_per_day);
  auto jd2 = jd1 + dur;
  eom::EcfEciSys f2i(jd1, jd2, dur, nullptr);

  Eigen::Matrix<double, 3, 1> r_gcrf {-5552.0, -2563.0, 3258.0};
  Eigen::Matrix<double, 3, 1> r_j2000_sofa = f2i.gcrf2j2000(r_gcrf);
  Eigen::Matrix<double, 3, 1> r_j2000_faa = fb.transpose()*r_gcrf;

  std::cout << '\n';

  std::cout << std::setprecision(10);
  std::cout << "\nGCRF:  " << r_gcrf;
  std::cout << "\nJ2000: " << r_j2000_sofa;
  std::cout << "\nJ2000: " << r_j2000_faa;
  std::cout << "\nGCRF vs. J2000:  " << (r_j2000_sofa - r_gcrf).norm();
  std::cout << "\nSOFA vs. FAA:    " << (r_j2000_sofa - r_j2000_faa).norm();

  std::cout << '\n';

}

