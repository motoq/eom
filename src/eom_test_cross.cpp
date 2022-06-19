/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <string>

#include <Eigen/Dense>

#include <mth_cross.h>

#include <eom_test.h>

namespace eom_app {

void eom_test_cross()
{

  std::cout << "\n\n  === Test:  N-Dimensional Cross Product ===";


  Eigen::Matrix<double, 3, 2> uMat = Eigen::Matrix<double, 3, 2>::Zero();
  uMat(2,0) = 1.0;
  uMat(0,1) = 1.0;

  std::cout << '\n' << uMat;
  
  Eigen::Matrix<double, 3, 1> v = eom::cross<double,3>(uMat);

  std::cout << '\n' << v;

  std::cout << "\n  === End Test:  N-Dimensional Cross Product ===\n\n";
}


}

