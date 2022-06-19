/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>

#include <Eigen/Dense>

#include <mth_cross.h>

#include <eom_test.h>

static void print_3cross(const Eigen::Matrix<double, 3, 2>&,
                         const Eigen::Matrix<double, 3, 1>&);

namespace eom_app {

void eom_test_cross()
{

  std::cout << "\n\n  === Test:  N-Dimensional Cross Product ===";


  Eigen::Matrix<double, 3, 2> uMat = Eigen::Matrix<double, 3, 2>::Zero();
  uMat(2,0) = 1.0;
  uMat(0,1) = 1.0;
  Eigen::Matrix<double, 3, 1> v = eom::cross<double,3>(uMat);
  print_3cross(uMat, v);
  uMat = Eigen::Matrix<double, 3, 2>::Zero();
  uMat(0,0) = 1.0;
  uMat(2,1) = 1.0;
  v = eom::cross<double,3>(uMat);
  std::cout << '\n';
  print_3cross(uMat, v);
  uMat = Eigen::Matrix<double, 3, 2>::Zero();
  uMat(0,0) = 1.0;
  uMat(1,1) = 1.0;
  v = eom::cross<double,3>(uMat);
  std::cout << '\n';
  print_3cross(uMat, v);

  std::cout << "\n  === End Test:  N-Dimensional Cross Product ===\n\n";

}


}


static void print_3cross(const Eigen::Matrix<double, 3, 2>& u,
                         const Eigen::Matrix<double, 3, 1>& v)
{
  for (int ii=0; ii<3; ++ii) {
    std::cout << "\n  " << u(ii,0) << " " << u(ii,1) << " " <<  v(ii);
  }
}
