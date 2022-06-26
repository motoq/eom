/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <cstdlib>
#include <ctime>

#include <Eigen/Dense>

#include <mth_cross_product.h>

#include <eom_test.h>

namespace {

template <typename T, unsigned int N>
T sumdot(Eigen::Matrix<T, N, N-1U>& mat, Eigen::Matrix<T, N, 1>& vec)
{
  T uv {0};
  for (unsigned int ii=0; ii<N-1U; ++ii) {
    uv += vec.transpose()*(mat.template block<N,1U>(0,ii));
  }
  return uv;
}

template <typename T>
T sumdot(Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& mat,
         Eigen::Matrix<T, Eigen::Dynamic, 1>& vec)
{
  T uv {0};
  for (unsigned int ii=0; ii<(mat.rows()-1); ++ii) {
    Eigen::Matrix<T, Eigen::Dynamic, 1> mcol = mat.block(0,ii,mat.rows(),1);
    uv += vec.dot(mcol);
  }
  return uv;
}

}

static void print_3cross(const Eigen::Matrix<double, 3, 2>&,
                         const Eigen::Matrix<double, 3, 1>&);

namespace eom_app {

void eom_test_cross()
{

  std::cout << "\n\n  === Test:  N-Dimensional Cross Product ===";


  eom::CrossProduct<double, 3> x3d;
  Eigen::Matrix<double, 3, 2> uMat = Eigen::Matrix<double, 3, 2>::Zero();
  uMat(2,0) = 1.0;
  uMat(0,1) = 1.0;
  Eigen::Matrix<double, 3, 1> v = x3d(uMat);
  print_3cross(uMat, v);
  uMat = Eigen::Matrix<double, 3, 2>::Zero();
  uMat(0,0) = 1.0;
  uMat(2,1) = 1.0;
  v = x3d(uMat);
  std::cout << '\n';
  print_3cross(uMat, v);
  uMat = Eigen::Matrix<double, 3, 2>::Zero();
  uMat(0,0) = 1.0;
  uMat(1,1) = 1.0;
  v = x3d(uMat);
  std::cout << '\n';
  print_3cross(uMat, v);

  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> ux = 
      Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(3,2);
  ux(2,0) = 1.0;
  ux(0,1) = 1.0;
  Eigen::Matrix<double, Eigen::Dynamic, 1> vx = eom::cross_product(ux);
  std::cout << '\n';
  print_3cross(ux, vx);

  ux = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(3,2);
  ux(0,0) = 1.0;
  ux(2,1) = 1.0;
  vx = eom::cross_product(ux);
  std::cout << '\n';
  print_3cross(ux, vx);

  ux = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>::Zero(3,2);
  ux(0,0) = 1.0;
  ux(1,1) = 1.0;
  vx = eom::cross_product(ux);
  std::cout << '\n';
  print_3cross(ux, vx);

  std::cout << '\n';
  std::srand((unsigned int) std::time(0));

  eom::CrossProduct<double, 2> x2d;
  Eigen::Matrix<double, 2, 1> u2x1 = Eigen::Matrix<double, 2, 1>::Random();
  Eigen::Matrix<double, 2, 1> v2 = x2d(u2x1);
  std::cout << "\n  Random 2Dd SumDot: " << sumdot<double, 2>(u2x1, v2);

  eom::CrossProduct<float, 2> x2f;
  Eigen::Matrix<float, 2, 1> u2x1f = Eigen::Matrix<float, 2, 1>::Random();
  Eigen::Matrix<float, 2, 1> v2f = x2f(u2x1f);
  std::cout << "\n  Random 2Df SumDot: " << sumdot<float, 2>(u2x1f, v2f);

  Eigen::Matrix<double, 3, 2> u3x2 = Eigen::Matrix<double, 3, 2>::Random();
  Eigen::Matrix<double, 3, 1> v3 = x3d(u3x2);
  std::cout << "\n  Random 3Dd SumDot: " << sumdot<double, 3>(u3x2, v3);

  eom::CrossProduct<float, 3> x3f;
  Eigen::Matrix<float, 3, 2> u3x2f = Eigen::Matrix<float, 3, 2>::Random();
  Eigen::Matrix<float, 3, 1> v3f = x3f(u3x2f);
  std::cout << "\n  Random 3Df SumDot: " << sumdot<float, 3>(u3x2f, v3f);

/*
  Eigen::Matrix<double, 4, 3> u4x3 = Eigen::Matrix<double, 4, 3>::Random();
  Eigen::Matrix<double, 4, 1> v4 = eom::cross<double,4>(u4x3);
  std::cout << "\n  Random 4Dd SumDot: " << sumdot<double, 4>(u4x3, v4);

  Eigen::Matrix<float, 4, 3> u4x3f = Eigen::Matrix<float, 4, 3>::Random();
  Eigen::Matrix<float, 4, 1> v4f = eom::cross<float,4>(u4x3f);
  std::cout << "\n  Random 4Df SumDot: " << sumdot<float, 4>(u4x3f, v4f);

  Eigen::Matrix<double, 5, 4> u5x4 = Eigen::Matrix<double, 5, 4>::Random();
  Eigen::Matrix<double, 5, 1> v5 = eom::cross<double,5>(u5x4);
  std::cout << "\n  Random 5Dd SumDot: " << sumdot<double, 5>(u5x4, v5);

  Eigen::Matrix<float, 5, 4> u5x4f = Eigen::Matrix<float, 5, 4>::Random();
  Eigen::Matrix<float, 5, 1> v5f = eom::cross<float,5>(u5x4f);
  std::cout << "\n  Random 5Df SumDot: " << sumdot<float, 5>(u5x4f, v5f);

  eom::CrossProduct<double, 6> x6d;
  Eigen::Matrix<double, 6, 5> u6x5 = Eigen::Matrix<double, 6, 5>::Random();
  Eigen::Matrix<double, 6, 1> v6 = x6d(u6x5);
  std::cout << "\n  Random 6Dd SumDot: " << sumdot<double, 6>(u6x5, v6);

  eom::CrossProduct<float, 6> x6f;
  Eigen::Matrix<float, 6, 5> u6x5f = Eigen::Matrix<float, 6, 5>::Random();
  Eigen::Matrix<float, 6, 1> v6f = x6f(u6x5f);
  std::cout << "\n  Random 6Df SumDot: " << sumdot<float, 6>(u6x5f, v6f);

  Eigen::Matrix<double, 8, 7> u8x7 = Eigen::Matrix<double, 8, 7>::Random();
  Eigen::Matrix<double, 8, 1> v8 = eom::cross<double,8>(u8x7);
  std::cout << "\n  Random 8Dd SumDot: " << sumdot<double, 8>(u8x7, v8);

  Eigen::Matrix<float, 8, 7> u8x7f = Eigen::Matrix<float, 8, 7>::Random();
  Eigen::Matrix<float, 8, 1> v8f = eom::cross<float,8>(u8x7f);
  std::cout << "\n  Random 8Df SumDot: " << sumdot<float, 8>(u8x7f, v8f);
*/

  Eigen::Matrix<double, 6, 5> u6x5 = Eigen::Matrix<double, 6, 5>::Random();
  Eigen::Matrix<float, 6, 5> u6x5f = Eigen::Matrix<float, 6, 5>::Random();
  Eigen::Matrix<double, 9, 8> u9x8 = Eigen::Matrix<double, 9, 8>::Random();

  ux = u6x5;
  vx = eom::cross_product<double>(ux);
  std::cout << "\n  Random 6Dd SumDot: " << sumdot<double>(ux, vx);

  Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> uxf = u6x5f;
  Eigen::Matrix<float, Eigen::Dynamic, 1> vxf =
                                          eom::cross_product<float>(uxf);
  std::cout << "\n  Random 6Df SumDot: " << sumdot<float>(uxf, vxf);

  ux = u9x8;
  vx = eom::cross_product<double>(ux);
  std::cout << "\n  Random 9Dd SumDot: " << sumdot<double>(ux, vx);

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
