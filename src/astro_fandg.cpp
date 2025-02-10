/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_fandg.h>

#include <cmath>
#include <memory>
#include <string>
#include <utility>

#include <Eigen/Dense>

#include <mth_util.h>
#include <cal_const.h>
#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_keplerian.h>

#include <iostream>

namespace {
  constexpr double f2inv {1.0/static_cast<double>(mth_util::factorial(2))};
  constexpr double f3inv {1.0/static_cast<double>(mth_util::factorial(3))};
  constexpr double f4inv {1.0/static_cast<double>(mth_util::factorial(4))};
  constexpr double f5inv {1.0/static_cast<double>(mth_util::factorial(5))};
  constexpr double f6inv {1.0/static_cast<double>(mth_util::factorial(6))};
  constexpr double f7inv {1.0/static_cast<double>(mth_util::factorial(7))};
  constexpr double f8inv {1.0/static_cast<double>(mth_util::factorial(8))};
  constexpr double f9inv {1.0/static_cast<double>(mth_util::factorial(9))};
    // Used only for d/dz
  constexpr double f10inv {1.0/static_cast<double>(mth_util::factorial(10))};
  constexpr double f11inv {1.0/static_cast<double>(mth_util::factorial(11))};

  constexpr double xtol {1.0e-10};
  constexpr int maxitr {15};
}

namespace eom {

FandG::FandG(const std::string& orbit_name,
             const JulianDate& epoch,
             const Eigen::Matrix<double, 6, 1>& xeci,
             std::shared_ptr<const EcfEciSys> ecfeciSys)
{
  m_name = orbit_name;
  m_jd0 = epoch;
  m_ecfeci = std::move(ecfeciSys);

  Keplerian kep {xeci};
  m_a = kep.getSemimajorAxis();
  m_r0_mag = xeci.block<3,1>(0,0).norm();

    // 2-body problem solved in plane spanned by r & v
  Eigen::Matrix<double, 3, 3> c_pi = kep.getEciToPerifocal();
  Eigen::Matrix<double, 3, 1> r0 = c_pi*xeci.block<3,1>(0,0);
  Eigen::Matrix<double, 3, 1> v0 = c_pi*xeci.block<3,1>(3,0);
  m_r0 = r0.block<2,1>(0,0);
  m_v0 = v0.block<2,1>(0,0);
    // Only need to convert from perifocal back to
    // computational ECI frame after initialization
  m_c_ip = c_pi.transpose();
}


std::string FandG::getName() const
{
  return m_name;
}


JulianDate FandG::getEpoch() const
{
  return m_jd0;
}


JulianDate FandG::getBeginTime() const
{
  return m_ecfeci->getBeginTime();
}


JulianDate FandG::getEndTime() const
{
  return m_ecfeci->getEndTime();
}


/*
 * Universal variable based solution to Kepler's Problem,
 * f and g implementation from Fundamentals of Astrodynamics,
 * Bate, Mueller, & White - simply known as BMW in the community...
 *
 * Implemented with elliptical orbits in mind only at this time as
 * eom is focused on earth gravity bound problems.
 */
Eigen::Matrix<double, 6, 1> FandG::getStateVector(const JulianDate& jd,
                                                   EphemFrame frame) const 
{
  // Units of DU and TU
  //   GM = 1.0

  const double dt {phy_const::tu_per_day*(jd - m_jd0)};
  const double ainv {1.0/m_a};

    // Universal variable - initial guess for elliptical is GM*dt/smaj
  double xx {dt*ainv};

    // Refine guess of x via Newton's method
  const double r0mag {m_r0_mag};
  const double r0dotv0 {m_r0.dot(m_v0)};
  bool converged {false};
  for (int ii=0; ii<maxitr; ++ii) {
    auto xx2 = xx*xx;
    auto zz = xx2*ainv;
    auto [cz, sz] = astro_fg_cands(zz);
    auto dtn = xx*xx2*sz + r0dotv0*xx2*cz + r0mag*xx*(1.0 - zz*sz);
    auto dtdx = xx2*cz + r0dotv0*xx*(1.0 - zz*sz) + r0mag*(1.0 - zz*cz);
    auto dx = (dt - dtn)/dtdx;
    xx += dx;
    if (std::abs(dx) < xtol) {
      //std::cerr << "\nConverged on iteration " << ii << "\n";
      converged = true;
      break;
    }
  }

  if (!converged) {
    std::cerr << "\n\nFandG didn't Converge\n";
  }

  auto xx2 = xx*xx;
  double zz {xx2*ainv};
  auto [cz, sz] = astro_fg_cands(zz);
  auto f = 1.0 - xx2*cz/r0mag;
  auto g = dt - xx*xx2*sz;
  Eigen::Matrix<double, 2, 1> r = f*m_r0 + g*m_v0;
  auto rmag = r.norm();
  auto gdot = 1.0 - xx2*cz/rmag;
  auto fdot = xx*(zz*sz - 1.0)/(r0mag*rmag);
  Eigen::Matrix<double, 2, 1> v = fdot*m_r0 + gdot*m_v0;

  Eigen::Matrix<double, 3, 1> rpef = Eigen::Matrix<double, 3, 1>::Zero();
  Eigen::Matrix<double, 3, 1> vpef = Eigen::Matrix<double, 3, 1>::Zero();
  rpef.block<2,1>(0,0) = r;
  vpef.block<2,1>(0,0) = v;
  Eigen::Matrix<double, 6, 1> xeci;
  xeci.block<3,1>(0,0) = m_c_ip*rpef;
  xeci.block<3,1>(3,0) = m_c_ip*vpef;

  if (frame == EphemFrame::ecf) {
    return m_ecfeci->eci2ecf(jd, xeci.block<3,1>(0,0), xeci.block<3,1>(3,0));
  }

  return xeci;
}


Eigen::Matrix<double, 3, 1> FandG::getPosition(const JulianDate& jd,
                                                EphemFrame frame) const
{

  Eigen::Matrix<double, 6, 1> xeci = this->getStateVector(jd, frame);

  if (frame == EphemFrame::ecf) {
    return m_ecfeci->eci2ecf(jd, xeci.block<3,1>(0,0));
  }
  return xeci.block<3,1>(0,0);
}

/*
 * BMW implementation including use of series expansion for small
 * values of z to avoid division by zero.
 */
std::pair<double, double> astro_fg_cands(double z)
{
  if (z > 0.1) {
    auto sqrtz = std::sqrt(z);
    return std::make_pair((1.0   - std::cos(sqrtz))/z,
                          (sqrtz - std::sin(sqrtz))/(z*sqrtz));
  }

  return std::make_pair(f2inv - z*(f4inv - z*(f6inv - z*f8inv)),
                        f3inv - z*(f5inv - z*(f7inv - z*f9inv)));
}

std::pair<double, double> astro_fg_dcands_dz(double z)
{
  auto [cz, sz] = astro_fg_cands(z);
  if (z > 0.1) {
    auto twozinv = 0.5/z;
    return std::make_pair(twozinv*(1.0 - z*sz - 2.0*cz),
                          twozinv*(cz - 3.0*sz));
  }

  return std::make_pair(f4inv + z*(2*f6inv - z*(3*f8inv - z*4*f10inv)),
                        f5inv + z*(2*f7inv - z*(3*f9inv - z*4*f11inv)));
}


}

