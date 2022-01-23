/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_ground_point.h>

#include <cmath>

#include <Eigen/Dense>

#include <utl_const.h>
#include <utl_math.h>
#include <phy_const.h>

namespace {
  constexpr int max_itr {10};
  constexpr double tol {1.0e-6};
  constexpr double ep {std::sqrt(1.0 - phy_const::ecc2)};
}

namespace eom {

// Standard geodetic to Cartesian using Radius of curvature 'n'
GroundPoint::GroundPoint(double lat, double lon, double alt)
{
    // Update this
  m_lat = lat;
  m_lon = lon;
  m_alt = alt;

  double slat {std::sin(lat)};
  double clat {std::cos(lat)};
  double n {phy_const::earth_smaj/
            std::sqrt(1.0 - phy_const::ecc2*slat*slat)};
  double nph {n + alt};
  m_xyz(0) = nph*clat*std::cos(lon);
  m_xyz(1) = nph*clat*std::sin(lon);
  m_xyz(2) = (n*(1.0 - phy_const::ecc2) + alt)*slat;
}


// The geodetic latitude method is based on an iterative method to
// maximize accuracy accross all conditions
//
// "Fast transform from geocentric to geodetic coordinates"
// T. Fukushima (National Astronomical Observatory, Japan)
// Journal of Geodesy (1999)
GroundPoint::GroundPoint(const Eigen::Matrix<double, 3, 1>& xyz)
{
    // Update this
  m_xyz = xyz;

  double ri {m_xyz(0)};
  double rj {m_xyz(1)};
  double rk {m_xyz(2)};
  
    // Check if on z-axis
  if (ri == 0.0  &&  rj == 0.0) {
    m_lon = 0.0;
    m_alt = std::fabs(rk) - phy_const::earth_smin;
    if (m_xyz(2) > 0.0) {
      m_lat = utl_const::pio2;
    } else {
      m_lat = -utl_const::pio2;
    }
  } else {
      // longitude
    m_lon = std::atan2(m_xyz(1), m_xyz(0));
      // latitude -  work solution in upper quadrant (positive z)
    auto lat_sgn = utl_math::sgn(rk);
    double p {std::sqrt(ri*ri + rj*rj)};
    double z {std::fabs(rk)};
    double zp {ep*z};
    double c {phy_const::earth_smaj*phy_const::ecc2};
    double u {2.0*(zp - c)};
    double v {2.0*(zp + c)};
    double tm {(c - zp)/p};
      // Get starting point for solution to quartic
      // case1 is the most common case for points on or above the
      // earth's surface (just below to well above the reference
      // ellipsoid).  For equatorial regions, case3b will get used.
    double t {1.0};
    if (tm <= 0.0) {
      t = (p - c + zp)/(p - c + 2*zp);
      fstarter = FukStarter::case1;
    } else if (tm >= 1.0) {
      t = p/(zp + c);
      fstarter = FukStarter::case2;
    } else {
      // 0 < tm < 1 is the equatorial region - evaluate quartic
      double fm {tm*(tm*tm*(p*tm + u) + v) - p};
      if (fm >= 0.0) {
        // tm >= 1 starting point
        t = p/(zp + c);
        fstarter = FukStarter::case3a;
      } else {
        // tm <= 0 starting point - this method sometimes gets used
        t = (p - c + zp)/(p - c + 2*zp);
        fstarter = FukStarter::case3b;
      }
    }
      // Newton's method
    itr = -1;
    for (int ii=0; ii<max_itr; ++ii) {
      double t2 {t*t};
      double num {p - (t*(t2*(p*t + u) + v))};
      double den {t2*(4.0*p*t + 3.0*u) + v};
      double dt {num/den};
      t += dt;
      if (std::fabs(dt) < tol) {
        itr = ii + 1;
        break;
      }
    }
      // Quartic zero t - now solve for geodetic latitude and altitude
    double t2 {t*t};
    double omt2 {1.0 - t2};
    double opt2 {1.0 + t2};
    m_lat = lat_sgn*std::atan2(omt2, 2.0*ep*t);
    m_alt = (2.0*p*ep*t + z*omt2 - phy_const::earth_smaj*ep*opt2)/
            std::sqrt(opt2*opt2 - 4.0*phy_const::ecc2*t2);
  }
}


}

