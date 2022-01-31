/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_EARTH_SURF_H
#define ASTRO_EARTH_SURF_H

#include <cmath>

#include <Eigen/Dense>

namespace eom {

/**
 * Indicates if the horizon point should be returned if there is no
 * intersection.  Also, with horizon_always, the horizon point will be
 * returned regardless of the pointing vector intersecting the earth.
 */
enum class HorizonMode {
  horizon_never,             ///< Never compute the horizon point
  horizon_miss,              ///< Compute horizon when no earth intersection
  horizon_always             ///< Always return a horizon point
};

/**
 * Computes either the intersection or horizon point to the earth's
 * oblate spheroid definition given a reference point outside of the
 * ellipsoid and a pointing vector.  The horizon point will be in the
 * plane formed by the position vector and the pointing vector, tangent
 * to the oblate spheroid.
 *
 * @author  Kurt Motekew, "Pointing Vector to Oblate Spheroid", July 24, 2014
 * @date    2022/01/29
 */
class EarthSurf {
public:
  /**
   * Initialize with desired horizon mode and optional altitude above
   * the ellipsoid
   *
   * @param  mode  Horizon mode
   * @param  alt   Altitude above ellipsoid, DU, default to zero
   */
  EarthSurf(HorizonMode mode, double alt = 0.0);

  /**
   * @param  pos  Position from which pointing vector emanates, ECEF, DU
   * @param  pnt  Pointing vector.
   */
  void setEarthSurf(const Eigen::Matrix<double, 3, 1>& pos,
                    const Eigen::Matrix<double, 3, 1>& pnt);

  /**
   * Indicates if the call to setEarthSurf resulted in an earth
   * intersection when in horizon modes horizon_never and horizon_miss.
   * When in horizon_always mode, true is always reterned.
   *
   * @return  If true, a horizon point is returned by getEarthSurf().
   */
  bool missed() const noexcept
  {
    return !earth_x;
  }

  /**
   * @return  Either earth intersection or horizon point, DU
   */
  Eigen::Matrix<double, 3, 1> getEarthSurf() const noexcept
  {
    return m_xyz;
  }

private:
  HorizonMode m_mode {HorizonMode::horizon_never};
  double a2 {1.0};
  double b2 {1.0};
  Eigen::Matrix<double, 3, 1> m_xyz;
  bool earth_x {false};
};

}

#endif

