/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_GROUND_POINT_H
#define ASTRO_GROUND_POINT_H

#include <cmath>

#include <Eigen/Dense>

namespace eom {

/**
 * Indicates which Starter method is employed by the Fukushima Cartesian
 * to geodetic algorithm.
 */
enum class FukStarter {
  none,                      ///< Fukushima lat to cart not used
  case1,                     ///< Most common, 99%
  case2,                     ///< Close to geocenter
  case3a,                    ///< Near equatorial a
  case3b                     ///< Near equatorial b
};

/**
 * An earth fixed point in Cartesian or geodetic coordinates.
 * Performs conversions between systems and supplies other useful
 * info.
 *
 * @author  Kurt Motekew
 * @date    2022/01/22
 */
class GroundPoint {
public:
  /**
   * Initialize with geodetic coordinates
   *
   * @param  lat  Geodetic latitude, radians
   * @param  lon  Longitude, radians
   * @param  alt  altitude, DU
   */
  GroundPoint(double lat, double lon, double alt);

  /**
   * Initialize with Cartesian earth fixed vector
   *
   * @param  xyz  Cartesian earth fixed coordinates, DU
   *
   * @throws  NonconvergenceException If somehow determination of
   *                                  geodetic latitude fails.
   */
  GroundPoint(const Eigen::Matrix<double, 3, 1>& xyz);

  /**
   * @return  Geodetic latitude, radians
   */
  double getLatitude() const noexcept
  {
    return m_lat;
  }

  /**
   * @return  Longitude, radians
   */
  double getLongitude() const noexcept
  {
    return m_lon;
  }

  /**
   * @return  Altitude above the ellipsoid, DU
   */
  double getAltitude() const noexcept
  {
    return m_alt;
  }

  /**
   * @return  Cartesian earth fixed vector, DU
   */
  Eigen::Matrix<double, 3, 1> getCartesian()
  {
    return m_xyz;
  }

  /**
   * @return  Number of iterations required to solve for geodetic
   *          latitude when initialized with Cartesian coordinates.
   *          A value of zero means the object was initialized
   *          with a geodetic coordinates.  A negative value indicates
   *          iteration failed to failed to converge ("should never
   *          happen").
   */
  int getItr() const noexcept
  {
    return itr;
  }
  
  /*
   * @return  Returns the Starter used by the Fukushima algorithm.  If
   *          non used (geodetic coord init), FukStarter::none is returned;
   */
  FukStarter getFukStarter()
  {
    return fstarter;
  }

private:
  double m_lat;
  double m_lon;
  double m_alt;
  Eigen::Matrix<double, 3, 1> m_xyz;

  int itr {0};
  FukStarter fstarter {FukStarter::none};
};


}

#endif

