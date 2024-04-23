/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_GROUND_POINT_H
#define ASTRO_GROUND_POINT_H

#include <utl_const.h>
#include <phy_const.h>
#include <obs_rng_az_sinel.h>

#include <Eigen/Dense>

#include <cmath>
#include <ostream>
#include <string>

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
   * @param  lat   Geodetic latitude, radians
   * @param  lon   Longitude, radians
   * @param  alt   altitude, DU
   * @param  name  Name of ground point, optional, defaults to empty
   *               string
   */
  GroundPoint(double lat, double lon, double alt,
              const std::string& name = "");

  /**
   * Initialize with Cartesian earth fixed vector
   *
   * @param  xyz   Cartesian earth fixed coordinates, DU
   * @param  name  Name of ground point, optional, defaults to empty
   *               string
   *
   * @throws  NonconvergenceException If somehow determination of
   *                                  geodetic latitude fails.
   */
  GroundPoint(const Eigen::Matrix<double, 3, 1>& xyz,
              const std::string& name = "");

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
  Eigen::Matrix<double, 3, 1> getCartesian() const noexcept
  {
    return m_xyz;
  }

  /**
   * @return  Name associated with this ground point.  Empty string
   *          if not set.
   */
  std::string getName() const noexcept
  {
    return m_name;
  }

  /**
   * Computes the sine of the elevation of an object w.r.t. the tangent
   * plane of the oblate spheroid at the location of this groud point.
   *
   * @param  pos  ECEF position, DU
   *
   * @return  Sine of elevation of pos w.r.t. this ground point
   */
  double getSinElevation(const Eigen::Matrix<double, 3, 1>& posF) const;

  /**
   * Computes the range, azimuth, and sine of the elevation of an object
   * w.r.t. the location of this groud point.
   *
   * @param  pos  ECEF position, DU
   *
   * @return  Range, az, el, DU, radians, and ND
   */
  rng_az_sinel<double>
  getRngAzSinEl(const Eigen::Matrix<double, 3, 1>& posF) const;

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
    return m_itr;
  }
  
  /*
   * @return  Returns the Starter used by the Fukushima algorithm.  If
   *          non used (geodetic coord init), FukStarter::none is returned;
   */
  FukStarter getFukStarter()
  {
    return m_fstarter;
  }

private:
    // Sets useful member variables after Cartesian and geodetic
  void finish();

    // Cartesian and geodetic ground point definitions
  std::string m_name;
  double m_lat;
  double m_lon;
  double m_alt;
  Eigen::Matrix<double, 3, 1> m_xyz;

    // Auxiliary variables of use set during construction
  double m_clat;
  double m_clon;
  double m_slat;
  double m_slon;

  int m_itr {0};
  FukStarter m_fstarter {FukStarter::none};
};


/**
 * Override output stream as non-member function.
 */
std::ostream& operator<<(std::ostream& out, const GroundPoint& gp);

}


#endif

