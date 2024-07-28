/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_ORBIT_DEF_H
#define ASTRO_ORBIT_DEF_H

#include <string>
#include <array>

#include <cal_julian_date.h>
#include <astro_propagator_config.h>
#include <astro_tle.h>

namespace eom {


/**
 * Supported input coordinate frame (orbital element) types
 */
enum class CoordType {
  cartesian,                      ///< position (DU), velocity (DU/TU)
  keplerian,                      ///< Orbital elements with true anomaly
  tle                             ///< SGP compatible
};

/**
 * Supported input reference frames
 */
enum class FrameType {
  itrf,                           ///< WGS 84 compatible
  gcrf,                           ///< GCRF (IAU 2000A/2006)
  teme                            ///< True Equator, Mean Equinox
};


/**
 * Holds parameters defining an orbit that can be used to create a
 * source of ephemeris (initial conditions, propagator type/settings,
 * epoch).  Successful creation of this object does not guarantee a
 * valid orbit definition.
 *
 * @author  Kurt Motekew
 * @date    202111
 */
class OrbitDef {
public:
  /**
   * Create orbit definition
   *
   * @param  orbit_name  Name (string identifieer) associated with orbit
   * @param  propConfig  Orbit propagator configuration
   * @param  epoch       Time of orbit state vector
   * @param  state       State vector, initial conditions at epoch time.
   *                     Units of DU, TU, and/or radians
   * @param  coord_type  State vector coordinate system type
   * @param  frame_type  State vector eference frame
   */
  OrbitDef(const std::string& orbit_name,
           const PropagatorConfig& propConfig,
           const JulianDate& epoch,
           const std::array<double, 6>& state,
           CoordType coord_type,
           FrameType frame_type);

  /**
   * Create orbit definition based on an SGP4 TLE
   *
   * @param  orbit_name  Name (string identifieer) associated with orbit
   * @param  tle1        TLE line 1
   * @param  tle2        TLE line 2
   */
  OrbitDef(const std::string& orbit_name,
           const std::string& tle1,
           const std::string& tle2);

  /**
   * @return  Name (string identifieer) associated with orbit
   */
  std::string getOrbitName() const noexcept { return m_name; }

  /**
   * @return  Orbit propagator configuration
   */
  PropagatorConfig getPropagatorConfig() const noexcept { return m_propCfg; }

  /**
   * @return  Time of orbit state vector
   */
  JulianDate getEpoch() const noexcept { return m_jd0; }

  /**
   * @return  State vector, initial conditions at epoch time.  Units of
   *          DU, TU, and/or radians
   */
  std::array<double, 6> getInitialState() const noexcept { return m_x0; }

  /**
   * @return  State vector coordinate system type
   */
  CoordType getCoordinateType() const noexcept { return m_coord; }

  /**
   * @return  State vector eference frame
   */
  FrameType getReferenceFrameType() const noexcept { return m_frame; }

  /**
   * @return TLE if an SGP4 propagator
   */
  Tle getTle() const noexcept;

private:
  std::string m_name {""};
  PropagatorConfig m_propCfg;
  CoordType m_coord;
  FrameType m_frame;
  JulianDate m_jd0;
  std::array<double, 6> m_x0;
    // SGP4/TLE only
  Tle m_tle;
};


}

#endif
