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

namespace eom {


/**
 * Supported input coordinate frame (orbital element) types
 */
enum class CoordType {
  cartesian                       ///< position (DU), velocity (DU/TU)
};

/**
 * Supported input reference frames
 */
enum class FrameType {
  itrf,                           ///< WGS 84 compatible
  gcrf                            ///< GCRF (IAU 2000A/2006)
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
  OrbitDef(std::string orbit_name, const PropagatorConfig& propConfig,
           const JulianDate& epoch, const std::array<double, 6>& state,
                            CoordType coord_type, FrameType frame_type);

  /**
   * @return  Name (string identifieer) associated with orbit
   */
  std::string getOrbitName() const noexcept { return name; }

  /**
   * @return  Orbit propagator configuration
   */
  PropagatorConfig getPropagatorConfig() const noexcept { return propCfg; }

  /**
   * @return  Time of orbit state vector
   */
  JulianDate getEpoch() const noexcept { return jd0; }

  /**
   * @return  State vector, initial conditions at epoch time.  Units of
   *          DU, TU, and/or radians
   */
  std::array<double, 6> getInitialState() const noexcept { return x0; }

  /**
   * @return  State vector coordinate system type
   */
  CoordType getCoordinateType() const noexcept { return coord; }

  /**
   * @return  State vector eference frame
   */
  FrameType getReferenceFrameType() const noexcept { return frame; }

private:
  std::string name {""};
  PropagatorConfig propCfg;
  CoordType coord;
  FrameType frame;
  JulianDate jd0;
  std::array<double, 6> x0;
};


}

#endif
