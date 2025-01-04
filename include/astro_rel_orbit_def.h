/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_REL_ORBIT_DEF_H
#define ASTRO_REL_ORBIT_DEF_H

#include <array>
#include <string>

#include <cal_julian_date.h>

namespace eom {


/**
 * Supported relative orbit definition types
 *   RTCT is defined by a radial, transverse, and cross-track offset
 *   that is symmetric about the reference (chief) orbit.  The 4th
 *   parameter is a transverse offset that is not symmetric about the
 *   reference orbit.
 */
enum class RelCoordType {
  rtct                  ///< Radial, transverse, cross-track, transverse-offset
};


/**
 * Holds parameters defining an orbit based on that of another orbit
 * definition.  Successful creation of this object does not guarantee a
 * valid orbit definition.
 *
 * @author  Kurt Motekew
 * @date    2022/01/03
 */
class RelOrbitDef {
public:
  /**
   * Create an orbit definition based on a another orbit definition.
   *
   * @param  orbit_name     Name (string identifieer) associated with
   *                        new orbit
   * @param  template_name  Name of orbit used to define this orbit
   * @param  rel_state      Parameters defining the state of this orbit
   *                        relative to the template orbit.  Unlike a
   *                        state vector, the number of parameters may
   *                        vary depending on the coordinate type.
   * @param  coord_type     Parameters defining a relative orbit.
   *                        Typically state vector or orbital element
   *                        differences.  Differences are added to the
   *                        template orbit.
   */
  RelOrbitDef(const std::string& orbit_name,
              const std::string& template_name,
              const std::array<double, 6>& rel_state,
              RelCoordType coord_type);

  /**
   * @return  Name (string identifieer) associated with new orbit
   */
  std::string getOrbitName() const noexcept;

  /**
   * @return  Name of template orbit
   */
  std::string getTemplateOrbitName() const noexcept;

  /**
   * @return  Relative orbit parameters or state vector, DU, TU,
   *          and/or radians
   */
  std::array<double, 6> getInitialState() const noexcept;

  /**
   * @return  State vector coordinate system type
   */
  RelCoordType getRelCoordinateType() const noexcept;

private:
  std::string m_name {""};
  std::string m_template_name {""};
  RelCoordType m_coord;
  std::array<double, 6> m_dx0;
};


}

#endif
