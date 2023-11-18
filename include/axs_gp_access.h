/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef AXS_GP_ACCESS_H
#define AXS_GP_ACCESS_H

#include <string>
#include <memory>
#include <vector>

#include <cal_julian_date.h>
#include <astro_ground_point.h>
#include <astro_ephemeris.h>
#include <axs_gp_constraints.h>
#include <axs_interval.h>

namespace eom {


/**
 * Performs access analysis between a ground point and ephemeris
 * resource.  The ephemeris resource is assumed to be a valid orbit.
 *
 * @author  Kurt Motekew
 * @date    20231027
 */
class GpAccess {
public:
  /**
   * Initialize with ground point, ephemeris, and static constraints
   *
   * @param  gp   Ground point definition
   * @param  eph  Ephemeris source
   * @param  xcs  Access constraints
   */
  GpAccess(const JulianDate& jdStart,
           const JulianDate& jdStop,
           const GroundPoint& gp,
           const Ephemeris& eph,
           const GpConstraints& xcs);

  /**
   * @return  Name (string identifieer) associated with ground point
   */
  std::string getGpName() const noexcept { return m_gp_name; }

  /**
   * @return  Name (string identifieer) associated with orbit
   */
  std::string getOrbitName() const { return m_eph_name; }

  /**
   * Provides constant iterator access to access interval structures
   *
   * @return  Beginning of range iterator for axs_interval structure
   */
  std::vector<axs_interval>::const_iterator cbegin() const {
    return m_intervals.cbegin();
  }

  /**
   * Provides constant iterator access to access interval structures.
   *
   * @return  End of range iterator for axs_interval structure
   */
  std::vector<axs_interval>::const_iterator cend() const {
    return m_intervals.cend();
  }

  /**
   * Provides constant iterator access to access interval structures.
   * Note constant implementation.
   *
   * @return  Beginning of range iterator for axs_interval structure
   */
  std::vector<axs_interval>::const_iterator begin() const {
    return m_intervals.cbegin();
  }

  /**
   * Provides constant iterator access to access interval structures.
   * Note constant implementation.
   *
   * @return  End of range iterator for axs_interval structure
   */
  std::vector<axs_interval>::const_iterator end() const {
    return m_intervals.cend();
  }

private:
  std::string m_gp_name;
  std::string m_eph_name;

  std::vector<axs_interval> m_intervals;
};


}

#endif
