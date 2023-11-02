/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef AXS_GP_ACCESS_DEF_H
#define AXS_GP_ACCESS_DEF_H

#include <string>
#include <memory>

#include <astro_ephemeris.h>
#include <astro_ground_point.h>
#include <axs_gp_constraints.h>

namespace eom {


/**
 * Holds parameters defining an access analysis request between an orbit
 * and a ground point.  Instantiation only sets the names of resources,
 * which must be added afterwards.
 *
 * @author  Kurt Motekew
 * @date    20230906
 */
class GpAccessDef {
public:
  /**
   * Create access definition from an ephemeris source to a ground
   * point.
   *
   * @param  orbit_name  Name of orbit generating ephemeris
   * @param  gp_name     Name of ground point for which access is to be
   *                     generated
   * @param  xcs         Access constraints
   */
  GpAccessDef(const std::string& orbit_name,
              const std::string& gp_name,
              const GpConstraints& xcs) : m_orbit_name {orbit_name},
                                          m_gp_name {gp_name},
                                          m_xcs {xcs}
  {
  }

  /**
   * @return  Name of orbit for which access is generated
   */
  std::string getOrbitName() const noexcept { return m_orbit_name; }

  /**
   * @return  Name of ground point for which access is generated
   */
  std::string getGpName() const noexcept { return m_gp_name; }

  /**
   * @return  Static (not dynamic) access constraints
   */
  GpConstraints getConstraints() const noexcept { return m_xcs; }

  /**
   * @return  Ground point definition used for access analysis
   */
  GroundPoint getGroundPoint() const { return *m_gp; }

  /**
   * @return  Ephemeris source used for access analysis
   */
  std::shared_ptr<Ephemeris> getEphemeris() const { return m_eph; }

  /**
   * Set resources needed by access analysis utilities.
   *
   * @param  gp   Ground point
   * @param  eph  Ephemeris source
   */
  void setResources(const GroundPoint& gp,
                    const std::shared_ptr<Ephemeris>& eph);

private:
  std::string m_orbit_name;
  std::string m_gp_name;
  GpConstraints m_xcs;

  std::unique_ptr<GroundPoint> m_gp {nullptr};
  std::shared_ptr<Ephemeris> m_eph {nullptr};
};


}

#endif
