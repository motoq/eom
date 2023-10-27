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

#include <astro_ground_point.h>
#include <astro_ephemeris.h>
#include <axs_gp_constraints.h>

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
  GpAccess(const GroundPoint& gp,
           const std::shared_ptr<eom::Ephemeris>& eph,
           const GpConstraints& xcs) : m_gp {gp},
                                       m_eph {eph},
                                       m_xcs {xcs}
  {
  }

  /**
   * @return  Name (string identifieer) associated with ground point
   */
  std::string getGpName() const noexcept { return m_gp.getName(); }

  /**
   * @return  Name (string identifieer) associated with orbit
   */
  std::string getOrbitName() const { return m_eph->getName(); }

private:
  GroundPoint m_gp;
  std::shared_ptr<eom::Ephemeris> m_eph;
  GpConstraints m_xcs;
};


}

#endif
