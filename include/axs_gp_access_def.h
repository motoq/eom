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

namespace eom {


/**
 * Holds parameters defining an access analysis request between an orbit
 * and a ground point.
 *
 * @author  Kurt Motekew
 * @date    20230906
 */
class GpAccessDef {
public:
  /**
   * Create access definition
   *
   * @param  orbit_name  Name (string identifieer) associated with orbit
   * @param  gp_name     Name (string identifieer) associated ground point
   */
  GpAccessDef(const std::string& orbit_name,
              const std::string& gp_name) : m_orbit_name(orbit_name),
                                            m_gp_name(gp_name)
  {
  }

  /**
   * @return  Name (string identifieer) associated with orbit
   */
  std::string getOrbitName() const noexcept { return m_orbit_name; }

  /**
   * @return  Name (string identifieer) associated with ground point
   */
  std::string getGptName() const noexcept { return m_gp_name; }

private:
  std::string m_orbit_name {""};
  std::string m_gp_name {""};
};


}

#endif
