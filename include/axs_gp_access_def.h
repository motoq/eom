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

#include <axs_gp_constraints.h>

namespace eom {

enum class AccessModel {
  std,                            ///< Standard model (robust, fairly quick)
  dbg                             ///< Debug model (very robust, slooow)
};

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
   * Create access definition from an ephemeris source to a ground
   * point.
   *
   * @param  orbit_name  Name of orbit generating ephemeris
   * @param  gp_name     Name of ground point for which access is to be
   *                     generated
   * @param  xcs         Access constraints
   * @param  mdl         Optional specifier of which algorithm to use
   *                     when locating and refining access intervals
   */
  GpAccessDef(const std::string& orbit_name,
              const std::string& gp_name,
              const GpConstraints& xcs,
              AccessModel  mdl = AccessModel::std) : m_orbit_name {orbit_name},
                                                     m_gp_name {gp_name},
                                                     m_xcs {xcs},
                                                     m_model {mdl}
  {
  }

  /**
   * @return  Name of orbit for which access is to be generated
   */
  std::string getOrbitName() const noexcept { return m_orbit_name; }

  /**
   * @return  Name of ground point for which access is to be generated
   */
  std::string getGpName() const noexcept { return m_gp_name; }

  /**
   * @return  The access model algorithm type to use
   */
  AccessModel getAccessModel() const noexcept { return m_model; }

  /**
   * @return  Static (not dynamic) access constraints
   */
  GpConstraints getConstraints() const noexcept { return m_xcs; }

private:
  std::string m_orbit_name;
  std::string m_gp_name;
  GpConstraints m_xcs;
  AccessModel m_model;
};


}

#endif
