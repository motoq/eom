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

struct aux_gp_constraints {
  bool use_max_sun_el {false};
  double max_sun_el {0.0};
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
   * @param  axcs        Auxiliary access constraints
   * @param  mdl         Optional specifier of which algorithm to use
   *                     when locating and refining access intervals
   */
  GpAccessDef(const std::string& orbit_name,
              const std::string& gp_name,
              const GpConstraints& xcs,
              const aux_gp_constraints& axcs,
              AccessModel  mdl = AccessModel::std);

  /**
   * @return  Name of orbit for which access is to be generated
   */
  std::string getOrbitName() const noexcept;

  /**
   * @return  Name of ground point for which access is to be generated
   */
  std::string getGpName() const noexcept;

  /**
   * @return  The access model algorithm type to use
   */
  AccessModel getAccessModel() const noexcept;

  /**
   * @return  Static (not dynamic) access constraints
   */
  GpConstraints getConstraints() const noexcept;

  /**
   * @return  True if auxiliary constraints need to be built
   */
  bool useAuxConstraints() const noexcept;

  /**
   * @return  Structure of auxiliary constraint settings
   *          that will require furhter construction
   */
  aux_gp_constraints getAuxConstraints() const noexcept;

private:
  std::string m_orbit_name;
  std::string m_gp_name;
  GpConstraints m_xcs;
  aux_gp_constraints m_axcs;
  AccessModel m_model;
};


}

#endif
