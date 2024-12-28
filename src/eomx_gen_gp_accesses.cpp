/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <string>
#include <utility>
#include <memory>
#include <vector>
#include <unordered_map>
#include <execution>

#include <astro_ecfeci_sys.h>
#include <astro_ground_point.h>
#include <axs_gp_access_def.h>
#include <axs_gp_access.h>
#include <axs_gp_access_std.h>
#include <axs_gp_access_debug.h>
#include <axs_gp_sun_constraint.h>

#include <eom_config.h>

#include <eomx_exception.h>
#include <eomx.h>

/**
 * See eomx.h
 */
std::unordered_map<std::string,std::shared_ptr<eom::GpAccess>>
eomx_gen_gp_accesses(
    const eom_app::EomConfig& cfg,
    const std::unordered_map<std::string,
                             std::shared_ptr<eom::GroundPoint>>& ground_points,
    const std::unordered_map<std::string,
                             std::shared_ptr<eom::Ephemeris>>& ephemerides,
    const std::vector<eom::GpAccessDef>& gp_access_defs,
    const std::shared_ptr<const eom::EcfEciSys>& f2iSys)
{

    // Create access analysis objects with resources
    // Error if resource name is not available in existing containers
  std::unordered_map<std::string,
                     std::shared_ptr<eom::GpAccess>> gp_accessors;
  for (const eom::GpAccessDef& axs : gp_access_defs) {
    bool first {true};
    try {
      auto gp_ptr = ground_points.at(axs.getGpName());
      first = false;
      auto eph_ptr = ephemerides.at(axs.getOrbitName());
      auto xcs = axs.getConstraints();
      if (axs.useAuxConstraints()) {
        eom::aux_gp_constraints axcs = axs.getAuxConstraints();
        if (axcs.use_max_sun_el) {
          auto sunx = std::make_shared<eom::GpSunConstraint>(*gp_ptr,
                                                             f2iSys);
          sunx->setMaxElevation(axcs.max_sun_el);
          xcs.addConstraint(sunx);
        }
      }
        // Select access determination algorithm
      if (axs.getAccessModel() == eom::AccessModel::dbg) {
        gp_accessors[gp_ptr->getName() + eph_ptr->getName()] =
            std::make_shared<eom::GpAccessDebug>(cfg.getStartTime(),
                                                 cfg.getStopTime(),
                                                 *gp_ptr,
                                                 xcs,
                                                 eph_ptr);
      } else {
        gp_accessors[gp_ptr->getName() + eph_ptr->getName()] =
            std::make_shared<eom::GpAccessStd>(cfg.getStartTime(),
                                               cfg.getStopTime(),
                                               *gp_ptr,
                                               xcs,
                                               eph_ptr);
      }
    } catch (const std::out_of_range& oor) {
      using namespace std::string_literals;
      if (first) {
        throw eom_app::EomXException("Error Finding GP Access Ground Point "s+
                                     oor.what());
      } else {
        throw eom_app::EomXException("Error Finding GP Access Ephemeris "s +
                                     oor.what());
      }
    }
  }

  //==>
    // Generate access times in parallel 
  std::for_each(std::execution::par,
                gp_accessors.begin(),
                gp_accessors.end(),
                [](auto& accessor) { accessor.second->findAllAccesses(); });
  //<==

  return gp_accessors;
}

