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

#include <axs_gp_access_def.h>
#include <axs_gp_visibility.h>
#include <axs_gp_access.h>

#include <eom_config.h>

#include <eomx_exception.h>
#include <eomx.h>

/**
 * See eomx.h
 */
std::unordered_map<std::string,std::shared_ptr<eom::GpVisibility>>
eomx_gen_gp_accesses(
    const eom_app::EomConfig& cfg,
    const std::unordered_map<std::string,
                             std::shared_ptr<eom::GroundPoint>>& ground_points,
    const std::unordered_map<std::string,
                             std::shared_ptr<eom::Ephemeris>>& ephemerides,
    const std::vector<eom::GpAccessDef>& gp_access_defs)
{

    // Create access analysis objects with resources
    // Error if resource name is not available in existing containers
  std::unordered_map<std::string,
                     std::shared_ptr<eom::GpVisibility>> gp_accessors;
  for (auto& axs : gp_access_defs) {
    bool first {true};
    try {
      auto gp_ptr = ground_points.at(axs.getGpName());
      first = false;
      auto eph_ptr = ephemerides.at(axs.getOrbitName());
      gp_accessors[gp_ptr->getName() + eph_ptr->getName()] =
          std::make_shared<eom::GpAccess>(cfg.getStartTime(),
                                          cfg.getStopTime(),
                                          *gp_ptr,
                                          axs.getConstraints(),
                                          eph_ptr);
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

