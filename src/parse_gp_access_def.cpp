/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_parse.h>

#include <array>
#include <string>
#include <deque>
#include <utility>
#include <stdexcept>

#include <astro_ground_point.h>
#include <axs_gp_constraints.h>
#include <axs_gp_access_def.h>

#include <eom_config.h>

static void parse_elevation_constraint(std::deque<std::string>& cnst_toks,
                                       const eom_app::EomConfig& cfg, 
                                       eom::GpConstraints& constraints);

namespace eom_app {

eom::GpAccessDef
parse_gp_access_def(std::deque<std::string>& tokens, const EomConfig& cfg)
{

  using namespace std::string_literals;
    // Need at least the name, coord type, and coordinates
  if (tokens.size() < 2) {
    throw std::invalid_argument("eom_app::parse_gp_access_def() "s +
                                "Minimum  of 2 tokens required vs. "s +
                                std::to_string(tokens.size()));
  }
  auto orbit_name = tokens[0];
  tokens.pop_front();
  auto gp_name = tokens[0];
  tokens.pop_front();

  int n_constraints {1};
  eom::GpConstraints xcs;
  for (int ii=0; ii<n_constraints; ++ii) {
    parse_elevation_constraint(tokens, cfg, xcs);
    if (tokens.size() == 0) {
        break;
    }
  }

  eom::GpAccessDef gpDef(orbit_name, gp_name, xcs);
  return gpDef;
}


}


// Add return of constraints structure
static void parse_elevation_constraint(std::deque<std::string>& cnst_toks,
                                       const eom_app::EomConfig& cfg, 
                                       eom::GpConstraints& constraints)
{
  using namespace std::string_literals;
  if (cnst_toks.size() > 1  &&  cnst_toks[0] == "MinimumElevation") {
    cnst_toks.pop_front();
    try {
      double rad_per_io {1.0/cfg.getIoPerRad()};
      constraints.setMinEl(rad_per_io*std::stod(cnst_toks[0]));
      cnst_toks.pop_front();
    } catch (const std::invalid_argument& ia) {
      throw std::invalid_argument("eom_app::parse_access_def() "s +
                                  "invalid Minimum Elevation"s);
    }
  }
}

