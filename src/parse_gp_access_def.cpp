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
#include <axs_gp_access_def.h>

#include <eom_config.h>

static void parse_elevation_constraint(std::deque<std::string>& cnst_toks);

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
  for (int ii=0; ii<n_constraints; ++ii) {
    parse_elevation_constraint(tokens);
    if (tokens.size() == 0) {
        break;
    }
  }

  eom::GpAccessDef gpDef(orbit_name, gp_name);
  return gpDef;

/*
  double rad_per_io {1.0/cfg.getIoPerRad()};
  double du_per_io {1.0/cfg.getIoPerDu()};
  if (coord_type == "LLA") {
    try {
      double lat {rad_per_io*std::stod(tokens[0])};
      tokens.pop_front();
      double lon {rad_per_io*std::stod(tokens[0])};
      tokens.pop_front();
      double alt {du_per_io*std::stod(tokens[0])};
      tokens.pop_front();
      eom::GroundPoint gp(lat, lon, alt);
      return std::make_pair(name, gp);
    } catch (const std::invalid_argument& ia) {
      throw std::invalid_argument("eom_app::parse_gp_access_def() "s +
                                  "invalid LLA parameter type"s);
    }
  }

  throw std::invalid_argument("eom_app::parse_gp_access_def() "s +
                              "Invalid coordinate type - "s + coord_type);

*/
}


}


// Add return of constraints structure
static void parse_elevation_constraint(std::deque<std::string>& cnst_toks)
{
  if (cnst_toks.size() > 1  &&  cnst_toks[0] == "MinimumElevation") {
    cnst_toks.pop_front();
      // Replace with double read and set constraints
    cnst_toks.pop_front();
  }
}

