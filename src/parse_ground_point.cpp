/*
 * Copyright 2021 Kurt Motekew
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

#include <eom_config.h>

namespace eom_app {

std::pair<std::string, eom::GroundPoint>
parse_ground_point(std::deque<std::string>& tokens, const EomConfig& cfg)
{
  using namespace std::string_literals;
    // Need at least the name, coord type, and coordinates
  if (tokens.size() < 5) {
    throw std::invalid_argument("eom_app::parse_ground_point() "s +
                                "5 tokens required vs. "s +
                                std::to_string(tokens.size()));
  }
  auto name = tokens[0];
  tokens.pop_front();
  auto coord_type = tokens[0];
  tokens.pop_front();

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
      throw std::invalid_argument("eom_app::parse_ground_point() "s +
                                  "invalid LLA parameter type"s);
    }
  }

  throw std::invalid_argument("eom_app::parse_ground_point() "s +
                              "Invalid coordinate type - "s + coord_type);

}


}

