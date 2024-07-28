/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_parse.h>

#include <deque>
#include <fstream>
#include <string>
#include <stdexcept>
#include <unordered_map>

#include <astro_ground_point.h>

namespace eom_app {

void parse_slr_snx_stations(
    std::deque<std::string>& tokens,
    std::unordered_map<std::string,
                       std::shared_ptr<eom::GroundPoint>>& ground_points)
{
  using namespace std::string_literals;
    // Need at least the name, coord type, and coordinates
  if (tokens.size() < 1) {
    throw std::invalid_argument("eom_app::parse_slr_snx_stations() "s +
                                "1 tokens required vs. "s +
                                std::to_string(tokens.size()));
  }
  auto file_name = tokens[0];
  tokens.pop_front();

  std::ifstream ifs(file_name);
  if (!ifs.is_open()) {
    throw std::runtime_error("parse_slr_snx_stations() Can't open "
                             + file_name);
  }
}


}

