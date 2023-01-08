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
#include <stdexcept>

#include <eom_config.h>
#include <astro_orbit_def.h>

namespace eom_app {

std::array<double, 6> parse_state_vector(std::deque<std::string>& tokens,
                                         const EomConfig& cfg,
                                         eom::CoordType& coord_type,
                                         eom::FrameType& frame_type)
{
  using namespace std::string_literals;
    // coordinate type, frame, and 6 elements/components
  if (tokens.size() < 8) {
    throw std::invalid_argument("eom_app::parse_state_vector() "s +
                                "8 tokens required vs. "s +
                                std::to_string(tokens.size()));
  }
  if (tokens[0] == "CART") {
    coord_type == eom::CoordType::cartesian;
    tokens.pop_front();
  } else {
    throw std::invalid_argument("eom_app::parse_state_vector() "s +
                                "Invalid coordinate system type: "s +
                                tokens[0]);
  }
  if (tokens[0] == "GCRF") {
    frame_type = eom::FrameType::gcrf;
    tokens.pop_front();
  } else if (tokens[0] == "ITRF") {
    frame_type = eom::FrameType::itrf;
    tokens.pop_front();
  } else {
    throw std::invalid_argument("eom_app::parse_state_vector() "s +
                                "Invalid reference frame type: "s +
                                tokens[0]);
  }

  double du_per_io {1.0/cfg.getIoPerDu()};
  double io_per_tu {cfg.getIoPerTu()};
  std::array<double, 6> xeci;
  try {
    for (unsigned int ii=0; ii<6; ++ii) {
      xeci[ii] = du_per_io*std::stod(tokens[0]);
      tokens.pop_front();
      if (ii > 2L) {
        xeci[ii] *= io_per_tu;
      }
    }
  } catch (const std::invalid_argument& ia) {
    throw std::invalid_argument("eom_app::parse_state_vector() "s +
                                "invalid parameter type"s);
  }
  return xeci;
}

}

