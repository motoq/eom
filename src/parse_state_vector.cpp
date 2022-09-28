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

namespace eom_app {

std::array<double, 6> parse_state_vector(std::deque<std::string>& tokens,
                                         const EomConfig& cfg)
{
  using namespace std::string_literals;
    // Need at least the name and type of orbit
  if (tokens.size() < 8) {
    throw std::invalid_argument("eom_app::parse_state_vector"s +
                                "  8 tokens required"s);
  }
  auto coord_type = tokens[0];
  tokens.pop_front();
  auto coord_sys = tokens[0];
  tokens.pop_front();

  double du_per_io {1.0/cfg.getIoPerDu()};
  double io_per_tu {cfg.getIoPerTu()};
  if (coord_type == "CART"  &&  coord_sys == "GCRF") {
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
      throw std::invalid_argument("eom_app::parse_state_vector"s +
                                  "  invalid parameter type"s);
    }
    return xeci;
  }

  throw std::invalid_argument("eom_app::parse_state_vector"s +
                              "  Invalid frame or system - "s +
                              coord_type + coord_sys);
}

}

