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

#include <eom_config.h>
#include <phy_const.h>

namespace eom_app {

std::array<double, 6> parse_state_vector(std::deque<std::string>& tokens,
                                                    const EomConfig& cfg)
{
    // Need at least the name and type of orbit
  if (tokens.size() < 8) {
    throw std::invalid_argument("8 tokens required by parse_state_vector");
  }
  auto coord_type = tokens[0];
  tokens.pop_front();
  auto coord_sys = tokens[0];
  tokens.pop_front();

  if (coord_type == "CART"  &&  coord_sys == "GCRF") {
    std::array<double, 6> xeci;
    try {
      for (unsigned int ii=0; ii<6; ++ii) {
        xeci[ii] = cfg.getToKilometers()*std::stod(tokens[0]);
        tokens.pop_front();
        if (ii > 2L) {
          cfg.getToSeconds();
          xeci[ii] /= cfg.getToSeconds();
        }
      }
    } catch(std::invalid_argument& ia) {
      throw std::invalid_argument("parse_state_vector: invalid parameter type");
    }
    for (unsigned int ii=0; ii<6; ++ii) {
      xeci[ii] *= phy_const::du_per_km;
      if (ii > 2L) {
        xeci[ii] *= phy_const::sec_per_tu;
      }
    }
    return xeci;
  }

  throw std::invalid_argument("Invalid parse_state_vector frame or system: " +
                                                        coord_type + coord_sys);
}

}

