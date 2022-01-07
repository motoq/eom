/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_parse.h>

#include <stdexcept>
#include <string>
#include <array>
#include <deque>

#include <astro_rel_orbit_def.h>

#include <eom_config.h>

namespace eom_app {

eom::RelOrbitDef parse_rel_orbit_def(std::deque<std::string>& tokens,
                                     const EomConfig& cfg)
{
  using namespace std::string_literals;
    // Need at least the name, template name, and type of definition
  if (tokens.size() < 3) {
     throw std::invalid_argument("eom_app::parse_rel_orbit_def:"s +
                 "  Invalid number of tokens to parse_rel_orbit"s);
  }
  auto name = tokens[0];
  tokens.pop_front();
  auto template_name = tokens[0];
  tokens.pop_front();
  auto frame = tokens[0];
  tokens.pop_front();

  if (frame == "RTCT"  &&  tokens.size() > 3 ) {
    std::array<double, 6> dx;
    dx[4] = 0.0;
    dx[5] = 0.0;
    try {
      double du_per_io {1.0/cfg.getIoPerDu()};
      for (int ii=0; ii<4; ++ii) {
        dx[ii] = du_per_io*std::stod(tokens[0]);
        tokens.pop_front();
      }
    } catch(std::invalid_argument& ia) {
      throw std::invalid_argument("eom_app::parse_rel_orbit_def"s +
                                  "  invalid relative orbit parameter type"s);
    }
    eom::RelOrbitDef orbit {name, template_name, dx, eom::RelCoordType::rtct};
    return orbit;
  }

  throw std::invalid_argument("eom_app::parse_rel_orbit_def:"s +
                              "  Invalid relative orbit type type: "s + frame);
}

}

