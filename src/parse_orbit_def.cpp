/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_parse.h>

#include <stdexcept>
#include <array>
#include <string>
#include <deque>

#include <cal_julian_date.h>
#include <astro_propagator_config.h>
#include <astro_orbit_def.h>

#include <eom_config.h>

namespace eom_app {

eom::OrbitDef parse_orbit_def(std::deque<std::string>& tokens,
                              const EomConfig& cfg)
{
  using namespace std::string_literals;
    // Need at least the name and type of orbit
  if (tokens.size() < 2) {
     throw std::invalid_argument("eom_app::parse_orbit_def:"s +
                                 "  Invalid number of tokens to parse_orbit"s);
  }
  auto name = tokens[0];
  tokens.pop_front();
  auto model = tokens[0];
  tokens.pop_front();

  if (model == "Kepler1"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::Kepler1};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> xeci = parse_state_vector(tokens, cfg);
    eom::OrbitDef orbit {name, propCfg, epoch, xeci,
                         eom::CoordType::cartesian, eom::FrameType::gcrf};
    return orbit;
  } else if (model == "Vinti6"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::Vinti6};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> xeci = parse_state_vector(tokens, cfg);
    eom::OrbitDef orbit {name, propCfg, epoch, xeci,
                         eom::CoordType::cartesian, eom::FrameType::gcrf};
    return orbit;
  } else if (model == "VintiJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::VintiJ2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> xeci = parse_state_vector(tokens, cfg);
    eom::OrbitDef orbit {name, propCfg, epoch, xeci,
                         eom::CoordType::cartesian, eom::FrameType::gcrf};
    return orbit;
#ifdef GENPL
  } else if (model == "SecJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::SecJ2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> xeci = parse_state_vector(tokens, cfg);
    eom::OrbitDef orbit {name, propCfg, epoch, xeci,
                         eom::CoordType::cartesian, eom::FrameType::gcrf};
    return orbit;
  } else if (model == "OscJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::OscJ2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> xeci = parse_state_vector(tokens, cfg);
    eom::OrbitDef orbit {name, propCfg, epoch, xeci,
                         eom::CoordType::cartesian, eom::FrameType::gcrf};
    return orbit;
#endif
  }

  throw std::invalid_argument("eom_app::parse_orbit_def:"s +
                              "  Invalid parse_orbit type: "s + model);
}

}

