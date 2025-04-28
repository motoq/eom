/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_parse.h>

#include <array>
#include <deque>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include <cal_julian_date.h>
#include <astro_propagator_config.h>
#include <astro_orbit_def.h>

#include <eom_config.h>

namespace eom_app {

eom::OrbitDef parse_orbit_def(std::deque<std::string>& tokens,
                              const EomConfig& cfg,
                              const std::unordered_map<
                                  std::string, eom::PropagatorConfig>& pcfgs)
{
  using namespace std::string_literals;
    // Need at least the name and type of orbit
  if (tokens.size() < 2) {
     throw std::invalid_argument("eom_app::parse_orbit_def() "s +
                                 "Invalid number of tokens to parse_orbit: "s +
                                 std::to_string(tokens.size()));
  }
  auto name = tokens[0];
  tokens.pop_front();
  auto model = tokens[0];
  tokens.pop_front();

  eom::CoordType coord_type {eom::CoordType::cartesian};
  eom::FrameType frame_type {eom::FrameType::gcrf};

  if (model == "SP") {
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
      // Set default propagator configuration now that state vector read
    eom::PropagatorConfig propCfg {eom::PropagatorType::sp};
      // And replace if config is specified
    if (tokens.size() > 0  &&  tokens[0] == "Propagator") {
      tokens.pop_front();
      try {
        propCfg = pcfgs.at(tokens[0]);
        tokens.pop_front();
      } catch (const std::out_of_range& oor) {
        throw std::invalid_argument("eom_app::parse_orbit_def() "s +
                                    "Invalid PropagatorConfig: "s +
                                    tokens[0]);
      }
    }
    propCfg.setStartStopTime(cfg.getStartTime(), cfg.getStopTime());
    return eom::OrbitDef {name, propCfg, epoch, state, coord_type, frame_type};
  } else if (model == "Kepler1"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::kepler1};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    return eom::OrbitDef {name, propCfg, epoch, state, coord_type, frame_type};
  } else if (model == "KeplerMod"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::kepler1mod};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    return eom::OrbitDef {name, propCfg, epoch, state, coord_type, frame_type};
  } else if (model == "FandG"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::fandg};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    return eom::OrbitDef {name, propCfg, epoch, state, coord_type, frame_type};
  } else if (model == "Vinti6"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::vinti6};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    return eom::OrbitDef {name, propCfg, epoch, state, coord_type, frame_type};
  } else if (model == "VintiJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::vinti_j2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    return eom::OrbitDef {name, propCfg, epoch, state, coord_type, frame_type};
  } else if (model == "VintiMod"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::vinti6mod};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    return eom::OrbitDef {name, propCfg, epoch, state, coord_type, frame_type};
#ifdef GENPL
  } else if (model == "SecJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::sec_j2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    return eom::OrbitDef {name, propCfg, epoch, state, coord_type, frame_type};
  } else if (model == "OscJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::osc_j2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    return eom::OrbitDef {name, propCfg, epoch, state, coord_type, frame_type};
#endif
  }

  throw std::invalid_argument("eom_app::parse_orbit_def() "s +
                              "Invalid parse_orbit type: "s + model);
}


}
