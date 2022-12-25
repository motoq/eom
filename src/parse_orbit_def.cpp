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
#include <astro_gravity_jn.h>

#include <eom_config.h>

static void parse_gravity_model(std::deque<std::string>&,
                                eom::PropagatorConfig&);

namespace eom_app {

eom::OrbitDef parse_orbit_def(std::deque<std::string>& tokens,
                              const EomConfig& cfg)
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

  if (model == "Sp"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::sp};
    propCfg.setStartStopTime(cfg.getStartTime(), cfg.getStopTime());
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> xeci = parse_state_vector(tokens, cfg);
      // Loop through enough times to support finding all
      // supported options - currently gravity model and 
      // integrator options.
    unsigned int sp_options {2};
    for (unsigned int ii=0; ii<sp_options; ++ii) {
      parse_gravity_model(tokens, propCfg);
    }
    eom::OrbitDef orbit {name, propCfg, epoch, xeci,
                         eom::CoordType::cartesian, eom::FrameType::gcrf};
    return orbit;
  } else if (model == "Kepler1"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::kepler1};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> xeci = parse_state_vector(tokens, cfg);
    eom::OrbitDef orbit {name, propCfg, epoch, xeci,
                         eom::CoordType::cartesian, eom::FrameType::gcrf};
    return orbit;
  } else if (model == "Vinti6"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::vinti6};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> xeci = parse_state_vector(tokens, cfg);
    eom::OrbitDef orbit {name, propCfg, epoch, xeci,
                         eom::CoordType::cartesian, eom::FrameType::gcrf};
    return orbit;
  } else if (model == "VintiJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::vinti_j2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> xeci = parse_state_vector(tokens, cfg);
    eom::OrbitDef orbit {name, propCfg, epoch, xeci,
                         eom::CoordType::cartesian, eom::FrameType::gcrf};
    return orbit;
#ifdef GENPL
  } else if (model == "SecJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::sec_j2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> xeci = parse_state_vector(tokens, cfg);
    eom::OrbitDef orbit {name, propCfg, epoch, xeci,
                         eom::CoordType::cartesian, eom::FrameType::gcrf};
    return orbit;
  } else if (model == "OscJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::osc_j2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> xeci = parse_state_vector(tokens, cfg);
    eom::OrbitDef orbit {name, propCfg, epoch, xeci,
                         eom::CoordType::cartesian, eom::FrameType::gcrf};
    return orbit;
#endif
  }

  throw std::invalid_argument("eom_app::parse_orbit_def() "s +
                              "Invalid parse_orbit type: "s + model);
}

}

// Incomplete parsing due to failing expectations for a particular
// gravity model will lead to failure to clear all tokens by the end of
// the parsing loop above
static void parse_gravity_model(std::deque<std::string>& grav_toks,
                                eom::PropagatorConfig& pCfg)
{
    // Minimum size is currently 3:  "GravityModel Jn 2"
    // Maximum size is currently 4:  "GravityMode  XX 12 12"
  if (grav_toks.size() > 2  &&  grav_toks[0] == "GravityModel") {
    grav_toks.pop_front();
    if (grav_toks.size() == 2  &&  grav_toks[0] == "Jn") {
      grav_toks.pop_front();
      try {
        int degree {std::stoi(grav_toks[0])};
        if (degree >= 0  &&  degree <= eom::GravityJn::getMaxDegree()) {
          grav_toks.pop_front();
          pCfg.setDegreeOrder(degree, 0);
        }
      } catch (const std::invalid_argument& ia) {
        ;
      }
    }
  }
}
