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

static void parse_propagator(std::deque<std::string>&,
                             eom::PropagatorConfig&);
static void parse_gravity_model(std::deque<std::string>&,
                                eom::PropagatorConfig&);
static void parse_sun_model(std::deque<std::string>&,
                            eom::PropagatorConfig&);
static void parse_moon_model(std::deque<std::string>&,
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

  eom::CoordType coord_type {eom::CoordType::cartesian};
  eom::FrameType frame_type {eom::FrameType::gcrf};

  if (model == "SP"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::sp};
    propCfg.setStartStopTime(cfg.getStartTime(), cfg.getStopTime());
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
      // Loop through enough times to support finding all
      // supported options:
      //   1. Earth gravity model
      //   2. Sun gravity model
      //   3. Moon gravity model
      //   4. Integrator options
    unsigned int sp_options {4};
    for (unsigned int ii=0; ii<sp_options; ++ii) {
      parse_gravity_model(tokens, propCfg);
      parse_sun_model(tokens, propCfg);
      parse_moon_model(tokens, propCfg);
      parse_propagator(tokens, propCfg);
      if (tokens.size() == 0) {
        break;
      }
    }
    eom::OrbitDef orbit {name, propCfg, epoch, state, coord_type, frame_type};
    return orbit;
  } else if (model == "Kepler1"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::kepler1};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    eom::OrbitDef orbit {name, propCfg, epoch, state, coord_type, frame_type};
    return orbit;
  } else if (model == "KeplerMod"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::kepler1mod};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    eom::OrbitDef orbit {name, propCfg, epoch, state, coord_type, frame_type};
    return orbit;
  } else if (model == "Vinti6"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::vinti6};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    eom::OrbitDef orbit {name, propCfg, epoch, state, coord_type, frame_type};
    return orbit;
  } else if (model == "VintiJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::vinti_j2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    eom::OrbitDef orbit {name, propCfg, epoch, state, coord_type, frame_type};
    return orbit;
  } else if (model == "VintiMod"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::vinti6mod};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    eom::OrbitDef orbit {name, propCfg, epoch, state, coord_type, frame_type};
    return orbit;
#ifdef GENPL
  } else if (model == "SecJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::sec_j2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    eom::OrbitDef orbit {name, propCfg, epoch, state, coord_type, frame_type};
    return orbit;
  } else if (model == "OscJ2"  &&  tokens.size() > 0 ) {
    eom::PropagatorConfig propCfg {eom::PropagatorType::osc_j2};
    eom::JulianDate epoch = parse_datetime(tokens);
    std::array<double, 6> state = parse_state_vector(tokens, cfg, coord_type,
                                                                  frame_type);
    eom::OrbitDef orbit {name, propCfg, epoch, state, coord_type, frame_type};
    return orbit;
#endif
  }

  throw std::invalid_argument("eom_app::parse_orbit_def() "s +
                              "Invalid parse_orbit type: "s + model);
}

}


//
// Incomplete parsing due to failing expectations for a particular
// gravity model or propagator will lead to failure to clear all
// tokens by the end of the parsing loop above.  This will indicate
// an invalid option in the input file, or incorrect options
//

static void parse_propagator(std::deque<std::string>& prop_toks,
                             eom::PropagatorConfig& pCfg)
{
    // "Propagator  Method Dur dt"
  if (prop_toks.size() > 3  &&  prop_toks[0] == "Propagator") {
    prop_toks.pop_front();
    if (prop_toks[0] == "RK4") {
      prop_toks.pop_front();
      pCfg.setPropagator(eom::Propagator::rk4);
    } else if (prop_toks[0] == "Adams4") {
      prop_toks.pop_front();
      pCfg.setPropagator(eom::Propagator::adams4);
#ifdef GENPL
    } else if (prop_toks[0] == "GJ") {
      prop_toks.pop_front();
      pCfg.setPropagator(eom::Propagator::gj);
    } else if (prop_toks[0] == "GJs") {
      prop_toks.pop_front();
      pCfg.setPropagator(eom::Propagator::gjs);
#endif
    }
    pCfg.setStepSize(eom_app::parse_duration(prop_toks));
  }
}


static void parse_gravity_model(std::deque<std::string>& grav_toks,
                                eom::PropagatorConfig& pCfg)
{
    // Minimum size is currently 3:  "GravityModel Jn 2"
    // Maximum size is currently 4:  "GravityMode  XX 12 12"
  if (grav_toks.size() > 2  &&  grav_toks[0] == "GravityModel") {
    grav_toks.pop_front();
    if (grav_toks.size() > 1  &&  grav_toks[0] == "Jn") {
      grav_toks.pop_front();
      pCfg.setGravityModel(eom::GravityModel::jn);
      try {
        int degree {std::stoi(grav_toks[0])};
        grav_toks.pop_front();
        pCfg.setDegreeOrder(degree, 0);
      } catch (const std::invalid_argument& ia) {
        ;
      }
    } else if (grav_toks.size() > 2  &&  grav_toks[0] == "Standard") {
      grav_toks.pop_front();
      pCfg.setGravityModel(eom::GravityModel::std);
        // Not all tokens will be consumed if an error is thrown,
        // signaling a parsing error
      try {
        int degree {std::stoi(grav_toks[0])};
        grav_toks.pop_front();
        int order {std::stoi(grav_toks[0])};
        grav_toks.pop_front();
        pCfg.setDegreeOrder(degree, order);
      } catch (const std::invalid_argument& ia) {
        ;
      }
#ifdef GENPL
    } else if (grav_toks.size() > 2  &&  grav_toks[0] == "Gravt") {
      grav_toks.pop_front();
      pCfg.setGravityModel(eom::GravityModel::gravt);
      try {
        int degree {std::stoi(grav_toks[0])};
        grav_toks.pop_front();
        int order {std::stoi(grav_toks[0])};
        grav_toks.pop_front();
        pCfg.setDegreeOrder(degree, order);
      } catch (const std::invalid_argument& ia) {
        ;
      }
#endif
    }
  }
}


static void parse_sun_model(std::deque<std::string>& sun_toks,
                            eom::PropagatorConfig& pCfg)
{
    // "SunGravity  Model"
  if (sun_toks.size() > 1  &&  sun_toks[0] == "SunGravity") {
    sun_toks.pop_front();
    if (sun_toks[0] == "Meeus") {
      sun_toks.pop_front();
      pCfg.setSunGravityModel(eom::SunGravityModel::meeus);
    }
  }
}


static void parse_moon_model(std::deque<std::string>& moon_toks,
                             eom::PropagatorConfig& pCfg)
{
    // "MoonGravity  Model"
  if (moon_toks.size() > 1  &&  moon_toks[0] == "MoonGravity") {
    moon_toks.pop_front();
    if (moon_toks[0] == "Meeus") {
      moon_toks.pop_front();
      pCfg.setMoonGravityModel(eom::MoonGravityModel::meeus);
    }
  }
}
