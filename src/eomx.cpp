/*
 * Copyright 2021, 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include <deque>
#include <unordered_map>
#include <execution>
#include <stdexcept>

#include <Eigen/Dense>

#include <eom_config.h>
#include <eom_parse.h>
#include <eom_command.h>
#include <eom_command_builder.h>
#include <eom_test.h>
#include <astro_orbit_def.h>
#include <astro_rel_orbit_def.h>
#include <astro_ephemeris_file.h>
#include <astro_ephemeris.h>
#include <astro_eop_sys.h>
#include <astro_ecfeci_sys.h>
#include <astro_build.h>
#include <astro_print.h>

#include <utl_const.h>
#include <phy_const.h>
#include <astro_keplerian.h>

/**
 * Equations of Motion:  An application focused on astrodynamics related
 * problems.  This program parses an input file building models and
 * commands to be applied to those models.
 *
 * @author  Kurt Motekew
 */
int main(int argc, char* argv[])
{
    // Check for filename
  if (argc < 2  ||  argc > 3) {
    std::cerr << "\nProper use is:  " << argv[0] << " <input_file_name> or";
    std::cerr << "\n                " << argv[0] << " <input_file_name>" <<
                                                    " <eop__file_name>\n";
    return 0;
  }
    // Try to open for input
  std::ifstream ifs(argv[1]);
  if (!ifs.is_open()) {
    std::cerr << "\nError opening " << argv[1] << '\n';
    return 0;
  }
  std::cout << "\nOpened " << argv[1] << '\n';

  //
  // Parse input file and generate the simulation configuration
  // parameters along with modeling component definitions (that will be
  // used to create the actual modeling components) and commands to be
  // applied to those models.
  //

    // General configuration parameter for the simulation
  eom_app::EomConfig cfg;
    // Orbit definitions, used to initialize propagators and/or generate
    // classes with buffered ephemeris
  std::vector<eom::OrbitDef> orbit_defs;
    // Orbit definitions based on other orbits.  As with orbit_defs,
    // will be used to initialize propagators and/or generate classes
    // with buffered ephemeris
  std::vector<eom::RelOrbitDef> rel_orbit_defs;
    // Ephemeris file definitions - not necessarily an orbit
  std::vector<eom::EphemerisFile> eph_file_defs;
    // Ephemeris objects.  Only the pointer is needed during parsing so
    // the source can be added to objects requiring ephemerides.  Ephemeris
    // objects are created after orbit_defs and rel_orbit_defs are
    // generated.
  const auto ephemerides =
      std::make_shared<std::unordered_map<std::string,
                                          std::shared_ptr<eom::Ephemeris>>>();
    // Earth fixed points
  std::unordered_map<std::string, eom::GroundPoint> ground_points;
    // A bucket of resources allowing for parsing and building of
    // commands to be applied to models during the simulation
  eom_app::EomCommandBuilder cmdBuilder(ephemerides);
    // The commands populated by cmdBuilder
  std::vector<std::shared_ptr<eom_app::EomCommand>> commands;
    // Read each line and pass to parser while tracking line number
    // Keep track of line number for error messages
  int line_number {0};
  std::string input_line;
  bool parse_tokens {false};
  bool input_error = {false};
  std::deque<std::string> tokens;
  while (std::getline(ifs,input_line)) {
    line_number++;
    std::string other_error {""};
    std::istringstream iss(input_line);
    std::string token;
    while (iss >> token  &&  !input_error) {
      if (token.front() == '#') {
        break;
      } else {
        if (token.back() == ';') {
          parse_tokens = true;
          if (token.size() > 1) {
            token.pop_back();
            tokens.push_back(token);
          }
        } else {
          tokens.push_back(token);
        }
        if (parse_tokens) {
          input_error = true;
          parse_tokens = false;
          if (tokens.size() > 0) {
            auto make = tokens[0];
            tokens.pop_front();
            // Start Input Types - cfg updates attempt to handle
            // exceptions.  Other inputs may require try blocks (see
            // "Orbit" as an example).
            if (make == "SimStart") {
              cfg.setStartTime(tokens);
              input_error = !cfg.isValid();
            } else if (make == "SimDuration") {
              cfg.setDuration(tokens);
              input_error = !cfg.isValid();
            } else if (make == "LeapSeconds") {
              cfg.setLeapSeconds(tokens);
              input_error = !cfg.isValid();
            } else if (make == "EcfEciRate") {
              cfg.setEcfEciRate(tokens);
              input_error = !cfg.isValid();
            } else if (make == "end") {
              input_error = false;
              ifs.seekg(0, std::ios::end);
            } else if (make == "AngleUnits") {
              cfg.setIoPerRad(tokens);
              input_error = !cfg.isValid();
            } else if (make == "DistanceUnits") {
              cfg.setIoPerDu(tokens);
              input_error = !cfg.isValid();
            } else if (make == "TimeUnits") {
              cfg.setIoPerTu(tokens);
              input_error = !cfg.isValid();
            } else if (make == "OutputRate") {
              cfg.setOutputRate(tokens);
              input_error = !cfg.isValid();
            } else if (make == "Orbit") {
              try {
                orbit_defs.push_back(eom_app::parse_orbit_def(tokens, cfg));
                input_error = false;
              } catch (const std::invalid_argument& ia) {
                std::string xerror = ia.what();
                other_error = "Invalid Orbit definition: " + xerror;
              }
            } else if (make == "RelativeOrbit") {
              try {
                rel_orbit_defs.push_back(eom_app::parse_rel_orbit_def(tokens,
                                                                      cfg));
                input_error = false;
              } catch (const std::invalid_argument& ia) {
                std::string xerror = ia.what();
                other_error = "Invalid Relative Orbit definition: " + xerror;
              }
            } else if (make == "EphemerisFile") {
              try {
                eph_file_defs.push_back(eom_app::parse_eph_file_def(tokens));
                input_error = false;
              } catch (const std::invalid_argument& ia) {
                std::string xerror = ia.what();
                other_error = "Invalid Ephemeris File definition: " + xerror;
              }
            } else if (make == "GroundPoint") {
              try {
                ground_points.insert(eom_app::parse_ground_point(tokens, cfg));
                input_error = false;
              } catch (const std::invalid_argument& ia) {
                std::string xerror = ia.what();
                other_error = "Invalid Ground Point definition: " + xerror;
              }
            } else if (make == "Command") {
              try {
                commands.push_back(cmdBuilder.buildCommand(tokens, cfg));
                input_error = false;
              } catch (const std::invalid_argument& ia) {
                std::string xerror = ia.what();
                other_error = "Invalid Command definition: " + xerror;
              }
            } else if (make == "Test") {
              try {
                eom_app::eom_test(tokens);
                input_error = false;
              } catch (const std::invalid_argument& ia) {
                std::string xerror = ia.what();
                other_error = "Invalid Test type: " + xerror;
              }
            } else {
              other_error = "Invalid input line type: " + make;
            }
            // End Input Types
          }
          if (tokens.size() > 0) {
            input_error = true;
            other_error += "\n  Did not use all tokens in last input record";
          }
        }
      }
    }
    if (input_error) {
      std::cerr << "\n\nError on line: " << line_number;
      std::cerr << '\n' << cfg.getError();
      std::cerr << "\nOther Error: " << other_error;
      std::cerr << '\n';
      break;
    }
  }
  ifs.close();
  if (tokens.size() > 0  &&  !input_error) {
    std::cerr << "\n\n=== Warning: Reached EOF non-empty que ===";
    std::cerr << "\n        (Probably left out a ';')";
    std::cerr << "\n        (Number of Tokens " << tokens.size() << " )";
    std::cerr << "\n        (Current Token " << tokens[0] << " )\n";
  }
  if (input_error) {
    std::cerr << "\n\nExiting on input error that can't be specified\n";
    return 0;
  }

    // Before moving on, check if all necessary template
    // orbits needed for the relative orbits exists
  for (const auto& relOrbit : rel_orbit_defs) {
    bool found {false};
    for (const auto& orbit : orbit_defs) {
      if (orbit.getOrbitName() == relOrbit.getTemplateOrbitName()) {
        found = true;
        break;
      }
    }
    if (!found) {
      std::cerr << "\nBad Relative Orbit Template Name: " <<
                   relOrbit.getTemplateOrbitName() << " - Exiting\n";
      return 0;
    }
  }

  //
  // Parsing complete - print scenario and generate models and services
  //

  cfg.print(std::cout);

    // Determine time span that must be supported by the simulation
    // based on the input scenario time and orbit epoch times.
  eom::JulianDate minJd = cfg.getStartTime();
  eom::JulianDate maxJd = cfg.getStopTime();
  for (const auto& orbit : orbit_defs) {
    if (orbit.getEpoch() < minJd) {
      minJd = orbit.getEpoch();
    }
    if (maxJd < orbit.getEpoch()) {
      maxJd = orbit.getEpoch();
    }
      // Backwards propagation for SP methods not currently supported
    if (orbit.getPropagatorConfig().getPropagatorType() ==
                                  eom::PropagatorType::sp) {
      if (!(orbit.getEpoch() - cfg.getStartTime()  <  phy_const::epsdt_days)) {
        std::cerr << "\n\nError:  SP orbit eopch for " <<
                      orbit.getOrbitName() <<
                     " must occur on or before the simulation start time.";
        std::cerr << "\nExiting\n";
        return 0;
      }
    }
  }

    // Ecf to Eci transformation service - immutable - pass as:
    //   const std::shared_ptr<const EcfEciSys>&
  std::shared_ptr<eom::EopSys> eopSys = nullptr;
  if (argc > 2) {
   eopSys = std::make_shared<eom::EopSys>(argv[2], minJd, maxJd);
  }
  auto f2iSys = std::make_shared<eom::EcfEciSys>(minJd,
                                                 maxJd,
                                                 cfg.getEcfEciRate(),
                                                 eopSys);

    // Parse interpolated ephemeris from files and
    // process sequentially
  for (const auto& ephFileDef : eph_file_defs) {
    (*ephemerides)[ephFileDef.getName()] = 
        eom::build_ephemeris(ephFileDef, cfg.getStartTime(),
                                         cfg.getStartTime(),
                                         f2iSys);
  }

  {//==>
    // Generate orbit definitions in parallel 
  std::vector<std::unique_ptr<eom::Ephemeris>> ephvec(orbit_defs.size());
  std::transform(std::execution::par,
                 orbit_defs.begin(), orbit_defs.end(), ephvec.begin(),
                 [f2iSys](const auto& orbit) {
                   return eom::build_orbit(orbit, f2iSys);
                 }
  );
    // Move ephemerides from temporary vector to ephemeris map
  for (unsigned int ii=0; ii<ephvec.size(); ++ii) {
    auto name = ephvec[ii]->getName();
    (*ephemerides)[name] = std::move(ephvec[ii]);
  }
  }//<==


  {//==>
    // Construct relative orbits - generate and append ephemerides.
    // Relative orbit definitions are based on primary orbit
    // definitions, not other relative orbit definitions (only
    // orbit_defs, not other rel_orbit_defs).
  std::vector<std::unique_ptr<eom::Ephemeris>> ephvec(rel_orbit_defs.size());
  std::transform(std::execution::par,
                 rel_orbit_defs.begin(), rel_orbit_defs.end(), ephvec.begin(),
                 [f2iSys, ephemerides, &orbit_defs](const auto& relOrbit) {
        // Find reference orbit - template names already validated
      std::unique_ptr<eom::Ephemeris> eph = nullptr;
      for (const auto& templateOrbit : orbit_defs) {
        if (templateOrbit.getOrbitName() == relOrbit.getTemplateOrbitName()) {
          std::shared_ptr<eom::Ephemeris> templateEph =
                               ephemerides->at(templateOrbit.getOrbitName());
          eph = eom::build_orbit(relOrbit, templateOrbit, templateEph, f2iSys);
        }
      }
      return eph;
    }
  );
    // Move ephemerides from temporary vector to ephemeris map
  for (unsigned int ii=0; ii<ephvec.size(); ++ii) {
    auto name = ephvec[ii]->getName();
    (*ephemerides)[name] = std::move(ephvec[ii]);
  }
  }//<==

    // Print derived orbits
  std::cout << '\n';
  for (const auto& relOrbit : rel_orbit_defs) {
    std::cout << "\n  " << relOrbit.getOrbitName() <<
                 "  derived from:  " <<
                 relOrbit.getTemplateOrbitName();
  }
    // Print all orbits as orbital elements
  std::cout << '\n';
  for (const auto& [name, eph] : *ephemerides) {
    std::cout << "\n  " << name;
    std::cout << "\n  " << eph->getEpoch().to_str();
    eom::Keplerian oeCart(eph->getStateVector(eph->getEpoch(),
                                              eom::EphemFrame::eci));
    oeCart.print(std::cout);
  }

    // Print ground points
  for (auto& nm_gp : ground_points) {
    std::cout << "\n  " << nm_gp.first;
    nm_gp.second.print(std::cout);
  }


  //
  // Model and command lists completed - no further modifications
  // Validate (exit on failure) & Execute Commands
  //

  for (auto& cmd : commands) {
    try {
      cmd->validate();
    } catch (const eom_app::CmdValidateException& cve) {
      std::cerr << "\n\nError Validating Command: " << cve.what() << '\n';
      return 0;
    }
  }

  for (auto& cmd : commands) {
    cmd->execute();
  }


  std::cout << "\n\n";

}

