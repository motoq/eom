/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include <deque>
#include <unordered_map>
#include <execution>

#include <astro_orbit_def.h>
#include <astro_rel_orbit_def.h>
#include <astro_ephemeris_file.h>
#include <axs_gp_access_def.h>

#include <eom_config.h>
#include <eom_parse.h>
#include <eom_command.h>
#include <eom_test.h>

#include <eomx_exception.h>
#include <eomx.h>

/*
 * See eomx.h
 */
void eomx_parse_input_file(const std::string& fname,
                           eom_app::EomConfig& cfg,
                           std::vector<eom::OrbitDef>& orbit_defs,
                           std::vector<eom::RelOrbitDef>& rel_orbit_defs,
                           std::vector<eom::EphemerisFile>& eph_file_defs,
                           std::unordered_map<std::string,
                               std::shared_ptr<
                                   const eom::GroundPoint>>& ground_points,
                           std::vector<eom::GpAccessDef>& gp_access_defs,
                           std::vector<
                               std::shared_ptr<eom_app::EomCommand>>& commands)
{
    // Open for input
  std::ifstream ifs(fname);
  if (!ifs.is_open()) {
    throw eom_app::EomXException("eomx::Error opening " + fname);
  }
  std::cout << "\nOpened " << fname << '\n';


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
            } else if (make == "CelestialEphemerides") {
              while (!tokens.empty()) {
                cfg.addCelestial(tokens[0]);
                tokens.pop_front();
                input_error = false;
              }
              if (input_error) {
                other_error =
                    "CelestialEphemerides:  No Celestial Bodies Listed";
              }
            } else if (make == "Orbit") {
              try {
                orbit_defs.push_back(eom_app::parse_orbit_def(tokens, cfg));
                cfg.addPendingOrbit(orbit_defs.back().getOrbitName());
                input_error = false;
              } catch (const std::invalid_argument& ia) {
                std::string xerror = ia.what();
                other_error = "Invalid Orbit definition: " + xerror;
              }
            } else if (make == "TLE") {
              if (tokens.size() > 0) {
                auto name = tokens[0];
                tokens.pop_front();
                std::string tle1;
                std::string tle2;
                bool read_two_lines {false};
                if (std::getline(ifs, tle1)) {
                  line_number++;
                  if (std::getline(ifs, tle2)) {
                    line_number++;
                    read_two_lines = true;
                  }
                }
                if (read_two_lines) {
                  orbit_defs.emplace_back(name, tle1, tle2);
                  cfg.addPendingOrbit(orbit_defs.back().getOrbitName());
                  input_error = false;
                }
              } else {
                other_error = "TLE command provided with no arguments";
              }
            } else if (make == "RelativeOrbit") {
              try {
                rel_orbit_defs.push_back(eom_app::parse_rel_orbit_def(tokens,
                                                                      cfg));
                cfg.addPendingOrbit(rel_orbit_defs.back().getOrbitName());
                input_error = false;
              } catch (const std::invalid_argument& ia) {
                std::string xerror = ia.what();
                other_error = "Invalid Relative Orbit definition: " + xerror;
              }
            } else if (make == "EphemerisFile") {
              try {
                eph_file_defs.push_back(eom_app::parse_eph_file_def(tokens));
                cfg.addPendingOrbit(eph_file_defs.back().getName());
                input_error = false;
              } catch (const std::invalid_argument& ia) {
                std::string xerror = ia.what();
                other_error = "Invalid Ephemeris File definition: " + xerror;
              }
            } else if (make == "GroundPoint") {
              try {
                eom::GroundPoint gp = eom_app::parse_ground_point(tokens, cfg);
                ground_points[gp.getName()] =
                    std::make_shared<const eom::GroundPoint>(gp);
                input_error = false;
              } catch (const std::invalid_argument& ia) {
                std::string xerror = ia.what();
                other_error = "Invalid Ground Point definition: " + xerror;
              }
            } else if (make == "SinexStations") {
              try {
                  // Insert multiple sites
                eom_app::parse_sinex_stations(tokens, ground_points);
                input_error = false;
              } catch (const std::invalid_argument& ia) {
                std::string xerror = ia.what();
                other_error = "Invalid SINEX station file format: " + xerror;
              }
            } else if (make == "Access") {
              if (tokens.size() > 0) {
                auto model = tokens[0];
                tokens.pop_front();
                if (model == "GroundPointAccess") {
                  try {
                    gp_access_defs.push_back(
                        eom_app::parse_gp_access_def(tokens, cfg));
                    input_error = false;
                  } catch (const std::invalid_argument& ia) {
                    std::string xerror = ia.what();
                    other_error =
                        "Invalid Ground Point Access definition: " + xerror;
                  }
                } else {
                  other_error = "Invalid Access command option: " + model;
                }
              } else {
                other_error = "Access command provided with no arguments";
              }
            } else if (make == "Command") {
              try {
                commands.push_back(eom_app::buildCommand(tokens, cfg));
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

    // Done with file - check status
  if (tokens.size() > 0  &&  !input_error) {
    std::cerr << "\n\n=== Warning: Reached EOF non-empty que ===";
    std::cerr << "\n        (Probably left out a ';')";
    std::cerr << "\n        (Number of Tokens " << tokens.size() << " )";
    std::cerr << "\n        (Current Token " << tokens[0] << " )\n";
  }
  if (input_error) {
    throw eom_app::EomXException(
        "eomx::Exiting on input error that can't be specified");
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
      throw eom_app::EomXException("eomx::Bad Relative Orbit Template Name: " +
                                   relOrbit.getTemplateOrbitName());
    }
  }

}
