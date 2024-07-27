/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_greg_date.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ephemeris_file.h>
#include <astro_sp3_chebyshev.h>
#include <astro_sp3_hermite.h>

#include <astro_build.h>

namespace eom {

std::unique_ptr<Ephemeris>
build_ephemeris(const EphemerisFile& efd,
                const JulianDate& startTime,
                const JulianDate& stopTime,
                const std::shared_ptr<const EcfEciSys>& ecfeciSys)
{
  std::vector<state_vector_rec> sp3_recs = 
      parse_sp3_file(efd.getEphFileName(), ecfeciSys->getBeginTime(),
                                           ecfeciSys->getEndTime());
  std::unique_ptr<Ephemeris> eph {nullptr};
  if (efd.getEphInterpMethod() == EphInterpType::hermite) {
    eph = std::make_unique<Sp3Hermite>(efd.getName(),
                                       sp3_recs,
                                       startTime,
                                       stopTime,
                                       ecfeciSys);
  } else {
    eph = std::make_unique<Sp3Chebyshev>(efd.getName(),
                                         sp3_recs,
                                         startTime,
                                         stopTime,
                                         ecfeciSys);
  }
  return eph;
}


std::vector<state_vector_rec> parse_sp3_file(const std::string& file_name,
                                             const JulianDate& jdStart,
                                             const JulianDate& jdStop)
{
  std::ifstream ifs(file_name);
  if (!ifs.is_open()) {
    throw std::runtime_error("parse_sp3_file() Can't open " + file_name);
  }

  std::string input_line;
  if (std::getline(ifs, input_line)) {
    if (input_line.size() < 51) {
      throw std::runtime_error(
          "parse_sp3_file() Invalid format, line 1 too short in file " +
          file_name + " and line " + input_line);
    }
    if (input_line[2] != 'V') {
      throw std::runtime_error(
          "parse_sp3_file() SP3 file must supply velocity; file " +
          file_name + " and line " + input_line);
    }
    auto frame = input_line.substr(48, 3);
    if (frame != "ECF") {
      throw std::runtime_error(
          "parse_sp3_file() Only ECF frame supported; file " +
          file_name + " and line " + input_line);
    }
  }

  for (int ii=0; ii<12; ++ii) {
    if (!std::getline(ifs, input_line)) {
      throw std::runtime_error(
          "parse_sp3_file() Incomplete header; file " + file_name);
    }
  }
  auto time_scale = input_line.substr(9, 3);
  if (time_scale != "UTC") {
    throw std::runtime_error(
        "parse_sp3_file() Only UTC time supported; file " +
          file_name + " and line " + input_line);
  }

    // Skip to line 23 for start of ephemeris
  for (int ii=0; ii<9; ++ii) {
    if (!std::getline(ifs, input_line)) {
      throw std::runtime_error(
          "parse_sp3_file() Ephemeris record start expected; file " +
          file_name + " and line " + input_line);
    }
  }

    // Attempt to read first ephemeris header record
  int line {1};
  JulianDate jd;
  Eigen::Matrix<double, 3, 1> pos;
  Eigen::Matrix<double, 3, 1> vel;
  std::vector<state_vector_rec> sp3_records;
  std::string sp3_name {""};
  bool done {false};
  while (std::getline(ifs,input_line)  &&  !done) {
    if (input_line == "EOF") {
      break;
    }
      // Skip extended state vector info - also handles bad EOF format
    if (input_line[0] == 'E') {
      continue;
    } 
    std::vector<std::string> tokens;
    std::istringstream iss(input_line);
    std::string token;
      // Tokenize the entire line
    while (iss >> token) {
      tokens.push_back(token);
    }
    switch (line) {
      case 1:
          // Ephemeris should start on line 23, but some sources
          // do not honor this - look for first '*'
        if (input_line.size() > 0  &&  input_line[0] != '*') {
          continue;
        }
        if (tokens.size() < 7  ||  input_line[0] != '*') {
          throw std::runtime_error(
              "parse_sp3_file() Invalid time record; file " +
              file_name + " and line " + input_line);
        }
        try {
          GregDate gd {tokens[1], tokens[2], tokens[3]};
          jd.set(gd, std::stoi(tokens[4]),
                     std::stoi(tokens[5]),
                     std::stod(tokens[6]));
        } catch (const std::invalid_argument& ia) {
          throw std::runtime_error(
              "parse_sp3_file() Error parsing date/time values; file " +
              file_name + " and line " + input_line);
        }
        line = 2;
        break;
      case 2:
        if (tokens.size() < 4  ||  input_line[0] != 'P') {
          throw std::runtime_error(
              "parse_sp3_file() Position record expected; file " +
              file_name + " and line " + input_line);
        }
        if (sp3_name.size() == 0) {
          sp3_name = input_line.substr(1, 3);
        } else if (sp3_name != input_line.substr(1, 3)) {
          throw std::runtime_error(
              "parse_sp3_file() Inconsistent satellite ID; file " +
              file_name + " and line " + input_line);
        }
        try {
          pos(0) = phy_const::du_per_km*std::stod(tokens[1]);
          pos(1) = phy_const::du_per_km*std::stod(tokens[2]);
          pos(2) = phy_const::du_per_km*std::stod(tokens[3]);
        } catch (const std::invalid_argument& ia) {
          throw std::runtime_error(
              "parse_sp3_file() Error parsing position values; file " +
              file_name + " and line " + input_line);
        }
        line = 3;
        break;
      case 3:
        if (tokens.size() < 4  ||  input_line[0] != 'V') {
          throw std::runtime_error(
              "parse_sp3_file() Velocity record expected; file " +
              file_name + " and line " + input_line);
        }
        if (sp3_name != input_line.substr(1, 3)) {
          throw std::runtime_error(
              "parse_sp3_file() Inconsistent satellite ID; file " +
              file_name + " and line " + input_line);
        }
        try {
          vel(0) = 1.0e-4*phy_const::sec_per_tu*
                          phy_const::du_per_km*std::stod(tokens[1]);
          vel(1) = 1.0e-4*phy_const::sec_per_tu*
                          phy_const::du_per_km*std::stod(tokens[2]);
          vel(2) = 1.0e-4*phy_const::sec_per_tu*
                          phy_const::du_per_km*std::stod(tokens[3]);
        } catch (const std::invalid_argument& ia) {
          throw std::runtime_error(
              "parse_sp3_file() Error parsing velocity values; file " +
              file_name + " and line " + input_line);
        }
        line = 1;
          // Can now insert full state vector
          // Skip if before ECF/ECI data
          // Done if beyond ECF/ECI data
        if (jd < jdStart) {
          break;
        }
        if (jd <= jdStop) {
          sp3_records.emplace_back(jd, pos, vel);
        } else {
          done = true;
        }
        break;
    }
  }

  return sp3_records;
}


}
