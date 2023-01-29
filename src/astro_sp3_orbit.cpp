/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_sp3_orbit.h>

#include <string>
#include <fstream>
#include <vector>
#include <memory>
#include <utility>
#include <stdexcept>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_duration.h>
#include <cal_greg_date.h>
#include <cal_julian_date.h>
#include <mth_hermite1.h>
#include <astro_ephemeris.h>
#include <mth_index_mapper.h>

namespace eom {

// Currently reads in all ephemeris within the file
Sp3Orbit::Sp3Orbit(const std::string& name,
                   const std::string& file_name,
                   const JulianDate& jdStart,
                   const JulianDate& jdStop,
                   const std::shared_ptr<const EcfEciSys>& ecfeciSys)
{
  m_name = name;
  m_ecfeciSys = ecfeciSys;

  std::ifstream ifs(file_name);
  if (!ifs.is_open()) {
    throw std::runtime_error("Sp3Orbit::Sp3Orbit() Can't open " + file_name);
  }

  std::string input_line;
  if (std::getline(ifs, input_line)) {
    if (input_line.size() < 51) {
      throw std::runtime_error(
          "Sp3Orbit::Sp3Orbit() Invalid format, line 1 too short in file " +
          file_name + " and line " + input_line);
    }
    if (input_line[2] != 'V') {
      throw std::runtime_error(
          "Sp3Orbit::Sp3Orbit() SP3 file must supply velocity; file " +
          file_name + " and line " + input_line);
    }
    auto frame = input_line.substr(48, 3);
    if (frame != "ECF") {
      throw std::runtime_error(
          "Sp3Orbit::Sp3Orbit() Only ECF frame supported; file " +
          file_name + " and line " + input_line);
    }
  }

  for (int ii=0; ii<12; ++ii) {
    if (!std::getline(ifs, input_line)) {
      throw std::runtime_error(
          "Sp3Orbit::Sp3Orbit() Incomplete header; file " + file_name);
    }
  }
  auto time_scale = input_line.substr(9, 3);
  if (time_scale != "UTC") {
    throw std::runtime_error(
        "Sp3Orbit::Sp3Orbit() Only UTC time supported; file " +
          file_name + " and line " + input_line);
  }

    // Skip to line 23 for start of ephemeris
  for (int ii=0; ii<9; ++ii) {
    if (!std::getline(ifs, input_line)) {
      throw std::runtime_error(
          "Sp3Orbit::Sp3Orbit() Ephemeris record start expected; file " +
          file_name + " and line " + input_line);
    }
  }

    // Attempt to read first ephemeris header record
  int line {1};
  JulianDate jd;
  Eigen::Matrix<double, 3, 1> pos;
  Eigen::Matrix<double, 3, 1> vel;
  std::vector<sp3_rec> sp3_records;
  std::string sp3_name {""};
  while (std::getline(ifs,input_line)) {
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
        if (tokens.size() < 7  ||  input_line[0] != '*') {
          throw std::runtime_error(
              "Sp3Orbit::Sp3Orbit() Invalid time record; file " +
              file_name + " and line " + input_line);
        }
        try {
          GregDate gd {tokens[1], tokens[2], tokens[3]};
          jd.set(gd, std::stoi(tokens[4]),
                     std::stoi(tokens[5]),
                     std::stod(tokens[6]));
        } catch (const std::invalid_argument& ia) {
          throw std::runtime_error(
              "Sp3Orbit::Sp3Orbit() Error parsing date/time values; file " +
              file_name + " and line " + input_line);
        }
        line = 2;
        break;
      case 2:
        if (tokens.size() < 4  ||  input_line[0] != 'P') {
          throw std::runtime_error(
              "Sp3Orbit::Sp3Orbit() Position record expected; file " +
              file_name + " and line " + input_line);
        }
        if (sp3_name.size() == 0) {
          sp3_name = input_line.substr(1, 3);
        } else if (sp3_name != input_line.substr(1, 3)) {
          throw std::runtime_error(
              "Sp3Orbit::Sp3Orbit() Inconsistent satellite ID; file " +
              file_name + " and line " + input_line);
        }
        try {
          pos(0) = phy_const::du_per_km*std::stod(tokens[1]);
          pos(1) = phy_const::du_per_km*std::stod(tokens[2]);
          pos(2) = phy_const::du_per_km*std::stod(tokens[3]);
        } catch (const std::invalid_argument& ia) {
          throw std::runtime_error(
              "Sp3Orbit::Sp3Orbit() Error parsing position values; file " +
              file_name + " and line " + input_line);
        }
        line = 3;
        break;
      case 3:
        if (tokens.size() < 4  ||  input_line[0] != 'V') {
          throw std::runtime_error(
              "Sp3Orbit::Sp3Orbit() Velocity record expected; file " +
              file_name + " and line " + input_line);
        }
        if (sp3_name != input_line.substr(1, 3)) {
          throw std::runtime_error(
              "Sp3Orbit::Sp3Orbit() Inconsistent satellite ID; file " +
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
              "Sp3Orbit::Sp3Orbit() Error parsing velocity values; file " +
              file_name + " and line " + input_line);
        }
        line = 1;
          // Can now insert full state vector
        sp3_records.emplace_back(jd, pos, vel);
        break;
    }
  }
  ifs.close();

  if (sp3_records.size() > 1) {
    m_jdEpoch = sp3_records.front().t;
    m_jdStart = m_jdEpoch;
    m_jdStop = sp3_records.back().t;
  } else {
    throw std::runtime_error(
        "Sp3Orbit::Sp3Orbit() Less than two eph records: " + file_name);
  }
  if (jdStart < m_jdStart) {
    throw std::runtime_error(
        "Sp3Orbit::Sp3Orbit() Ephemeris begins too late: " + file_name);
  }
  if (m_jdStop < jdStop) {
    throw std::runtime_error(
        "Sp3Orbit::Sp3Orbit() Ephemeris ends too early: " + file_name);
  }

  std::vector<std::pair<JulianDate, JulianDate>> times;
    // Generate and store Hermite interpolation objects
  for (unsigned int ii=1U; ii<sp3_records.size(); ++ii) {
    sp3_rec& r1 = sp3_records[ii-1U];
    sp3_rec& r2 = sp3_records[ii];
    double dt_tu {phy_const::tu_per_day*(r2.t - r1.t)};
    Hermite1<double, 3> hItp(dt_tu,
                             r1.p, r1.v,
                             r2.p, r2.v,
                             phy_const::epsdt);
    m_eph_interpolators.emplace_back(r1.t, r2.t, hItp);
    times.emplace_back(r1.t, r2.t);
  }
  m_ndxr = std::make_unique<IndexMapper<JulianDate>>(times);
}


Eigen::Matrix<double, 6, 1> Sp3Orbit::getStateVector(const JulianDate& jd,
                                                        EphemFrame frame) const
{
  unsigned long ndx {};
  try {
    ndx = m_ndxr->getIndex(jd);
  } catch (const std::invalid_argument& ia) {
    throw std::invalid_argument("Sp3Orbit::getStateVector() - bad time");
  }
  const auto& irec = m_eph_interpolators[ndx];
  double dt_tu {phy_const::tu_per_day*(jd - irec.jd1)};
  Eigen::Matrix<double, 6, 1> xecf;
  xecf.block<3,1>(0,0) = irec.hItp.getPosition(dt_tu);
  xecf.block<3,1>(3,0) = irec.hItp.getVelocity(dt_tu);

  if (frame == EphemFrame::eci) {
    return m_ecfeciSys->ecf2eci(jd, xecf.block<3,1>(0,0), xecf.block<3,1>(3,0));
  }

  return xecf;
}


Eigen::Matrix<double, 3, 1> Sp3Orbit::getPosition(const JulianDate& jd,
                                                  EphemFrame frame) const
{
  unsigned long ndx {};
  try {
    ndx = m_ndxr->getIndex(jd);
  } catch (const std::invalid_argument& ia) {
    throw std::invalid_argument("Sp3Orbit::getPosition - bad time");
  }
  const auto& irec = m_eph_interpolators[ndx];
  double dt_tu {phy_const::tu_per_day*(jd - irec.jd1)};
  Eigen::Matrix<double, 3, 1> xecf = irec.hItp.getPosition(dt_tu);

  if (frame == EphemFrame::eci) {
    return m_ecfeciSys->ecf2eci(jd, xecf);
  }

  return xecf;
}


}
