/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_parse.h>

#include <algorithm>
#include <deque>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <unordered_map>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ground_point.h>

#include <iostream>

namespace {
  constexpr double day_per_year {365.25};
    // Uninitialized position is 1 AU
  constexpr double bad_pos {phy_const::m_per_du*phy_const::du_per_au};
    // Uninitialized velocity is the circumference of the earth per year
  constexpr double bad_vel {utl_const::tpi*phy_const::m_per_du};

  struct snx_rec {
    double x {bad_pos};
    double y {bad_pos};
    double z {bad_pos};
    double dx {bad_vel};
    double dy {bad_vel};
    double dz {bad_vel};
    eom::JulianDate epoch;
    int soln {0};
    std::string code;
  };
}

namespace eom_app {

void parse_sinex_stations(
    std::deque<std::string>& tokens,
    std::unordered_map<std::string,
                       std::shared_ptr<eom::GroundPoint>>& ground_points)
{
  using namespace std::string_literals;
    // Need at least the filename
  if (tokens.size() < 1) {
    throw std::invalid_argument("eom_app::parse_sinex_stations() "s +
                                "1 tokens required vs. "s +
                                std::to_string(tokens.size()));
  }
  auto file_name = tokens[0];
  tokens.pop_front();
  
    // Parse time for which to compute station location
  eom::JulianDate jd;
  try {
    jd = parse_datetime(tokens);
  } catch (const std::invalid_argument& ia) {
    throw std::invalid_argument("eom_app::parse_sinex_stations() "s +
                                "invalid time for station evaluation " +
                                ia.what());
  }

  std::ifstream ifs(file_name);
  if (!ifs.is_open()) {
    throw std::invalid_argument("parse_slr_snx_stations() Can't open "
                                + file_name);
  }

  std::string input_line;
    // Locate start of station location parameters
  while (std::getline(ifs, input_line)) {
    if (input_line.find("+SOLUTION/ESTIMATE") != std::string::npos) {
        // Get header line
      std::getline(ifs, input_line);
      break;
    }
  }
    // Next line is expected to be header
  if (!ifs  ||  !(input_line.find("*INDEX") != std::string::npos)) {
    throw std::invalid_argument(
        "parse_slr_snx_stations() Missing SOLUTION/ESTIMATE header " +
        file_name);
  }
    // Collect column labels
  std::unordered_map<std::string, unsigned int> col_labels;
  std::istringstream header_stream(input_line);
  std::string token;
  unsigned int ndx {0};
  while (header_stream >> token) {
    col_labels[token] = ndx++;
  }
  std::cout << '\n' << col_labels.size();
  for (const auto& [lbl, ndx2] : col_labels) {
    std::cout << '\n' << ndx2 << ":" << lbl;
  }
    // Resolve columns of interest
  unsigned int type_ndx;                         // Pos & vel components
  unsigned int code_ndx;                         // Station designator
  unsigned int soln_ndx;                         // Solution ID
  unsigned int epoch_ndx;                        // YY:doy:sssss (sec into day)
  unsigned int unit_ndx;                         // Must be "m" or "m/y"
  unsigned int value_ndx;
  try {
    type_ndx = col_labels.at("TYPE__");
    code_ndx = col_labels.at("CODE");
    soln_ndx = col_labels.at("SOLN");
    epoch_ndx = col_labels.at("_REF_EPOCH__");
    unit_ndx = col_labels.at("UNIT");
  } catch (const std::out_of_range& oor) {
    throw std::invalid_argument(
        "parse_slr_snx_stations() Bad SNX file headers " + file_name);
  }
    // File seems to be missing '_' in estimated value header, so search
    // for the first part.  But, support more consistently formatted.
    // The value comes near the end, so should not affect placement of
    // other values
  try {
    // What the header looks like
    value_ndx = col_labels.at("__ESTIMATED");
  } catch (const std::out_of_range& oor) {
    try {
      // What the header should probably look like
      value_ndx = col_labels.at("__ESTIMATED_VALUE____");
    } catch (const std::out_of_range& oor_inner) {
    throw std::invalid_argument(
        "parse_slr_snx_stations() Bad SNX ESTIMATED header " + file_name);
    }
  }
  auto num_ndxs = std::max({type_ndx,
                            code_ndx,
                            soln_ndx,
                            epoch_ndx,
                            unit_ndx,
                            value_ndx});
    // Read solutions, compute current location, insert, until
    // end of data indicator
  std::unordered_map<std::string, snx_rec> station_recs;
  while (std::getline(ifs, input_line)) {
    // Comment line
    if (input_line.size() < 1  ||  input_line.front() == '*') {
      continue;
    }
    // End of data
    if (input_line.find("-SOLUTION/ESTIMATE") != std::string::npos) {
      break;
    }
    std::vector<std::string> snx_tokens;
    std::istringstream snx_stream(input_line);
    while (snx_stream >> token) {
      snx_tokens.push_back(token);
    }
    if (snx_tokens.size() < num_ndxs) {
      throw std::invalid_argument(
          "parse_slr_snx_stations() Bad SNX record " + input_line);
    }
      // Outer try is for parsing strings - failure terminates parsing
      // Inner try is to see if a station matching the currently parsed
      // code already exists - create a new record if it doesn't and
      // move on.
    try {
      auto code = snx_tokens[code_ndx];
      auto soln = std::stoi(snx_tokens[soln_ndx]);
      try {
        auto& srec = station_recs.at(code);
        if (srec.soln < soln) {
          snx_rec new_rec;
          new_rec.code = code;
          new_rec.soln = soln;
          // Add date
          srec = new_rec;
        }
      } catch (const std::out_of_range& oor) {
        snx_rec new_rec;
        new_rec.code = code;
        new_rec.soln = soln;
          // Add date
        station_recs[code] = new_rec;
      }
        // At this point, a record with "code" exists with soln and
        // epoch set.  Only add if soln matches current record.  Error
        // if epoch is not the same.
      auto& srec = station_recs[code];
      if (srec.soln == soln) {
          ;  // If epoch difference less than a TU
      }
    } catch (const std::invalid_argument& ia) {
      throw std::invalid_argument("parse_slr_snx_stations(): " + input_line);
    }
  }
  std::cout << "\nMade " << station_recs.size() << " entries";
}




}

/*
    double x {bad_pos};
    double y {bad_pos};
    double z {bad_pos};
    double dx {bad_vel};
    double dy {bad_vel};
    double dz {bad_vel};
    eom::JulianDate epoch;
    int soln {0};
    std::string code;
*/
