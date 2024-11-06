/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_parse.h>

#include <algorithm>
#include <cmath>
#include <deque>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <unordered_map>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_greg_date.h>
#include <cal_julian_date.h>
#include <astro_ground_point.h>

#include <iostream>

/*
 * @param  dts  Date time string, YY:doy:sssss, where YY is a 2 digit
 *              year, doy is the day of the year (Jan 1 = 1), and sssss
 *              is seconds into the day (0 to 86400).
 */
static eom::JulianDate get_sinex_date_time(std::string& dts);

namespace {
  constexpr double year_per_day {1.0/365.25};
  constexpr double year_per_tu {year_per_day*phy_const::day_per_tu};
    // Uninitialized position is 1 AU
  constexpr double bad_pos {phy_const::m_per_du*phy_const::du_per_au};
  constexpr double max_pos {0.1*bad_pos};
    // Uninitialized velocity is the circumference of the earth per year
  constexpr double bad_vel {utl_const::tpi*phy_const::m_per_du};
  constexpr double max_vel {0.1*bad_vel};

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
          new_rec.epoch = get_sinex_date_time(snx_tokens[epoch_ndx]);

          srec = new_rec;
        }
      } catch (const std::out_of_range& oor) {
        snx_rec new_rec;
        new_rec.code = code;
        new_rec.soln = soln;
        new_rec.epoch = get_sinex_date_time(snx_tokens[epoch_ndx]);

        station_recs[code] = new_rec;
      }
        // At this point, a record with "code" exists with soln and
        // epoch set.  Only add if soln matches current record.  Error
        // if epoch is not the same.
      auto& srec = station_recs[code];
      if (srec.soln == soln) {
        eom::JulianDate recJd = get_sinex_date_time(snx_tokens[epoch_ndx]);
        if (std::abs(srec.epoch - recJd) > phy_const::epsdt_days) {
          throw std::invalid_argument(" inconsistent epoch ");
        }
        auto component_type = snx_tokens[type_ndx];
        if (component_type == "STAX") {
          srec.x = phy_const::du_per_m*std::stod(snx_tokens[value_ndx]);
        } else if (component_type == "STAY") {
          srec.y = phy_const::du_per_m*std::stod(snx_tokens[value_ndx]);
        } else if (component_type == "STAZ") {
          srec.z = phy_const::du_per_m*std::stod(snx_tokens[value_ndx]);
        } else if (component_type == "VELX") {
          srec.dx = phy_const::du_per_m*std::stod(snx_tokens[value_ndx])*
                               year_per_tu;
        } else if (component_type == "VELY") {
          srec.dy = phy_const::du_per_m*std::stod(snx_tokens[value_ndx])*
                               year_per_tu;;
        } else if (component_type == "VELZ") {
          srec.dz = phy_const::du_per_m*std::stod(snx_tokens[value_ndx])*
                               year_per_tu;
        }
      }
    } catch (const std::invalid_argument& ia) {
      std::string estr(ia.what());
      throw std::invalid_argument("parse_slr_snx_stations(): " +
                                  estr + " " + input_line);
    }
  }
  int num_bad_recs {0};
  for (const auto& [k, v] : station_recs) {
    if (std::abs(v.x) > max_pos  ||
        std::abs(v.y) > max_pos  ||
        std::abs(v.z) > max_pos  ||
        std::abs(v.dx) > max_vel  ||
        std::abs(v.dy) > max_vel  ||
        std::abs(v.dz) > max_vel) {
      num_bad_recs++;
    }
  }
  std::cout << "\nMade " << station_recs.size() << " entries";
  std::cout << "  " << num_bad_recs << " bad entries";
}


}

static eom::JulianDate get_sinex_date_time(std::string& dts)
{
  std::string token;
  std::vector<std::string> gd_tokens;
  std::stringstream gd_stream(dts);
  while (std::getline(gd_stream, token, ':')) {
    gd_tokens.push_back(token);
  }
  if (gd_tokens.size() != 3) {
    throw std::invalid_argument("Invalid number of tokens");
  }
  int year = eom::yy_to_yyyy(std::stoi(gd_tokens[0]));
  int doy = std::stoi(gd_tokens[1]);
  double sec = std::stod(gd_tokens[2]);
  eom::GregDate gd(year, 1, 1);
  eom::JulianDate jd(gd);
  jd += sec/86400.0 - 1.0 + doy;

  return jd;
}

