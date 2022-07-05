/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_eop_sys.h>

#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <stdexcept>

#include <cal_const.h>
#include <cal_julian_date.h>

namespace {
  constexpr double nsec {1.0e-9*cal_const::day_per_sec};
}

namespace eom {

EopSys::EopSys(std::string fname,
         const JulianDate& startTime, const JulianDate& stopTime)
{
    // Pad by an extra day for interpolation options
  mjd_first = static_cast<unsigned int>(startTime.getMjd() - 1.0);
  mjd_last = static_cast<unsigned int>(stopTime.getMjd() + 1.0);

    // Open IERS EOP file
  std::ifstream ifs(fname);
  if (!ifs.is_open()) {
    throw std::runtime_error("EopSys::EopSys Can't open " + fname);
  }
   // Collect column labels
  std::unordered_map<std::string, unsigned int> col_labels;
  std::string input_line;
  if (std::getline(ifs, input_line)) {
    std::stringstream header_stream(input_line);
    std::string token;
    unsigned int ndx {0};
    while (std::getline(header_stream, token, ';')) {
      col_labels[token] = ndx++;
    }
  }
    // Resolve columns of interest
  unsigned int mjd_ndx;
  unsigned int xp_ndx;
  unsigned int yp_ndx;
  unsigned int ut1mutc_ndx;
  unsigned int lod_ndx;
  unsigned int dx_ndx;
  unsigned int dy_ndx;
  try {
    mjd_ndx = col_labels.at("MJD");
    xp_ndx = col_labels.at("x_pole");
    yp_ndx = col_labels.at("y_pole");
    ut1mutc_ndx = col_labels.at("UT1-UTC");
    lod_ndx = col_labels.at("LOD");
    dx_ndx = col_labels.at("dX");
    dy_ndx = col_labels.at("dY");
  } catch (const std::out_of_range& oor) {
    throw std::runtime_error("EopSys::EopSys Bad EOP file headers");
  }

    // Force assumption of MJD being first to speed locating proper
    // in IERS EOP file (not necessary given map, but speeds things up).
  if (mjd_ndx != 0) {
    throw std::runtime_error("EopSys::EopSys MJD expected to be first column");
  }

    // Locate first line and then read data
  unsigned long mjd {0};
  std::string token;
  bool found_start {false};
  bool done {false};
  while (!done  &&  std::getline(ifs,input_line)) {
    std::stringstream eop_stream(input_line);
    if (std::getline(eop_stream, token, ';')) {
      mjd = std::stol(token);
      if (mjd == mjd_first) {
        found_start = true;
      }
      if (mjd == mjd_last) {
        done = true;
      }
    }
    if (found_start) {
      eop_record eop; 
      eop.mjd = mjd;
      std::vector<std::string> eop_tokens;
      eop_tokens.push_back("");
      while (std::getline(eop_stream, token, ';')) {
        eop_tokens.push_back(token);
      }
        // Parse available data - default to zero if "" string.
      try {
        eop.xp = std::stod(eop_tokens[xp_ndx]);
        eop.yp = std::stod(eop_tokens[yp_ndx]);
        eop.ut1mutc = std::stod(eop_tokens[ut1mutc_ndx]);
      } catch(std::invalid_argument const& ex) {
        ;
      }
      try {
        eop.lod = std::stod(eop_tokens[lod_ndx]);
      } catch(std::invalid_argument const& ex) {
        ;
      }
      try {
        eop.dx = std::stod(eop_tokens[dx_ndx]);
        eop.dy = std::stod(eop_tokens[dy_ndx]);
      } catch(std::invalid_argument const& ex) {
        ;
      }
    }
  }

  if (!found_start) {
    throw std::runtime_error("EopSys::EopSys Can't find start MJD ");
  } else if (!done) {
    throw std::runtime_error("EopSys::EopSys Can't find end MJD ");
  }

}


eop_record EopSys::getEop(const JulianDate& jd) const
{
    // No EOP data loaded
  if (eopData.size() == 0) {
    eop_record eop;
    return eop;
  }

  auto mjd = jd.getMjd();
  auto days = mjd - mjd_first;

  if (days < -nsec) {
   throw std::out_of_range("EopSys::getEop: " + std::to_string(mjd) +
                           " < " + std::to_string(mjd_first));
  }
  if (mjd_last - mjd < -nsec) {
   throw std::out_of_range("EopSys::getEop: " + std::to_string(mjd) +
                           " > " + std::to_string(mjd_last));
  }

  long ndx1 = static_cast<long>(days);
  if (days - ndx1 < nsec) {
    return eopData[ndx1]; 
  }

  //long ndx2 = ndx1 + 1L;
    // compute slope for each over 1 day
  //auto dmjd = mjd - eopData[ndx1].mjd;
  //auto dxp_dmjd = eopData[ndx2].xp - eopData[ndx1].xp;

  return eopData[ndx1]; 

}


}


