/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_range_printer.h>

#include <memory>
#include <string>
#include <vector>
#include <deque>
#include <stdexcept>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_print.h>

#include <eom_parse.h>

namespace eom_app {


EomRangePrinter::EomRangePrinter(std::deque<std::string>& tokens,
       const eom::JulianDate& jdEphStart, const eom::JulianDate& jdEphStop,
       const std::shared_ptr<std::unordered_map<std::string, int>>& ephem_ndxs,
       const std::shared_ptr<std::vector<std::shared_ptr<eom::Ephemeris>>>&
                                                                  ephem_list)
{
    // Read orbit name, output frame, and output filename
  using namespace std::string_literals;
  if (tokens.size() != 5) {
    throw std::invalid_argument("EomRangePrinter::EomRangePrinter:"s +
                                " PrintRange requires 4 arguments"s +
                                " vs. input "s +
                                std::to_string(tokens.size()));
  }
  for (int ii=0; ii<2; ++ii) {
    auto orbit_name = tokens[0];
    tokens.pop_front();
    try {
      endxs[ii] = ephem_ndxs->at(orbit_name);
    } catch (std::out_of_range& oor) {
      throw std::invalid_argument("EomRangePrinter::EomRangePrinter:"s +
                                  " Invalid orbit name in PrintRange: "s +
                                    orbit_name);
    }
  }
  dt = parse_duration(tokens);
  file_name = tokens[0];
  tokens.pop_front();
  jdStart = jdEphStart;
  jdStop = jdEphStop;
  ephemerides = ephem_list;
}

// Add functionality to determine proper output rate given orbit type
// Currently defaulting to a 60 second output rate
void EomRangePrinter::execute() const
{
  auto jd = jdStart;
  while (jd < jdStop) {
   ;
    jd += dt;
  }
}

}
