/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_ephem_printer.h>

#include <memory>
#include <string>
#include <deque>
#include <stdexcept>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_print.h>

#include <eom_command.h>

namespace eom_app {


EomEphemPrinter::EomEphemPrinter(std::deque<std::string>& tokens,
    const eom::JulianDate& jdEphStart, const eom::JulianDate& jdEphStop,
    const std::shared_ptr<std::unordered_map<std::string,
                          std::shared_ptr<eom::Ephemeris>>>& ephemerides)
{
  m_ephemerides = ephemerides;

    // Read orbit name, output frame, and output filename
  using namespace std::string_literals;
  if (tokens.size() != 3) {
    throw std::invalid_argument("EomEphemPrinter::EomEphemPrinter:"s +
                                " PrintEphemeris requires 3 arguments"s +
                                " vs. input "s +
                                std::to_string(tokens.size()));
  }
  orbit_name = tokens[0];
  tokens.pop_front();
  auto frame_name = tokens[0];
  tokens.pop_front();
  if (frame_name == "GCRF") {
    frame = eom::EphemFrame::eci;
  } else if (frame_name == "ITRF") {
    frame = eom::EphemFrame::ecf;
  } else {
    throw std::invalid_argument("EomEphemPrinter::EomEphemPrinter:"s +
                                " Invalid frame type in PrintEphemeris: "s +
                                frame_name);
  }
  file_name = tokens[0];
  tokens.pop_front();
  jdStart = jdEphStart;
  jdStop = jdEphStop;
}

void EomEphemPrinter::validate()
{
  using namespace std::string_literals;
  try {
    eph = m_ephemerides->at(orbit_name);
  } catch (std::out_of_range& oor) {
    throw CmdValidateException("EomEphemPrinter::EomEphemPrinter:"s +
                               " Invalid orbit name in PrintEphemeris: "s +
                               orbit_name);
  }
}

// Add functionality to determine proper output rate given orbit type
// Currently defaulting to a 60 second output rate
void EomEphemPrinter::execute() const
{
  eom::print_ephemeris(file_name, jdStart, jdStop,
                       {60.0, phy_const::tu_per_sec}, frame, eph);
}

}
