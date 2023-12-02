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
#include <unordered_map>
#include <stdexcept>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_print.h>

#include <eom_command.h>
#include <eom_config.h>

namespace eom_app {


EomEphemPrinter::EomEphemPrinter(std::deque<std::string>& tokens,
                                 const EomConfig& cfg)
{
    // Read orbit name, output frame, and output filename
  using namespace std::string_literals;
  if (tokens.size() != 3) {
    throw std::invalid_argument("EomEphemPrinter::EomEphemPrinter() "s +
                                "PrintEphemeris requires 3 arguments " +
                                "vs. input " +
                                std::to_string(tokens.size()));
  }
  m_orbit_name = tokens[0];
  tokens.pop_front();
  if (!cfg.pendingOrbit(m_orbit_name)) {
    throw std::invalid_argument("EomEphemPrinter::EomEphemPrinter() "s +
                                "Invalid orbit name " + m_orbit_name);
  }
  auto frame_name = tokens[0];
  tokens.pop_front();
  if (frame_name == "GCRF") {
    m_frame = eom::EphemFrame::eci;
  } else if (frame_name == "ITRF") {
    m_frame = eom::EphemFrame::ecf;
  } else {
    throw std::invalid_argument("EomEphemPrinter::EomEphemPrinter() "s +
                                "Invalid frame type in PrintEphemeris: " +
                                frame_name);
  }
  m_file_name = tokens[0];
  tokens.pop_front();
  m_jdStart = cfg.getStartTime();
  m_jdStop = cfg.getStopTime();
}

void EomEphemPrinter::validate(const std::unordered_map<
    std::string, std::shared_ptr<eom::Ephemeris>>& ephemerides)
{
  using namespace std::string_literals;
  try {
    m_eph = ephemerides.at(m_orbit_name);
  } catch (const std::out_of_range& oor) {
    throw CmdValidateException("EomEphemPrinter::validate() "s +
                               "Invalid orbit name in PrintEphemeris: " +
                               m_orbit_name);
  }
}

// Add functionality to determine proper output rate given orbit type
// Currently defaulting to a 60 second output rate
void EomEphemPrinter::execute() const
{
  eom::print_ephemeris(m_file_name, m_jdStart, m_jdStop,
                       {60.0, phy_const::tu_per_sec}, m_frame, m_eph);
}

}
