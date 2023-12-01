/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_command_builder.h>

#include <memory>
#include <utility>
#include <string>
#include <deque>
#include <unordered_map>
#include <stdexcept>

#include <eom_config.h>
#include <eom_command.h>
#include <eom_ephem_printer.h>
#include <eom_orbit_printer.h>
#include <eom_range_printer.h>
#include <eom_rtc_printer.h>

#include <astro_orbit_def.h>
#include <astro_ephemeris.h>

namespace eom_app {

EomCommandBuilder::
EomCommandBuilder(std::shared_ptr<
                  std::unordered_map<
                  std::string, std::shared_ptr<eom::Ephemeris>>> ephemerides)
{
  m_ephemerides = std::move(ephemerides);
}


std::unique_ptr<EomCommand>
EomCommandBuilder::buildCommand(std::deque<std::string>& tokens,
                                const EomConfig& cfg,
                                const std::vector<eom::OrbitDef>& orbit_defs)
{
  if (tokens.size() < 1) {
    throw std::invalid_argument("EomCommandBuilder::buildCommand() No tokens");
  }
  auto command_str = tokens[0];
  tokens.pop_front();
  if (command_str == "PrintEphemeris") {
    std::unique_ptr<EomCommand> command =
        std::make_unique<EomEphemPrinter>(tokens, cfg, orbit_defs);
    return command;
  } else if (command_str == "PrintOrbit") {
    std::unique_ptr<EomCommand> command =
        std::make_unique<EomOrbitPrinter>(tokens, cfg, orbit_defs);
    return command;
  } else if (command_str == "PrintRange") {
    std::unique_ptr<EomCommand> command =
        std::make_unique<EomRangePrinter>(tokens, cfg, orbit_defs);
    return command;
  } else if (command_str == "PrintRangeSpectrum") {
    std::unique_ptr<EomCommand> command =
        std::make_unique<EomRangePrinter>(tokens, cfg, orbit_defs, true);
    return command;
  } else if (command_str == "PrintRTC") {
    std::unique_ptr<EomCommand> command =
        std::make_unique<EomRtcPrinter>(tokens, cfg, m_ephemerides);
    return command;
  } else {
    throw std::invalid_argument(
        "EomCommandBuilder::buildCommand Invalid command type: " + command_str);
  }
}


}
