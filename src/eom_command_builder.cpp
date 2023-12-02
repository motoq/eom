/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_command.h>

#include <memory>
#include <utility>
#include <string>
#include <deque>
#include <stdexcept>

#include <eom_config.h>
#include <eom_ephem_printer.h>
#include <eom_orbit_printer.h>
#include <eom_range_printer.h>
#include <eom_rtc_printer.h>

#include <astro_orbit_def.h>

namespace eom_app {

std::unique_ptr<EomCommand> buildCommand(std::deque<std::string>& tokens,
                                         const EomConfig& cfg)
{
  if (tokens.size() < 1) {
    throw std::invalid_argument("EomCommandBuilder::buildCommand() No tokens");
  }
  auto command_str = tokens[0];
  tokens.pop_front();
  if (command_str == "PrintEphemeris") {
    std::unique_ptr<EomCommand> command =
        std::make_unique<EomEphemPrinter>(tokens, cfg);
    return command;
  } else if (command_str == "PrintOrbit") {
    std::unique_ptr<EomCommand> command =
        std::make_unique<EomOrbitPrinter>(tokens, cfg);
    return command;
  } else if (command_str == "PrintRange") {
    std::unique_ptr<EomCommand> command =
        std::make_unique<EomRangePrinter>(tokens, cfg);
    return command;
  } else if (command_str == "PrintRangeSpectrum") {
    std::unique_ptr<EomCommand> command =
        std::make_unique<EomRangePrinter>(tokens, cfg, true);
    return command;
  } else if (command_str == "PrintRTC") {
    std::unique_ptr<EomCommand> command =
        std::make_unique<EomRtcPrinter>(tokens, cfg);
    return command;
  } else {
    throw std::invalid_argument(
        "EomCommandBuilder::buildCommand Invalid command type: " + command_str);
  }
}


}
