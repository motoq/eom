/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_parse.h>

#include <string>
#include <deque>
#include <stdexcept>

#include <utl_units.h>
#include <cal_duration.h>

namespace eom_app {

eom::Duration parse_duration(std::deque<std::string>& tokens)
{
  using namespace std::string_literals;
  using namespace utl_units;

    // Only single value duration types at this time
  if (tokens.size() < 2) {
    throw std::invalid_argument("eom_app::parse_duration() "s +
                                "Invalid number of tokens: "s +
                                std::to_string(tokens.size()));
  }

  auto model = tokens[0];
  tokens.pop_front();
  double dur {0.0};
  try {
    dur = std::stod(tokens[0]);
    tokens.pop_front();
  } catch (const std::invalid_argument& ia) {
    throw std::invalid_argument("eom_app::parse_duration()  Invalid Duration");
  }
  if (model == "Days") {
    return {dur, 1.0_day};
  } else if (model == "Minutes") {
    return {dur, 1.0_min};
  } else if (model == "Seconds") {
    return {dur, 1.0_sec};
  } else {
    throw std::invalid_argument("eom_app::parse_duration() "s +
                                "Invalid Duration Units Type "s + model);
  }
}


}

