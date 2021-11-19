/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_parse.h>

#include <stdexcept>
#include <string>
#include <deque>

#include <cal_greg_date.h>

namespace eom_app {

eom::JulianDate parse_datetime(std::deque<std::string>& tokens)
{
    // Need at least the type of date/time format and the value(s)
  if (tokens.size() < 2) {
    throw std::invalid_argument(
               "Invalid number of arguments to parse_date_time");
  }
  auto model = tokens[0];
  tokens.pop_front();

    // Parse Gregorian Date input format - must have YYYY MM DD HH MM SS.S
  if (model == "GD") {
    if (tokens.size() >= 6) {
      auto year_str = tokens[0];
      tokens.pop_front();
      auto month_str = tokens[0];
      tokens.pop_front();
      auto day_str = tokens[0];
      tokens.pop_front();
      try {
        eom::GregDate gd(year_str, month_str, day_str);
        int hours {std::stoi(tokens[0])};
        tokens.pop_front();
        int minutes {std::stoi(tokens[0])};
        tokens.pop_front();
        double seconds {std::stod(tokens[0])};
        tokens.pop_front();
        eom::JulianDate jd(gd, hours, minutes, seconds);
        return jd;
      } catch(std::invalid_argument& ia) {
        throw std::invalid_argument("parse_date_time: error parsing GD values");
      }
    } else {
      throw std::invalid_argument(
                 "parse_date_time GD type requires 6 arguments");
    }
  }
  throw std::invalid_argument("Invalid parse_date_time type: " + model);
}

}

