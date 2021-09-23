/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <ostream>
#include <stdexcept>
#include <string>
#include <deque>


#include <eom_config.h>


#include <cal_greg_date.h>
#include <cal_julian_date.h>

namespace eom {


void EomConfig::setStartTime(std::deque<std::string> tokens)
{
  valid = false;
  if (tokens.size() < 1) {
    return;
  }

  auto model = tokens[0];
  tokens.pop_front();

    // Parse Gregorian Date input format - must have at least YYYY MM DD
  if (model == "GD") {
    if (tokens.size() > 2) {
      auto year_str = tokens[0];
      tokens.pop_front();
      auto month_str = tokens[0];
      tokens.pop_front();
      auto day_str = tokens[0];
      tokens.pop_front();
      try {
        GregDate gd(year_str, month_str, day_str);
        jdStart.set(gd);
        valid = true;
      } catch(std::invalid_argument& ia) {
        error_string = ia.what();
        return;
      }
    } else {
      return;
    }
  }
  

}


void EomConfig::setDuration(std::deque<std::string> tokens)
{
  valid = false;
  if (tokens.size() < 1) {
    return;
  }

  auto model = tokens[0];
  tokens.pop_front();

    // Days
  if (model == "Days") {
    if (tokens.size() == 1) {
      try {
        auto days = std::stod(tokens[0]);
        jdStop = jdStart + days;
        valid = true;
      } catch(std::invalid_argument& ia) {
        error_string = "Invalid Days";
        return;
      }
    } else {
      return;
    }
  }
  

}



void EomConfig::print(std::ostream& stream) const
{
  stream << "\nSimulation Start Time: " << getStartTime().to_str();
  stream << "\nSimulation Stop Time:  " << getStopTime().to_str();
}

}
