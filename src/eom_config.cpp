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
#include <array>
#include <deque>


#include <eom_config.h>


#include <utl_units.h>
#include <cal_greg_date.h>
#include <cal_duration.h>

static eom::Duration parse_duration(std::deque<std::string>&);

namespace eom_app {


void EomConfig::setStartTime(std::deque<std::string>& tokens)
{
  valid = false;
  if (epoch_set) {
    error_string = "Error:  Simulation Start Time Set More Than Once!";
    return;
  }
  if (tokens.size() < 1) {
    error_string = "Invalid Number of Arguments for Start Time";
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
      eom::GregDate gd;
      try {
        gd.set(year_str, month_str, day_str);
        int hours {0};
        int minutes {0};
        double seconds {0.0};
        if (tokens.size() > 0) {
          hours = std::stoi(tokens[0]);
          tokens.pop_front();
        }
        if (tokens.size() > 0) {
          minutes = std::stoi(tokens[0]);
          tokens.pop_front();
        }
        if (tokens.size() > 0) {
          seconds = std::stod(tokens[0]);
          tokens.pop_front();
        }
        if (tokens.size() == 0) {
          jdStart.set(gd, hours, minutes, seconds);
          valid = true;
          epoch_set = true;
        }
      } catch(std::invalid_argument& ia) {
        error_string = ia.what();
        return;
      }
    } else {
      return;
    }
  }
  

}


void EomConfig::setDuration(std::deque<std::string>& tokens)
{
  valid = false;
  try {
    auto dur = parse_duration(tokens);
    jdStop = jdStart + dur;
    valid = true;
  } catch(std::invalid_argument& ia) {
    error_string = ia.what();
    error_string += "  EomConfig::setDuration";
  }
}


void EomConfig::setEcfEciRate(std::deque<std::string>& tokens)
{
  valid = false;
  try {
    dtEcfEci = parse_duration(tokens);
    valid = true;
    f2i_rate_set = true;
  } catch(std::invalid_argument& ia) {
    error_string = ia.what();
    error_string += "  EomConfig::setEcfEciRate";
  }
}


void EomConfig::print(std::ostream& stream) const
{
  stream << "\nSimulation Start Time: " << jdStart.to_str();
  stream << "\nSimulation Stop Time:  " << jdStop.to_str();
  stream << "\nEcfEci Output Rate is " << 
            cal_const::MIN_PER_DAY*dtEcfEci.getDays() << " minutes";
}


}

/**
 * Parses two tokens, the first a string representing the units of time,
 * and the second a string that will be converted to a double.  If there
 * are not exactly two arguments, or if a double can't be parsed from
 * the 2nd entry, an exception is thrown.  The input list is modified -
 * it should be empty upon return from this function.
 */
static eom::Duration parse_duration(std::deque<std::string>& tokens)
{
  using namespace utl_units;

  if (tokens.size() != 2) {
    throw std::invalid_argument("Invalid number of arguments to Duration");
  }

  auto model = tokens[0];
  tokens.pop_front();

    // All single value duration types
  double dur {0.0};
  try {
    dur = std::stod(tokens[0]);
  } catch(std::invalid_argument& ia) {
    throw std::invalid_argument("Invalid Duration");
  }
  if (model == "Days") {
    return {dur, 1.0_day};
  } else if (model == "Minutes") {
    return {dur, 1.0_min};
  } else if (model == "Seconds") {
    return {dur, 1.0_sec};
  } else {
    throw std::invalid_argument("Invalid Duration Units Type");
  }
    
}
