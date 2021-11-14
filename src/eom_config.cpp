/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_config.h>

#include <ostream>
#include <stdexcept>
#include <string>
#include <deque>

#include <phy_const.h>
#include <eom_parse.h>
#include <cal_leap_seconds.h>

namespace eom_app {


void EomConfig::setStartTime(std::deque<std::string>& tokens)
{
  valid = false;
  if (epoch_set) {
    error_string = "Error:  Simulation Start Time Set More Than Once!";
    return;
  }
  try {
    jdStart = eom::parse_datetime(tokens);
    jdStop = jdStart;
    valid = true;
    epoch_set = true;
  } catch(std::invalid_argument& ia) {
    error_string = ia.what();
    return;
  }
}


void EomConfig::setDuration(std::deque<std::string>& tokens)
{
  valid = false;
  if (!epoch_set) {
    error_string = "Error:  Must set Simulation Start Time before Duration!";
    return;
  }
  try {
    auto dur = eom::parse_duration(tokens);
    jdStop = jdStart + dur;
    valid = true;
  } catch(std::invalid_argument& ia) {
    error_string = ia.what();
    error_string += "  EomConfig::setDuration";
  }
}


void EomConfig::setLeapSeconds(std::deque<std::string>& tokens)
{
  valid = false;
  if (tokens.size() != 1) {
    error_string = "Invalid number of parameters EomConfig::setLeapSeconds";
    return;
  }
  try {
    auto seconds = std::stod(tokens[0]);
    tokens.pop_front();
    eom::LeapSeconds& ls = eom::LeapSeconds::getInstance();
    ls.setTai_Utc(seconds);
    valid = true;
    leapsec_set = true;
  } catch(std::invalid_argument& ia) {
    error_string = ia.what();
    error_string += " EomConfig::setLeapSeconds";
  }
}


void EomConfig::setEcfEciRate(std::deque<std::string>& tokens)
{
  valid = false;
  if (tokens.size() != 2) {
    error_string = "Invalid number of parameters EomConfig::setEcfEciRate";
    return;
  }
  try {
    dtEcfEci = eom::parse_duration(tokens);
    valid = true;
    f2i_rate_set = true;
  } catch(std::invalid_argument& ia) {
    error_string = ia.what();
    error_string += "  EomConfig::setEcfEciRate";
  }
}


void EomConfig::setToKilometers(std::deque<std::string>& tokens)
{
  valid = false;
  if (tokens.size() != 1) {
    error_string = "Invalid number of parameters EomConfig::setToKilometers";
    return;
  }
  try {
    to_km = std::stod(tokens[0]);
    tokens.pop_front();
    valid = true;
  } catch(std::invalid_argument& ia) {
    error_string = ia.what();
    error_string += "  EomConfig::setToKilometers";
  }
}


void EomConfig::setToSeconds(std::deque<std::string>& tokens)
{
  valid = false;
  if (tokens.size() != 1) {
    error_string = "Invalid number of parameters EomConfig::setToSeconds";
    return;
  }
  try {
    to_sec = std::stod(tokens[0]);
    tokens.pop_front();
    valid = true;
  } catch(std::invalid_argument& ia) {
    error_string = ia.what();
    error_string += "  EomConfig::setToSeconds";
  }
}


void EomConfig::print(std::ostream& stream) const
{
  stream << "\nSimulation Start Time: " << jdStart.to_str();
  stream << "\nSimulation Stop Time:  " << jdStop.to_str();
  stream << "\nEcfEci Output Rate is " << 
             cal_const::min_per_day*dtEcfEci.getDays() << " minutes";
  eom::LeapSeconds& ls = eom::LeapSeconds::getInstance();
  stream << "\nLeap Seconds (TAI - UTC): " << ls.getTai_Utc();
  stream << "\nUsing dt eps: " <<
             phy_const::epsdt*phy_const::sec_per_tu << " seconds";
}


}
