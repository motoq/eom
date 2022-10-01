/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_config.h>

#include <ostream>
#include <string>
#include <deque>
#include <stdexcept>

#include <phy_const.h>
#include <utl_units.h>
#include <cal_leap_seconds.h>

#include <eom_parse.h>

/*
 * The parsing functions here should generally indicate an error
 * condition if all of the input tokens are not consumed.  This is in
 * contrast to the more generalized parsing fuctions from eom_parse.h
 * that could be a subset of a larger block of text to parse.
 */
namespace eom_app {

void EomConfig::setStartTime(std::deque<std::string>& tokens)
{
  valid = false;
  if (epoch_set) {
    error_string = "Error:  Simulation Start Time Set More Than Once!";
    return;
  }
  try {
    jdStart = parse_datetime(tokens);
    jdStop = jdStart;
    valid = true;
    epoch_set = true;
  } catch (const std::invalid_argument& ia) {
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
    auto dur = parse_duration(tokens);
    jdStop = jdStart + dur;
    valid = true;
  } catch (const std::invalid_argument& ia) {
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
  } catch (const std::invalid_argument& ia) {
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
    dtEcfEci = parse_duration(tokens);
    valid = true;
    f2i_rate_set = true;
  } catch (const std::invalid_argument& ia) {
    error_string = ia.what();
    error_string += "  EomConfig::setEcfEciRate";
  }
}


void EomConfig::setIoPerRad(std::deque<std::string>& tokens)
{
  using namespace std::string_literals;
  valid = false;
  if (tokens.size() != 1) {
    error_string = "Invalid number of parameters EomConfig::setIoPerRad";
    return;
  }
  units_angle = tokens[0];
  tokens.pop_front();
  try {
    io_per_rad = utl_units::per_rad.at(units_angle);
    valid = true;
  } catch (const std::out_of_range& oor) {
    throw std::invalid_argument("EomConfig::setIoPerRad() "s +
                                "Invalid units type: "s + units_angle);
  }
}


void EomConfig::setIoPerDu(std::deque<std::string>& tokens)
{
  using namespace std::string_literals;
  valid = false;
  if (tokens.size() != 1) {
    error_string = "Invalid number of parameters EomConfig::setIoPerDu";
    return;
  }
  units_distance = tokens[0];
  tokens.pop_front();
  try {
    io_per_du = utl_units::per_du.at(units_distance);
    valid = true;
  } catch (const std::out_of_range& oor) {
    throw std::invalid_argument("EomConfig::setIoPerDu() "s +
                                "Invalid units type: "s + units_distance);
  }
}


void EomConfig::setIoPerTu(std::deque<std::string>& tokens)
{
  using namespace std::string_literals;
  valid = false;
  if (tokens.size() != 1) {
    error_string = "Invalid number of parameters EomConfig::setIoPerTu";
    return;
  }
  units_time = tokens[0];
  tokens.pop_front();
  try {
    io_per_tu = utl_units::per_tu.at(units_time);
    valid = true;
  } catch (const std::out_of_range& oor) {
    throw std::invalid_argument("EomConfig::setIoPerTu() "s +
                                "Invalid units type: "s + units_time);
  }
}


void EomConfig::setOutputRate(std::deque<std::string>& tokens)
{
  valid = false;
  if (tokens.size() != 2) {
    error_string = "Invalid number of parameters EomConfig::setOutputRate";
    return;
  }
  try {
    dtOut = parse_duration(tokens);
    valid = true;
  } catch (const std::invalid_argument& ia) {
    error_string = ia.what();
    error_string += "  EomConfig::setOutputRate";
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
  stream << "\nOne ER is " <<
            static_cast<long>(phy_const::m_per_er) << " meters";
  stream << "\nOne DU is " <<
            static_cast<long>(phy_const::m_per_du) << " meters";
  stream << "\nOne TU is " << phy_const::sec_per_tu << " seconds";
}


}
