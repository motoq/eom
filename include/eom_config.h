/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_CONFIG_H
#define EOM_CONFIG_H

#include <ostream>
#include <string>
#include <deque>

#include <cal_julian_date.h>
#include <cal_duration.h>

namespace eom_app {


/**
 * A class containing EOM scenario configuration parameters and methods
 * for parsing them from tokenized string representations.
 *
 * @author  Kurt Motekew
 * @date    20210922
 */
class EomConfig {
public:
  /**
   * Initialize with default values.  The default settings are valid,
   * but not very useful.
   */
  EomConfig() {}

  /**
   * @param  tokens  Tokenized parameters representing the start time of
   *                 the simulation
   */
  void setStartTime(std::deque<std::string>& tokens);

  /**
   * @return  The start time of the simulation
   */
  eom::JulianDate getStartTime() const noexcept { return jdStart; }

  /**
   * @return  The end time of the simulation
   */
  eom::JulianDate getStopTime() const noexcept { return jdStop; }

  /**
   * @param  tokens  Tokenized parameters indicating the duration of the
   *                 scenario
   */
  void setDuration(std::deque<std::string>& tokens);

  /**
   * @param  tokens  Tokenized parameters indicating the number of assumed leap
   *                 seconds
   */
  void setLeapSeconds(std::deque<std::string>& tokens);

  /**
   * @param  tokens  Tokenized parameters representing the time update rate
   *                 that will be used when storing reduction parameters
   */
  void setEcfEciRate(std::deque<std::string>& tokens);

  /**
   * @return  The duration between stored reduction parameters
   */
  eom::Duration getEcfEciRate() const noexcept { return dtEcfEci; }

  /**
   * @param  tokens  Change the distance conversion factor used for values
   *                 being parsed from tokenized inputs.  This factor will
   *                 convert input distance units to kilometers.
   */
  void setToKilometers(std::deque<std::string>& tokens);

  /**
   * @return  Conversion factor to kilometers for distance based values
   *          being parsed from tokenized inputs.  This factor converts
   *          input distance units to kilometers.
   */
  double getToKilometers() const noexcept { return to_km; }

  /**
   * @param  tokens  Change the time conversion factor to be used for values
   *                 being parsed from tokenized inputs.  This factor will
   *                 convert input time units to seconds.
   */
  void setToSeconds(std::deque<std::string>& tokens);

  /**
   * @return  Conversion factor to seconds for time based values being
   *          parsed from tokenized inputs.  This factor converts
   *          input times units to seconds.
   */
  double getToSeconds() const noexcept { return to_sec; }

  /**
   * @return  If an error was encountered while building the scenario,
   *          the return value is false.  Call getError().
   */
  bool isValid() const noexcept { return valid; }

  /**
   * @return  A string indicating an error state.  The string should be
   *          empty if no errors are present.
   */
  std::string getError() { return error_string; }

  /**
   * Output configuration settings to the supplied stream
   *
   * @param  stream  The stream to output the scenario setting to
   */
  void print(std::ostream& stream) const;

private:
  double to_km {1.0};
  double to_sec {1.0};
  bool valid {true};
  bool epoch_set {false};
  bool f2i_rate_set {false};
  bool leapsec_set {false};
  std::string error_string {""};
  eom::JulianDate jdStart;
  eom::JulianDate jdStop;
  eom::Duration dtEcfEci;
};


}

#endif
