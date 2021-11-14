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
 * A 
 *
 * @author  Kurt Motekew
 * @date    20210922
 */
class EomConfig {
public:
  /**
   * Initialize with default values.
   */
  EomConfig() {}

  void setStartTime(std::deque<std::string>& tokens);

  eom::JulianDate getStartTime() const noexcept { return jdStart; }

  eom::JulianDate getStopTime() const noexcept { return jdStop; }

  void setDuration(std::deque<std::string>& tokens);

  void setLeapSeconds(std::deque<std::string>& tokens);

  void setEcfEciRate(std::deque<std::string>& tokens);

  eom::Duration getEcfEciRate() const noexcept { return dtEcfEci; }

  void setToKilometers(std::deque<std::string>& tokens);

  void setToSeconds(std::deque<std::string>& tokens);

  bool isValid() const noexcept { return valid; }

  std::string getError() { return error_string; }

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
