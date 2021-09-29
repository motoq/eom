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

  void setStartTime(std::deque<std::string> tokens);

  eom::JulianDate getStartTime() const noexcept { return jdStart; }

  void setDuration(std::deque<std::string> tokens);

  eom::JulianDate getStopTime() const noexcept { return jdStop; }

  bool isValid() const noexcept { return valid; }

  void print(std::ostream& stream) const;

private:
  bool valid {true};
  std::string error_string {""};
  eom::JulianDate jdStart;
  eom::JulianDate jdStop;
};


}

#endif
