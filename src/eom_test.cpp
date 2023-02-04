/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <string>
#include <deque>
#include <stdexcept>

#include <eom_test.h>

namespace eom_app {

void eom_test(std::deque<std::string>& tokens)
{
  if (tokens.size() != 1) {
    throw std::invalid_argument("eom_test Expected 1 token vs. " +
                                std::to_string(tokens.size()));
  }
  auto test_str = tokens[0];
  tokens.pop_front();
  if (test_str == "EarthX") {
    eom_test_earth_x();
  } else if (test_str == "SunEph") {
    eom_test_sun();
  } else {
    throw std::invalid_argument("eom_test Invalid test type: " + test_str);
  }
}


}
