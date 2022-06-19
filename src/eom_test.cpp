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
    throw std::invalid_argument("eom_test: Expected 1 token");
  }
  auto test_str = tokens[0];
  tokens.pop_front();
  if (test_str == "EarthX") {
    eom_test_earth_x();
  } else if (test_str == "GroundPoint") {
    eom_test_ground_point();
  } else if (test_str == "Cross") {
    eom_test_cross();
  } else {
    throw std::invalid_argument("Invalid test type: " + test_str);
  }
}


}
