/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_parse.h>

#include <string>
#include <deque>
#include <stdexcept>

#include <astro_ephemeris_file.h>


namespace eom_app {

eom::EphemerisFile parse_eph_file_def(std::deque<std::string>& tokens)
{
  using namespace std::string_literals;
    // Require name, format, interpolator, and filename
  if (tokens.size() != 4) {
     throw std::invalid_argument("eom_app::parse_eph_file_def() "s +
                                 "Invalid number of tokens to parse: "s +
                                 std::to_string(tokens.size()));
  }
  auto name = tokens[0];
  tokens.pop_front();
  auto model = tokens[0];
  if (model == "SP3c") {
    tokens.pop_front();
  } else {
     throw std::invalid_argument("eom_app::parse_eph_file_def() "s +
                                 "Invalid ephemeris file type: "s +
                                 model);
  }
  auto interp = tokens[0];
  tokens.pop_front();
  eom::EphInterpType interp_type {};
  if (interp == "Chebyshev") {
    interp_type = eom::EphInterpType::chebyshev;
  } else if (interp == "Hermite") {
    interp_type = eom::EphInterpType::hermite;
  } else {
     throw std::invalid_argument("eom_app::parse_eph_file_def() "s +
                                 "Invalid ephemeris interpolation type: "s +
                                 interp);
  }

  auto file_name = tokens[0];
  tokens.pop_front();
  eom::EphemerisFile efd(name,
                         file_name,
                         eom::EphFileFormat::sp3c,
                         interp_type);
  return efd;
}


}
