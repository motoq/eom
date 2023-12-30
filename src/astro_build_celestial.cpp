/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <string>
#include <memory>
#include <vector>
#include <stdexcept>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_null_ephemeris.h>

#include <astro_build.h>

namespace eom {

std::unique_ptr<Ephemeris>
build_celestial(const std::string& name_prefix,
                const JulianDate& startTime,
                const JulianDate& stopTime,
                const std::shared_ptr<const EcfEciSys>& ecfeciSys)
{
  std::unique_ptr<Ephemeris> eph = std::make_unique<NullEphemeris>();

  return eph;
}


}
