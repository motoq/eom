/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <memory>

#include <cal_julian_date.h>
#include <astro_ephemeris_file.h>
#include <astro_sp3_ephem.h>
#include <astro_sp3_orbit.h>
#include <astro_ephemeris.h>

#include <astro_build.h>

namespace eom {

std::unique_ptr<Ephemeris>
build_ephemeris(const EphemerisFile& efd,
                const JulianDate& startTime,
                const JulianDate& stopTime,
                const std::shared_ptr<const EcfEciSys>& ecfeciSys)
{
  std::unique_ptr<Ephemeris> eph {nullptr};
  if (efd.getEphInterpMethod() == EphInterpType::hermite) {
    eph = std::make_unique<Sp3Orbit>(efd.getName(),
                                     efd.getEphFileName(),
                                     startTime,
                                     stopTime,
                                     ecfeciSys);
  } else {
    eph = std::make_unique<Sp3Ephem>(efd.getName(),
                                     efd.getEphFileName(),
                                     startTime,
                                     stopTime,
                                     ecfeciSys);
  }
  return eph;
}


}
