/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_PRINT_H
#define ASTRO_PRINT_H

#include <memory>
#include <string>

#include <cal_duration.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>

namespace eom {

void print_ephemeris(std::string file_name, 
                     const JulianDate& jdStart, const JulianDate& jdStop, 
                     const Duration& dtout, EphemFrame frame,
                     const std::shared_ptr<const Ephemeris>& orbit);

}

#endif
