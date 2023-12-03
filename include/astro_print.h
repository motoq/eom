/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_PRINT_H
#define ASTRO_PRINT_H

#include <string>

#include <cal_duration.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>

namespace eom {

/**
 * Print utility function declarations
 */

/**
 * Write .e format ephemeris to a file.  Position and velocity are
 * output in meters and meters/sec
 * 
 * @param  file_name  Name of the output file to be created.  The
 *                    full filename should be used including the desired
 *                    extension as .e is not assumed.
 * @param  jdStart    Time of first ephemeris output
 * @param  jdStop     Time of final ephemeris output
 * @param  dtout      Output rate
 * @param  frame      Output reference frame
 * @param  orbit      Ephemeris source
 */
void print_ephemeris(std::string file_name, 
                     const JulianDate& jdStart, const JulianDate& jdStop, 
                     const Duration& dtout, EphemFrame frame,
                     const Ephemeris& orbit);

}

#endif
