/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_ecfeci_sys.h>

#include <sofa.h>

#include <cal_julian_date.h>
#include <cal_greg_date.h>
#include <cal_duration.h>

namespace eom {

EcfEciSys::EcfEciSys(JulianDate startTime, JulianDate stopTime, Duration dt) :
                     jdStart {startTime}, jdStop {stopTime},
                                          dt_days {dt.getDays()}
{

    // Convert to TT
  double x {0.0};                 // Celestial Intermediate Pole x coordinate
  double y {0.0};                 // Celestial Intermediate Pole y coordinate
  double s {0.0};                 // CIO locator, radians
  iauXys06a(jdStart.getJdHigh(), jdStart.getJdLow(), &x, &y, &s);
}

}
