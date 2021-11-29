/*
 * Copyright 2016, 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef CAL_CONST_H
#define CAL_CONST_H

/**
 * Constants relating to dates and standard (non-canonical) time
 * transformations.
 *
 * @author  Kurt Motekew
 */
namespace cal_const {

    // Useful Julian dates
  constexpr double j2000 {2451545.0};       ///< Jan 1, 2000 12:00:00 TT
  constexpr double j1900 {2415021.0};       ///< Jan 1, 1900 12:00:00 UT1
  constexpr double gps0  {2444244.5};       ///< Jan 6, 1980 00:00:00 UT
  constexpr double mjd   {2400000.5};       ///< Subtract from JD to get MJD
    // Time scale conversions
  constexpr double ttmtai {32.184};         ///< TT - TAI, seconds
    // Time unit conversions
  constexpr double hr_per_day  {24.0};
  constexpr double day_per_hr  {1.0/hr_per_day};
  constexpr double min_per_day {1440.0};
  constexpr double day_per_min {1.0/min_per_day};
  constexpr double sec_per_day {86400.0};
  constexpr double day_per_sec {1.0/sec_per_day};
  constexpr double sec_per_min {60.0};
  constexpr double min_per_sec {1.0/sec_per_min};

}

#endif
