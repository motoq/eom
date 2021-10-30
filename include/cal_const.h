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
 */
namespace cal_const {

    // Useful Julian dates
  constexpr double J2000 {2451545.0};       ///< Jan 1, 2000 12:00:00 TT
  constexpr double J1900 {2415021.0};       ///< Jan 1, 1900 12:00:00 UT1
  constexpr double GPS0  {2444244.5};       ///< Jan 6, 1980 00:00:00 UT
  constexpr double MJD   {2400000.5};       ///< Subtract from JD to get MJD
    // Time scale conversions
  constexpr double ttmtai {32.184};         ///< TT - TAI, seconds
    // Time unit conversions
  constexpr double HR_PER_DAY  {24.0};
  constexpr double DAY_PER_HR  {1.0/HR_PER_DAY};
  constexpr double MIN_PER_DAY {1440.0};
  constexpr double DAY_PER_MIN {1.0/MIN_PER_DAY};
  constexpr double SEC_PER_DAY {86400.0};
  constexpr double DAY_PER_SEC {1.0/SEC_PER_DAY};
  constexpr double SEC_PER_MIN {60.0};
  constexpr double MIN_PER_SEC {1.0/SEC_PER_MIN};

}

#endif
