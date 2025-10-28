/*
 * Copyright 2016, 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef UTL_CONST_H
#define UTL_CONST_H

/**
 * Mathematical constants such as pi or conversion factors based on standards
 * (vs. emperical "constants" such as GM).
 *
 * @author  Kurt Motekew  2016/03/14
 *                        2025/02/17  Moved redundant values from
 *                                    cal_const.h to here.
 */
namespace utl_const {

constexpr double pi   {3.141592653589793238462643383279502884};
constexpr double pio2 {pi/2.0};
constexpr double tpi  {2.0*pi};

constexpr double sec_per_day {86400.0};
constexpr double day_per_sec {1.0/sec_per_day};
constexpr double sec_per_min {60.0};
constexpr double min_per_sec {1.0/sec_per_min};
constexpr double min_per_day {1440.0};
constexpr double day_per_min {1.0/min_per_day};
constexpr double hr_per_day {24.0};
constexpr double day_per_hr {1/hr_per_day};

constexpr double rad_per_deg {pi/180.0};
constexpr double deg_per_rad {180.0/pi};
constexpr double arcsec_per_rad {180.0*3600.0/pi};       //! Per Radian
constexpr double rad_per_arcsec {1.0/arcsec_per_rad};
constexpr double mas_per_rad {1000.0*arcsec_per_rad};    //! Milliarcsecond
constexpr double rad_per_mas {1.0/mas_per_rad};

}

#endif

