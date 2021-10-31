/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef UTL_UNITS_H
#define UTL_UNITS_H

#include <phy_const.h>

/**
 * Utility functions allowing for conversion of external to
 * computational (internal) units of TU (time units) and DU
 * (distance units).
 *
 * @author  Kurt Motekew
 * @date    20210925
 */
namespace utl_units {

/**
 * @param  Duration in days
 *
 * @return  Duration in TU
 */
constexpr double operator"" _day(long double days) noexcept
{
  return phy_const::tu_per_day*days;
}

/**
 * @param  Duration in minutes
 *
 * @return  Duration in TU
 */
constexpr double operator"" _min(long double minutes) noexcept
{
  return phy_const::tu_per_min*minutes;
}

}

/**
 * @param  Duration in seconds
 *
 * @return  Duration in TU
 */
constexpr double operator"" _sec(long double seconds) noexcept
{
  return phy_const::tu_per_sec*seconds;
}

#endif

