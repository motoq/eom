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
 * Utility functions allowing for conversion of external to internal
 * units.
 *
 * @author  Kurt Motekew
 * @date    20210925
 */
namespace utl_units {

constexpr double operator"" _day(long double days)
{
  return phy_const::TU_PER_DAY*days;
}

constexpr double operator"" _min(long double minutes)
{
  return phy_const::TU_PER_MIN*minutes;
}

}

#endif

