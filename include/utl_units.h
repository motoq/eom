/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef UTL_UNITS_H
#define UTL_UNITS_H

#include <unordered_map>

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

const std::unordered_map<std::string, double> per_tu = {
  {"TU", 1.0},
  {"Seconds", phy_const::sec_per_tu},
  {"Minutes", phy_const::min_per_tu}
};

const std::unordered_map<std::string, double> per_du = {
  {"DU", 1.0},
  {"Meters",     phy_const::m_per_du},
  {"Kilometers", phy_const::km_per_du}
};

/**
 * @param  kilometers  Distance, in kilometers
 *
 * @return  Distance in DU
 */
constexpr double operator"" _km(long double kilometers) noexcept
{
  return phy_const::du_per_km*kilometers;
}

/**
 * @param  km_sec  Velocity, in kilometers/sec
 *
 * @return  Velocity in DU/TU
 */
constexpr double operator"" _kms(long double km_sec) noexcept
{
  return phy_const::du_per_km*km_sec*phy_const::sec_per_tu;
}

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

/**
 * @param  Duration in seconds
 *
 * @return  Duration in TU
 */
constexpr double operator"" _sec(long double seconds) noexcept
{
  return phy_const::tu_per_sec*seconds;
}


}

#endif

