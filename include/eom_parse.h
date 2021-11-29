/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_PARSE_H
#define EOM_PARSE_H

#include <array>
#include <string>
#include <deque>

#include <cal_julian_date.h>
#include <cal_duration.h>
#include <astro_orbit_def.h>

#include <eom_config.h>

namespace eom_app {

/**
 * Parses a list of string tokens representing a date and time.  The
 * first token indicates which type.  Gregorian date with time of day is
 * the only available input format at this time:
 *   GD YYYY MM DD HH MM SS.S
 *
 * @param  List of tokens to be parsed.  This list is modified such that
 *         all parsed values are consumed (pop_front()).
 *
 * @return  A point in time defined by a Julian date.
 *
 * @throws  An invalid_argument exception if parsing fails.  No error is
 *          thrown if the list of tokens is not empty upon completion.
 */
eom::JulianDate parse_datetime(std::deque<std::string>& tokens);

/**
 * Parses two tokens to form a duration in time.
 *
 * @param  Two tokens to parse, the first a string representing the units of
 *         time, and the second a string that will be converted to a double.
 *         This list is modified such that all parsed values are consumed
 *         (pop_front()).
 *
 * @return  A duration
 *
 * @throws  An invalid_argument exception if parsing fails.  No error is
 *          thrown if the list of tokens is not empty upon completion.
 */
eom::Duration parse_duration(std::deque<std::string>& tokens);

/**
 * Parses an orbit definition
 *
 * @param  Tokens consisting of an orbit name, type, epoch, and state
 *         vector.  This list is modified such that all parsed values
 *         are consumedi (pop_front()).
 *
 * @return  An orbit definition, used in the generation of an orbit
 *          model
 *
 * @throws  An invalid_argument exception if parsing fails.  No error is
 *          thrown if the list of tokens is not empty upon completion.
 */
eom::OrbitDef parse_orbit_def(std::deque<std::string>& tokens,
                              const EomConfig& cfg);

/**
 * Parses an orbit state vector
 *
 * @param  Tokens consisting of coordinat system time, reference frame,
 *         and 6 components of the state vector.
 *
 * @return  Orbit state vector.
 *
 * @throws  An invalid_argument exception if parsing fails.  No error is
 *          thrown if the list of tokens is not empty upon completion.
 */
std::array<double, 6> parse_state_vector(std::deque<std::string>& tokens,
                                         const EomConfig& cfg);


}

#endif
