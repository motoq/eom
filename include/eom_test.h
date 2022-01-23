/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_TEST_H
#define EOM_TEST_H

#include <string>
#include <deque>

namespace eom_app {

/**
 * Location for validity checking or intermediate algorithm testing
 */

/**
 * Primary entry into test functions - reads type and calls the
 * appropriate function
 *
 * @param  tokens  A single value indicating which self contained test
 *                 to run.  The token is consumed.
 *
 * @throws  An invalid_argument exception if parsing fails.  An error is
 *          thrown if the token is not valid, or if more than one is
 *          supplied.
 */
void eom_test(std::deque<std::string>& tokens);

/**
 * Runs validity tests for the eom::GroundPoint type.
 */
void eom_test_ground_point();


}

#endif
