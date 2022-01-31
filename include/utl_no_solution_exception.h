/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef UTL_NO_SOLUTION_EXCEPTION_H
#define UTL_NO_SOLUTION_EXCEPTION_H

#include <stdexcept>
#include <string>

namespace eom {

/**
 * Exception to be thrown when an algorithm encounters a condition under
 * which a solution does not exist
 */
class NoSolutionException : public std::runtime_error {
public:
  NoSolutionException() : std::runtime_error("NoSolutionException") { }

  NoSolutionException(const std::string& msg) : std::runtime_error(msg) { }
};

}

#endif
