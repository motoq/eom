/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef UTL_NONCONVERGENCE_EXCEPTION_H
#define UTL_NONCONVERGENCE_EXCEPTION_H

#include <stdexcept>
#include <string>

namespace eom {

/**
 * Exception to be thrown when an algorithm fails to converge
 */
class NonconvergenceException : public std::runtime_error {
public:
  NonconvergenceException() : std::runtime_error("NonconvergenceException") { }

  NonconvergenceException(const std::string& msg) : std::runtime_error(msg) { }
};

}

#endif
