/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOMX_EXCEPTION_H
#define EOMX_EXCEPTION_H

#include <stdexcept>
#include <string>

namespace eom_app {

/**
 * Exception to be thrown when there is an exit() worthy eomx failure
 */
class EomXException : public std::runtime_error {
public:
  EomXException() : std::runtime_error("EomXException") { }

  EomXException(const std::string& msg) : std::runtime_error(msg) { }
};

}

#endif
