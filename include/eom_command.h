/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_COMMAND_H
#define EOM_COMMAND_H

#include <stdexcept>
#include <string>

namespace eom_app {

/**
 * Exception to be thrown when the validate command fails
 */
class CmdValidateException : public std::runtime_error {
public:
  CmdValidateException() : std::runtime_error("CmdValidateException") { }

  CmdValidateException(const std::string& msg) : std::runtime_error(msg) { }
};

/**
 * Interface for an EOM application command function.
 *
 * @author  Kurt Motekew
 * @date    202111xx
 */
class EomCommand {
public:
  virtual ~EomCommand() {}

  /**
   * Perform final validation of command before potentially
   * computationally intensive processing takes place.
   */
  virtual void validate() = 0;

  /**
   * Carry out the functionality of the implementing class
   */
  virtual void execute() const = 0;
};


}

#endif
