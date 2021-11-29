/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_COMMAND_H
#define EOM_COMMAND_H

namespace eom_app {

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
   * Carry out the functionality of the implementing class
   */
  virtual void execute() const=0;
};


}

#endif
