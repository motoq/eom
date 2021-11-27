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

//enum class EphemFrame {
//  eci,                            ///< GCRF (IAU 2000A/2006)
//  ecf                             ///< ITRF (~WGS 84)
//};

class EomCommand {
public:
  virtual ~EomCommand() {}

  virtual void execute() const=0;
};


}

#endif
