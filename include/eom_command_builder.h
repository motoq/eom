/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_COMMAND_BUILDER_H
#define EOM_COMMAND_BUILDER_H

#include <memory>
#include <string>
#include <deque>
#include <unordered_map>

#include <eom_config.h>
#include <eom_command.h>

#include <astro_ephemeris.h>

namespace eom_app {

/**
 * A utility used to build EOM commands given a tokenized string and a
 * bucket of resources.
 *
 * @author  Kurt Motekew
 * @date    20211126
 */
class EomCommandBuilder {
public:
  /**
   * Initialize with resources needed to build commands.  All input
   * collections are expected to grow.
   *
   * @param  ephemerides  Ephemeris source list
   */
  EomCommandBuilder(std::shared_ptr<
                    std::unordered_map<
                    std::string, std::shared_ptr<eom::Ephemeris>>> ephemerides);

  /**
   * Builds a command given tokens, stored resources, and the simulation
   * configuration.
   *
   * @param  tokens  List of strings describing a command that should be
   *                 created and executed
   * @param  cfg     EOM simulation configuration
   */
  std::unique_ptr<EomCommand> buildCommand(std::deque<std::string>& tokens,
                                           const EomConfig& cfg);

private:
  std::shared_ptr<std::unordered_map<std::string,
                  std::shared_ptr<eom::Ephemeris>>> m_ephemerides;
};


}

#endif
