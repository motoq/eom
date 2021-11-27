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
#include <vector>
#include <deque>
#include <unordered_map>

#include <eom_config.h>
#include <eom_command.h>

#include <astro_orbit_def.h>
#include <astro_ephemeris.h>

namespace eom_app {


/**
 * A 
 *
 * @author  Kurt Motekew
 * @date    20211126
 */
class EomCommandBuilder {
public:
  EomCommandBuilder(const std::shared_ptr<
                          std::unordered_map<std::string, int>>& orbit_ids,
                    const std::shared_ptr<
                          std::vector<eom::OrbitDef>>& orbit_definitions,
                    const std::shared_ptr<
                          std::vector<std::shared_ptr<eom::Ephemeris>>>&
                                                                  orbit_ephems);

  /**
   * Consumes tokens
   */
  std::unique_ptr<EomCommand> buildCommand(std::deque<std::string>& tokens,
                                           const EomConfig& cfg);

private:
  std::shared_ptr<std::unordered_map<std::string, int>> vids;
  std::shared_ptr<std::vector<eom::OrbitDef>> orbit_defs;
  std::shared_ptr<std::vector<std::shared_ptr<eom::Ephemeris>>> orbits;
};


}

#endif
