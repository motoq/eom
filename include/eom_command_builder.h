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
   * collections are expected to grow.  Additions to each list can be
   * made, but the order of objects must not be altered after being
   * added to any list, with the exception of ephem_nids, which maps
   * orbit names to NIDS.  Orbit and other definitions must be
   * be present before a command can be processed making use of such
   * models and resources.
   *
   * @param  ephem_nids         List of currently existing orbit Numeric IDs.
   *                            Each ID represents the indext into the
   *                            Ephemeris vector and any other resource
   *                            lists generated during EOM initialization.
   * @param  orbit_definitions  List of orbit definitions that will be
   *                            used to generate ephemeris providers.
   * @param  ephem_list         List of Ephemeris definitions      
   */
  EomCommandBuilder(const std::shared_ptr<
                          std::unordered_map<std::string, int>>& ephem_nids,
                    const std::shared_ptr<
                          std::vector<eom::OrbitDef>>& orbit_definitions,
                    const std::shared_ptr<std::vector<
                          std::shared_ptr<eom::Ephemeris>>>& ephem_list);

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
    // NIDS maps orbit names to orbit_def and ephemerides indexing
  std::shared_ptr<std::unordered_map<std::string, int>> eph_nids;
    // Order may not be changed without updating eph_nids.  Containers may
    // grow due to actions external to this class.  
  std::shared_ptr<std::vector<eom::OrbitDef>> orbit_defs;
  std::shared_ptr<std::vector<std::shared_ptr<eom::Ephemeris>>> ephemerides;
};


}

#endif
