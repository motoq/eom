/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOMX_H
#define EOMX_H

#include <string>
#include <vector>
#include <memory>
#include <unordered_map>

#include <eom_config.h>
#include <astro_orbit_def.h>
#include <astro_rel_orbit_def.h>
#include <astro_ephemeris_file.h>
#include <astro_ground_point.h>
#include <axs_gp_access_def.h>
#include <eom_command.h>

/**
 * Parses an eomx input file.  All areguments except fname are modified.
 * Generates the simulation configuration parameters along with modeling
 * component definitions (that will be used to create the actual modeling
 * components) and commands to be applied to those models.
 *
 * @param  fname           Name of file to parse
 * @param  cfg             Scenario configuration
 * @param  orbit_defs      Orbit definitions based on an initial state
 * @param  rel_orbit_defs  Orbit definitions based on another orbit
 * @param  eph_file_defs   Ephemeris file definition
 * @param  gp_access_defs  Access definition to a ground point
 * @param  commands        Actions to take on above definitions
 */
void eomx_parse_input_file(const std::string& fname,
                           eom_app::EomConfig& cfg,
                           std::vector<eom::OrbitDef>& orbit_defs,
                           std::vector<eom::RelOrbitDef>& rel_orbit_defs,
                           std::vector<eom::EphemerisFile>& eph_file_defs,
                           std::unordered_map<std::string,std::shared_ptr<
                                              eom::GroundPoint>>& ground_points,
                           std::vector<eom::GpAccessDef>& gp_access_defs,
                           std::vector<
                               std::shared_ptr<eom_app::EomCommand>>& commands);


#endif
