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
#include <utility>
#include <memory>
#include <unordered_map>

#include <cal_julian_date.h>
#include <astro_orbit_def.h>
#include <astro_rel_orbit_def.h>
#include <astro_ephemeris_file.h>
#include <astro_ecfeci_sys.h>
#include <astro_ground_point.h>
#include <axs_gp_access_def.h>
#include <axs_gp_access.h>

#include <eom_config.h>
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
 *
 * @throws  eom_app::EomXException if there is an error parsing the file
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

/**
 * Determine the time span that must be supported by the simulation
 * resources based on input scenario time and orbit epoch times.
 *
 * param  cfg             Scenario configuration
 * @param  orbit_defs      Orbit definitions based on an initial state
 *
 * @return  [Minimum, Maximum] required Julian dates to support the
 *          simulation.
 * @throws  eom_app::EomXException if the simulation time can't be
 *          supported
 */
std::pair<eom::JulianDate, eom::JulianDate>
eomx_simulation_time(const eom_app::EomConfig& cfg,
                     const std::vector<eom::OrbitDef>& orbit_defs);

/**
 * Generate orbital ephemeris objects.  This is done either by reading
 * in an ephemeris file to be interpolated, initializing a general
 * perturbations method, or propagating and storing ephemeris via
 * special perturbation methods.
 *
 * @param  cfg             Scenario configuration
 * @param  orbit_defs      Orbit definitions based on an initial state
 * @param  rel_orbit_defs  Orbit definitions based on another orbit
 * @param  eph_file_defs   Ephemeris file definition
 * @param  f2iSys          ECF/ECI transformation service that will be
 *                         copied into each ephemeris type created.
 *
 * @return  Map of ephemerides indexed by orbit name
 */
std::unordered_map<std::string,std::shared_ptr<eom::Ephemeris>>
eomx_gen_ephemerides(const eom_app::EomConfig& cfg,
                     const std::vector<eom::OrbitDef>& orbit_defs,
                     const std::vector<eom::RelOrbitDef>& rel_orbit_defs,
                     const std::vector<eom::EphemerisFile>& eph_file_defs,
                     const std::shared_ptr<eom::EcfEciSys>& f2iSys);

/**
 * Given access analysis definitions, assign resources and run analysis.
 *
 * @param  cfg             Scenario configuration
 * @param  ground_points   Available ground point definitions
 * @param  ephemerides     Available ephemeris resources
 * @param  gp_access_defs  Access analysis requests to process
 *
 * @throws  eom_app::EomXException if necessary resources are not
 *          available
 */
std::vector<eom::GpAccess> eomx_gen_gp_accesses(
    const eom_app::EomConfig& cfg,
    const std::unordered_map<std::string,
                             std::shared_ptr<eom::GroundPoint>>& ground_points,
    const std::unordered_map<std::string,
                             std::shared_ptr<eom::Ephemeris>>& ephemerides,
    const std::vector<eom::GpAccessDef>& gp_access_defs);

#endif
