/*
 * Copyright 2021-2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <string>
#include <utility>
#include <memory>
#include <vector>
#include <unordered_map>
#include <execution>
#include <stdexcept>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_keplerian.h>
#include <astro_orbit_def.h>
#include <astro_rel_orbit_def.h>
#include <astro_ephemeris_file.h>
#include <astro_ephemeris.h>
#include <astro_eop_sys.h>
#include <astro_ecfeci_sys.h>
#include <astro_build.h>
#include <astro_ground_point.h>
#include <axs_gp_access_def.h>
#include <axs_gp_access.h>

#include <eom_config.h>
#include <eom_command.h>

#include <eomx_exception.h>
#include <eomx.h>

/**
 * Equations of Motion:  An application focused on astrodynamics related
 * problems.  This program parses an input file building models and
 * commands to be applied to those models.
 *
 * @author  Kurt Motekew
 */
int main(int argc, char* argv[])
{
    // Check for filename
  if (argc < 2  ||  argc > 3) {
    std::cerr << "\nProper use is:  " << argv[0] << " <input_file_name> or";
    std::cerr << "\n                " << argv[0] << " <input_file_name>" <<
                                                    " <eop__file_name>\n";
    return 0;
  }

    // General configuration parameter for the simulation
  eom_app::EomConfig cfg;
    // Orbit definitions, used to initialize propagators and/or generate
    // classes with buffered ephemeris
  std::vector<eom::OrbitDef> orbit_defs;
    // Orbit definitions based on other orbits.  As with orbit_defs,
    // will be used to initialize propagators and/or generate classes
    // with buffered ephemeris
  std::vector<eom::RelOrbitDef> rel_orbit_defs;
    // Ephemeris file definitions - not necessarily an orbit
  std::vector<eom::EphemerisFile> eph_file_defs;
    // Earth fixed points (ground points)
  std::unordered_map<std::string,
                     std::shared_ptr<eom::GroundPoint>> ground_points;
    // Definitions of orbit to ground access analysis requests and
    // access analysis producers
  std::vector<eom::GpAccessDef> gp_access_defs;
    // The commands populated by cmdBuilder
  std::vector<std::shared_ptr<eom_app::EomCommand>> commands;

    // Parse input file
  try {
    eomx_parse_input_file(argv[1],
                          cfg,
                          orbit_defs, rel_orbit_defs, eph_file_defs,
                          ground_points, gp_access_defs,
                          commands);
  } catch (const eom_app::EomXException& exe) {
    std::cerr << "\nError parsing input file:  " << exe.what() << '\n';
    return 0;
  }
    // And print scenario
  cfg.print(std::cout);

    // Determine time span that must be supported by the simulation
    // based on the input scenario time and orbit epoch times.
    // Then create earth fixed/inertial reference frame transformation
    // service
  std::shared_ptr<eom::EcfEciSys> f2iSys {nullptr};
  try {
    auto [minJd, maxJd] = eomx_simulation_time(cfg, orbit_defs);

    std::shared_ptr<eom::EopSys> eopSys = nullptr;
    if (argc > 2) {
     eopSys = std::make_shared<eom::EopSys>(argv[2], minJd, maxJd);
    }
    f2iSys = std::make_shared<eom::EcfEciSys>(minJd,
                                              maxJd,
                                              cfg.getEcfEciRate(),
                                              eopSys);
  } catch (const eom_app::EomXException& exe) {
    std::cerr << "\nSimulatio Time Error:  " << exe.what() << '\n';
    return 0;
  }

    // Generate ephemerides
  auto ephemerides = eomx_gen_ephemerides(cfg,
                                          orbit_defs,
                                          rel_orbit_defs,
                                          eph_file_defs,
                                          f2iSys);

    // Create access analysis objects with resources
    // Error if resource name is not available in existing containers
  std::vector<eom::GpAccess> gp_accessors;
  for (auto& axs : gp_access_defs) {
    bool first {true};
    try {
      auto gp_ptr = ground_points.at(axs.getGpName());
      first = false;
      auto eph_ptr = ephemerides.at(axs.getOrbitName());
      gp_accessors.emplace_back(cfg.getStartTime(),
                                cfg.getStopTime(),
                                *gp_ptr,
                                axs.getConstraints(),
                                eph_ptr);
    } catch (const std::out_of_range& oor) {
      if (first) {
        std::cerr << "\n\nError Assigning GP Access Ground Point: ";
      } else {
        std::cerr << "\n\nError Assigning GP Access Ephemeris: ";
      }
      std::cerr << oor.what() << '\n';
      return 0;
    }
  }

  //==>
    // Generate access times in parallel 
  std::for_each(std::execution::par,
                gp_accessors.begin(),
                gp_accessors.end(),
                [](auto& accessor) { accessor.findAllAccesses(); });
  //<==

  //
  // Print inputs
  //

    // Print derived orbits
  std::cout << '\n';
  for (const auto& relOrbit : rel_orbit_defs) {
    std::cout << "\n  " << relOrbit.getOrbitName() <<
                 "  derived from:  " <<
                 relOrbit.getTemplateOrbitName();
  }
    // Print all orbits as orbital elements
  std::cout << '\n';
  for (const auto& [name, eph] : ephemerides) {
    std::cout << "\n  " << name;
    std::cout << "\n  " << eph->getEpoch().to_str() << "    GCRF";
    eom::Keplerian oeCart(eph->getStateVector(eph->getEpoch(),
                                              eom::EphemFrame::eci));
    oeCart.print(std::cout);
  }

    // Print ground points
  for (const auto& [name, gp] : ground_points) {
    std::cout << "\n  " << name;
    gp->print(std::cout);
  }

    // Print Access Requests - print access here for now also
  for (const auto& axses : gp_accessors) {
    std::cout << "\n  Computing access for " << axses.getOrbitName() <<
                 " against " << axses.getGpName();
    for (const auto& axs : axses) {
      std::cout << '\n' << axs.rise.to_str() <<
                   "  " << axs.set.to_str() <<
                   "    {" <<
                   utl_const::deg_per_rad*std::asin(axs.sinel_rise) <<
                   ", " <<
                   utl_const::deg_per_rad*std::asin(axs.sinel_set) <<
                   "} deg Elevation";
    }
  }

  //
  // Model and command lists completed - no further modifications
  // Validate (exit on failure) & Execute Commands
  //

  for (auto& cmd : commands) {
    try {
      cmd->validate(ephemerides);
    } catch (const eom_app::CmdValidateException& cve) {
      std::cerr << "\n\nError Validating Command: " << cve.what() << '\n';
      return 0;
    }
  }

  for (auto& cmd : commands) {
    cmd->execute();
  }


  std::cout << "\n\n";

}

