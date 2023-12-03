/*
 * Copyright 2021-2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eomx.h>

#include <iostream>
#include <string>
#include <utility>
#include <memory>
#include <vector>
#include <unordered_map>
#include <execution>
#include <stdexcept>

#include <eom_config.h>
#include <eom_command.h>
#include <eomx_exception.h>
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

#include <utl_const.h>
#include <phy_const.h>
#include <astro_keplerian.h>

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
  cfg.print(std::cout);

    // Determine time span that must be supported by the simulation
    // based on the input scenario time and orbit epoch times.
  eom::JulianDate minJd = cfg.getStartTime();
  eom::JulianDate maxJd = cfg.getStopTime();
  for (const auto& orbit : orbit_defs) {
    if (orbit.getEpoch() < minJd) {
      minJd = orbit.getEpoch();
    }
    if (maxJd < orbit.getEpoch()) {
      maxJd = orbit.getEpoch();
    }
      // Backwards propagation for SP methods not currently supported
    if (orbit.getPropagatorConfig().getPropagatorType() ==
                                  eom::PropagatorType::sp) {
      if (!(orbit.getEpoch() - cfg.getStartTime()  <  phy_const::epsdt_days)) {
        std::cerr << "\n\nError:  SP orbit eopch for " <<
                      orbit.getOrbitName() <<
                     " must occur on or before the simulation start time.";
        std::cerr << "\nExiting\n";
        return 0;
      }
    }
  }

    // Ecf to Eci transformation service - immutable - pass as:
    //   const std::shared_ptr<const EcfEciSys>&
  std::shared_ptr<eom::EopSys> eopSys = nullptr;
  if (argc > 2) {
   eopSys = std::make_shared<eom::EopSys>(argv[2], minJd, maxJd);
  }
  auto f2iSys = std::make_shared<eom::EcfEciSys>(minJd,
                                                 maxJd,
                                                 cfg.getEcfEciRate(),
                                                 eopSys);

    // Ephemeris objects - build file based, then initial state based,
    // then relative orbits
  std::unordered_map<std::string,
                     std::shared_ptr<eom::Ephemeris>> ephemerides;

    // Parse interpolated ephemeris from files and
    // process sequentially
  for (const auto& ephFileDef : eph_file_defs) {
    ephemerides[ephFileDef.getName()] = 
        eom::build_ephemeris(ephFileDef, cfg.getStartTime(),
                                         cfg.getStartTime(),
                                         f2iSys);
  }

  {//==>
    // Generate orbit definitions in parallel 
  std::vector<std::unique_ptr<eom::Ephemeris>> ephvec(orbit_defs.size());
  std::transform(std::execution::par,
                 orbit_defs.begin(), orbit_defs.end(), ephvec.begin(),
                 [f2iSys](const auto& orbit) {
                   return eom::build_orbit(orbit, f2iSys);
                 }
  );
    // Move ephemerides from temporary vector to ephemeris map
  for (unsigned int ii=0; ii<ephvec.size(); ++ii) {
    auto name = ephvec[ii]->getName();
    ephemerides[name] = std::move(ephvec[ii]);
  }
  }//<==


  {//==>
    // Construct relative orbits - generate and append ephemerides.
    // Relative orbit definitions are based on primary orbit
    // definitions, not other relative orbit definitions (only
    // orbit_defs, not other rel_orbit_defs).
  std::vector<std::unique_ptr<eom::Ephemeris>> ephvec(rel_orbit_defs.size());
  std::transform(std::execution::par,
                 rel_orbit_defs.begin(), rel_orbit_defs.end(), ephvec.begin(),
                 [f2iSys, &ephemerides, &orbit_defs](const auto& relOrbit) {
        // Find reference orbit - template names already validated
      std::unique_ptr<eom::Ephemeris> eph = nullptr;
      for (const auto& templateOrbit : orbit_defs) {
        if (templateOrbit.getOrbitName() == relOrbit.getTemplateOrbitName()) {
          std::shared_ptr<eom::Ephemeris> templateEph =
                               ephemerides.at(templateOrbit.getOrbitName());
          eph = eom::build_orbit(relOrbit, templateOrbit, templateEph, f2iSys);
        }
      }
      return eph;
    }
  );
    // Move ephemerides from temporary vector to ephemeris map
  for (unsigned int ii=0; ii<ephvec.size(); ++ii) {
    auto name = ephvec[ii]->getName();
    ephemerides[name] = std::move(ephvec[ii]);
  }
  }//<==

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

