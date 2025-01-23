/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <string>
#include <utility>
#include <memory>
#include <vector>
#include <unordered_map>
#include <execution>

#include <astro_orbit_def.h>
#include <astro_rel_orbit_def.h>
#include <astro_ephemeris_file.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <astro_build.h>

#include <eomx.h>

/**
 * See eomx.h
 */
std::unordered_map<std::string,std::shared_ptr<eom::Ephemeris>>
eomx_gen_ephemerides(const eom_app::EomConfig& cfg,
                     const std::vector<eom::OrbitDef>& orbit_defs,
                     const std::vector<eom::RelOrbitDef>& rel_orbit_defs,
                     const std::vector<eom::EphemerisFile>& eph_file_defs,
                     const std::shared_ptr<const eom::EcfEciSys>& f2iSys)
{
    // Celestial Ephemeris objects - read ephemerides from files
  std::unordered_map<std::string,
                     std::vector<eom::state_vector_rec>> celestials;
  std::vector<std::string> celestial_names = cfg.getCelestials();
  for (const auto& name : celestial_names) {
    celestials[name] = eom::build_celestial(name, cfg.getStartTime(),
                                                  cfg.getStopTime());
  }

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
                 [f2iSys, &celestials](const auto& orbit) {
                   return eom::build_orbit(orbit, f2iSys, celestials);
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
                 [f2iSys,
                  &ephemerides,
                  &orbit_defs,
                  &celestials](const auto& relOrbit) {
        // Find reference orbit - template names already validated
      std::unique_ptr<eom::Ephemeris> eph = nullptr;
      for (const auto& templateOrbit : orbit_defs) {
        if (templateOrbit.getOrbitName() == relOrbit.getTemplateOrbitName()) {
          std::shared_ptr<eom::Ephemeris> templateEph =
                               ephemerides.at(templateOrbit.getOrbitName());
          eph = eom::build_orbit(relOrbit,
                                 templateOrbit,
                                 *templateEph,
                                 f2iSys, celestials);
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


  return ephemerides;
}

