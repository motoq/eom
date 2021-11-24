/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_build.h>

#include <memory>
#include <array>

#include <Eigen/Dense>

#include <astro_propagator_config.h>
#include <astro_orbit_def.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <astro_kepler.h>
#include <astro_vinti.h>

namespace eom {

std::unique_ptr<Ephemeris> 
build_orbit(const OrbitDef& orbitParams,
            const std::shared_ptr<const EcfEciSys>& ecfeciSys)
{
    // Once more than one CoordType and FrameType are added
    // Convert xeciVec to Cartesian GCRF here
    // For now, just copy elements into vector
  std::array<double, 6> xeci_array = orbitParams.getInitialState();
  Eigen::Matrix<double, 6, 1> xeciVec;
  for (int ii=0; ii<6; ++ii) {
    xeciVec(ii) = xeci_array[ii];
  }

  PropagatorConfig pCfg = orbitParams.getPropagatorConfig();
  if (pCfg.getPropagatorType() == PropagatorType::Kepler1) {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<Kepler>(orbitParams.getEpoch(), xeciVec, ecfeciSys);
    return orbit;
  } else if (pCfg.getPropagatorType() == PropagatorType::Vinti6) {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<Vinti>(orbitParams.getEpoch(), xeciVec, ecfeciSys);
    return orbit;
  } else {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<Kepler>(orbitParams.getEpoch(), xeciVec, ecfeciSys);
    return orbit;
  }

}

}

