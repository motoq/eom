/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <memory>
#include <array>
#include <algorithm>

#include <Eigen/Dense>

#include <phy_const.h>
#include <astro_propagator_config.h>
#include <astro_orbit_def.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <astro_keplerian.h>
#include <astro_kepler.h>
#include <astro_vinti.h>
#ifdef GENPL
#include <astro_oscj2.h>
#include <astro_secj2.h>
#endif

#include <astro_build.h>

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
  } else if (pCfg.getPropagatorType() == PropagatorType::VintiJ2) {
    std::unique_ptr<Ephemeris> orbit =
         std::make_unique<Vinti>(orbitParams.getEpoch(), xeciVec, ecfeciSys,
                                                    VintiPertModel::J2_ONLY);
    return orbit;
#ifdef GENPL
  } else if (pCfg.getPropagatorType() == PropagatorType::SecJ2) {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<SecJ2>(orbitParams.getEpoch(), xeciVec, ecfeciSys);
    return orbit;
  } else if (pCfg.getPropagatorType() == PropagatorType::OscJ2) {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<OscJ2>(orbitParams.getEpoch(), xeciVec, ecfeciSys);
    return orbit;
#endif
  } else {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<Kepler>(orbitParams.getEpoch(), xeciVec, ecfeciSys);
    return orbit;
  }

}


// Relative orbit definition
std::unique_ptr<Ephemeris>
build_orbit(const RelOrbitDef& relOrbit, const OrbitDef& refOrbit,
            const std::shared_ptr<eom::Ephemeris>& refEph,
            const std::shared_ptr<const EcfEciSys>& ecfeciSys)
{
    // Only a single relative orbit definition in RelCoordType exists so
    // no decisions to make.
    // Algorithm for computing orbital elements given radial,
    // transverse, cross-track, and transverse offset distances derived
    // based on "Spacecraft Relative Orbit Geometry Description Through
    // Orbit Element Differences" by Hanspeter Schaub.  This method of
    // defining a bounding box and offset automatically guarantees the
    // energy matching constraint.
  //if (relOrbit.getRelCoordType == rtct)
  Keplerian refOe(refEph->getStateVector(refOrbit.getEpoch(),
                                         eom::EphemFrame::eci));
  std::array<double, 6> oe = refOe.getOrbitalElements();
    // Update OE with distance offsets
  double h {refOe.getAngularMomentum()};
  double r90 {h*h/phy_const::gm};
  double eta2 {1.0 - oe[1]*oe[1]};
  std::array<double, 6> dx = relOrbit.getInitialState();
  double de_r {dx[0]/oe[0]};
  double de_t {0.5*eta2*dx[1]/r90};
  double di {dx[2]/r90};
  double dm {std::sqrt(eta2)*dx[3]/r90};

  oe[1] += std::max(de_r, de_t);
  oe[2] += di;
  Keplerian newOe(oe);
  newOe.setWithMeanAnomaly(newOe.getMeanAnomaly() + dm);

  Eigen::Matrix<double, 6,1> xvec {newOe.getCartesian()};
  std::array<double, 6> xarr = {xvec(0), xvec(1), xvec(2),
                                xvec(3), xvec(4), xvec(5)};
  OrbitDef newOrbit(relOrbit.getOrbitName(),
                    refOrbit.getPropagatorConfig(),
                    refOrbit.getEpoch(),
                    xarr,
                    eom::CoordType::cartesian, eom::FrameType::gcrf);
  return build_orbit(newOrbit, ecfeciSys);
}


}

