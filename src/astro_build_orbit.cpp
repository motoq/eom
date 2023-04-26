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
#include <cal_julian_date.h>
#include <astro_propagator_config.h>
#include <astro_orbit_def.h>
#include <astro_ecfeci_sys.h>
#include <astro_deq.h>
#include <mth_ode_solver.h>
#include <astro_rk4.h>
#include <astro_adams_4th.h>
#include <astro_gravity.h>
#include <astro_gravity_jn.h>
#include <astro_gravity_std.h>
#include <astro_force_model.h>
#include <astro_keplerian.h>
#include <astro_ephemeris.h>
#include <astro_sp_ephemeris.h>
#include <astro_kepler.h>
#include <astro_kepler_prop.h>
#include <astro_vinti.h>
#include <astro_vinti_prop.h>
#include <astro_sun_meeus.h>
#include <astro_moon_meeus.h>
#include <astro_third_body_gravity.h>
#ifdef GENPL
#include <astro_gravt.h>
#include <astro_oscj2.h>
#include <astro_secj2.h>
#include <astro_gauss_jackson.h>
#include <astro_gj_lite.h>
#endif

#include <astro_build.h>

namespace eom {

std::unique_ptr<Ephemeris> 
build_orbit(const OrbitDef& orbitParams,
            const std::shared_ptr<const EcfEciSys>& ecfeciSys)
{
  std::array<double, 6> xeci_array = orbitParams.getInitialState();
  Eigen::Matrix<double, 6, 1> xeciVec;
    // Once more reference frame besides GCRF is added for Keplerian
    // elements, the transformation to GCRF should be done here.
  if (orbitParams.getCoordinateType() == CoordType::keplerian) {
    Keplerian kep(xeci_array);
    xeciVec = kep.getCartesian();
  } else {
    for (int ii=0; ii<6; ++ii) {
      xeciVec(ii) = xeci_array[ii];
    }
  }
    // ITRF to GCRF for Cartesian
  if (orbitParams.getCoordinateType() == CoordType::cartesian  &&
      orbitParams.getReferenceFrameType() == FrameType::itrf) {
    xeciVec = ecfeciSys->ecf2eci(orbitParams.getEpoch(),
                                 xeciVec.block<3, 1>(0, 0),
                                 xeciVec.block<3, 1>(3, 0));
  }

    // Build orbit definition based on propagator configuration
    // Default options are two-body systems if PropagatorConfig
    // values fall out if sync with options checked here.
  PropagatorConfig pCfg = orbitParams.getPropagatorConfig();
  if (pCfg.getPropagatorType() == PropagatorType::sp) {
      // Force model must always include central body
    std::unique_ptr<Gravity> forceModel {nullptr};
    if (pCfg.getGravityModel() == GravityModel::jn) {
      forceModel = std::make_unique<GravityJn>(pCfg.getDegree());
    } else if (pCfg.getGravityModel() == GravityModel::std) {
      forceModel = std::make_unique<GravityStd>(pCfg.getDegree(),
                                                pCfg.getOrder());
#ifdef GENPL
    } else if (pCfg.getGravityModel() == GravityModel::gravt) {
      forceModel = std::make_unique<Gravt>(pCfg.getDegree(), pCfg.getOrder());
#endif
    } else {
      forceModel = std::make_unique<GravityJn>(0);
    }
    auto deq = std::make_unique<Deq>(std::move(forceModel), ecfeciSys);
      // Additional force models
    if (pCfg.getSunGravityModel() == SunGravityModel::meeus) {
      std::unique_ptr<Ephemeris> sunEph = std::make_unique<SunMeeus>(ecfeciSys);
      std::unique_ptr<ForceModel> sunGrav =
              std::make_unique<ThirdBodyGravity>(phy_const::gm_sun,
                                                 std::move(sunEph));
      deq->addForceModel(std::move(sunGrav));
    }
    if (pCfg.getMoonGravityModel() == MoonGravityModel::meeus) {
      std::unique_ptr<Ephemeris> moonEph =
              std::make_unique<MoonMeeus>(ecfeciSys);
      std::unique_ptr<ForceModel> moonGrav =
              std::make_unique<ThirdBodyGravity>(phy_const::gm_moon,
                                                 std::move(moonEph));
      deq->addForceModel(std::move(moonGrav));
    }
      // Integrator
    std::unique_ptr<OdeSolver<JulianDate, double, 6>> sp {nullptr};
    if (pCfg.getPropagator() == Propagator::rk4) {
      sp = std::make_unique<Rk4>(std::move(deq),
                                 pCfg.getStepSize(),
                                 orbitParams.getEpoch(),
                                 xeciVec);
    } else if (pCfg.getPropagator() == Propagator::adams4) {
      sp = std::make_unique<Adams4th>(std::move(deq),
                                      pCfg.getStepSize(),
                                      orbitParams.getEpoch(),
                                      xeciVec);
#ifdef GENPL
    } else if (pCfg.getPropagator() == Propagator::gj) {
      sp = std::make_unique<GaussJackson>(std::move(deq),
                                          orbitParams.getEpoch(),
                                          xeciVec);
    } else if (pCfg.getPropagator() == Propagator::gjs) {
      sp = std::make_unique<GjLite>(std::move(deq),
                                    orbitParams.getEpoch(),
                                    xeciVec);
#endif
    } else {
      sp = std::make_unique<Rk4>(std::move(deq),
                                 pCfg.getStepSize(),
                                 orbitParams.getEpoch(),
                                 xeciVec);
    }
      // Ready to generate ephemeris
    std::unique_ptr<Ephemeris> orbit =
        std::make_unique<SpEphemeris>(orbitParams.getOrbitName(),
                                      pCfg.getStartTime(),
                                      pCfg.getStopTime(),
                                      ecfeciSys,
                                      std::move(sp));
    return orbit;
  } else if (pCfg.getPropagatorType() == PropagatorType::kepler1) {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<Kepler>(orbitParams.getOrbitName(),
                                    orbitParams.getEpoch(),
                                    xeciVec,
                                    ecfeciSys);
    return orbit;
  } else if (pCfg.getPropagatorType() == PropagatorType::kepler1mod) {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<KeplerProp>(orbitParams.getOrbitName(),
                                        orbitParams.getEpoch(),
                                        xeciVec,
                                        ecfeciSys);
    return orbit;
  } else if (pCfg.getPropagatorType() == PropagatorType::vinti6) {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<Vinti>(orbitParams.getOrbitName(),
                                   orbitParams.getEpoch(),
                                   xeciVec,
                                   ecfeciSys);
    return orbit;
  } else if (pCfg.getPropagatorType() == PropagatorType::vinti_j2) {
    std::unique_ptr<Ephemeris> orbit =
         std::make_unique<Vinti>(orbitParams.getOrbitName(),
                                 orbitParams.getEpoch(),
                                 xeciVec,
                                 ecfeciSys,
                                 VintiPertModel::J2_ONLY);
    return orbit;
  } else if (pCfg.getPropagatorType() == PropagatorType::vinti6mod) {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<VintiProp>(orbitParams.getOrbitName(),
                                       orbitParams.getEpoch(),
                                       xeciVec,
                                       ecfeciSys);
    return orbit;
#ifdef GENPL
  } else if (pCfg.getPropagatorType() == PropagatorType::sec_j2) {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<SecJ2>(orbitParams.getOrbitName(),
                                   orbitParams.getEpoch(),
                                   xeciVec,
                                   ecfeciSys);
    return orbit;
  } else if (pCfg.getPropagatorType() == PropagatorType::osc_j2) {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<OscJ2>(orbitParams.getOrbitName(),
                                   orbitParams.getEpoch(),
                                   xeciVec,
                                   ecfeciSys);
    return orbit;
#endif
  } else {
    std::unique_ptr<Ephemeris> orbit =
           std::make_unique<Kepler>(orbitParams.getOrbitName(),
                                    orbitParams.getEpoch(),
                                    xeciVec,
                                    ecfeciSys);
    return orbit;
  }

}


// Relative orbit definition
std::unique_ptr<Ephemeris>
build_orbit(const RelOrbitDef& relOrbit,
            const OrbitDef& refOrbit,
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
  Keplerian refOe(refEph->getStateVector(refEph->getEpoch(),
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
                    refEph->getEpoch(),
                    xarr,
                    eom::CoordType::cartesian, eom::FrameType::gcrf);
  return build_orbit(newOrbit, ecfeciSys);
}


}

