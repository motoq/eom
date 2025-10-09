/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */


#include <algorithm>
#include <array>
#include <memory>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <mth_ode_solver.h>
#include <astro_adams_4th.h>
#include <astro_deq.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_fandg.h>
#include <astro_force_model.h>
#include <astro_gravity.h>
#include <astro_gravity_jn.h>
#include <astro_gravity_std.h>
#include <astro_hermite1_eph.h>
#include <astro_hermite1_tc_eph.h>
#include <astro_kepler.h>
#include <astro_keplerian.h>
#include <astro_kepler_prop.h>
#include <astro_moon_meeus.h>
#include <astro_orbit_def.h>
#include <astro_propagator_config.h>
#include <astro_rk4.h>
#include <astro_rk4s.h>
#include <astro_secular_j2.h>
#include <astro_gpx.h>
#include <astro_sgp4.h>
#include <astro_sp_ephemeris.h>
#include <astro_srp_spherical.h>
#include <astro_sun_meeus.h>
#include <astro_third_body_gravity.h>
#include <astro_vinti.h>
#include <astro_vinti_prop.h>
#ifdef GENPL
#include <astro_gauss_jackson.h>
#include <astro_gj_lite.h>
#include <astro_gravt.h>
#include <astro_oscj2.h>
#include <astro_secj2.h>
#endif

#include <astro_build.h>

namespace eom {

std::unique_ptr<Ephemeris> 
build_orbit(const OrbitDef& orbitParams,
            const std::shared_ptr<const EcfEciSys>& ecfeciSys,
            const std::unordered_map<std::string,
                                     std::vector<eom::state_vector_rec>>& ceph)
{
    // Use of NAVSPASUR element sets would change this but they
    // should be restricted to OLEs (as SGP4 is to TLEs)
  if (orbitParams.getCoordinateType() == CoordType::keplerian  &&
      orbitParams.getReferenceFrameType() == FrameType::itrf) {
    throw std::invalid_argument("Orbital elements not compatible with ITRF");
  }
  std::array<double, 6> xeci_array = orbitParams.getInitialState();
  Eigen::Matrix<double, 6, 1> xeciVec;
    // Keplerian to Cartesian - ensured not ITRF above...
  if (orbitParams.getCoordinateType() == CoordType::keplerian) {
    Keplerian kep(xeci_array);
    xeciVec = kep.getCartesian();
  } else {
    for (int ii=0; ii<6; ++ii) {
      xeciVec(ii) = xeci_array[ii];
    }
  }
    // TEME to ITRF (then to GCRF)
  bool teme2itrf {false};
  if (orbitParams.getReferenceFrameType() == FrameType::teme) {
    teme2itrf = true;
    xeciVec = ecfeciSys->teme2ecf(orbitParams.getEpoch(),
                                  xeciVec.block<3, 1>(0, 0),
                                  xeciVec.block<3, 1>(3, 0));
  }
    // ITRF to GCRF - everything is Cartesian by this point unless
    // working with some form of xLE (e.g., a TLE that will be parsed
    // and sent to a specialized propagator.
  if (teme2itrf  ||
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
    } else if (pCfg.getSunGravityModel() == SunGravityModel::eph) {
      std::unique_ptr<Ephemeris> sunEph =
          std::make_unique<Hermite1Eph>("sun",
                                        ceph.at("sun"),
                                        pCfg.getStartTime(),
                                        pCfg.getStopTime(),
                                        ecfeciSys);
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
    } else if (pCfg.getMoonGravityModel() == MoonGravityModel::eph) {
      std::unique_ptr<Ephemeris> moonEph =
          std::make_unique<Hermite1Eph>("moon",
                                        ceph.at("moon"),
                                        pCfg.getStartTime(),
                                        pCfg.getStopTime(),
                                        ecfeciSys);
      std::unique_ptr<ForceModel> moonGrav =
              std::make_unique<ThirdBodyGravity>(phy_const::gm_moon,
                                                 std::move(moonEph));
      deq->addForceModel(std::move(moonGrav));
    }
    if (pCfg.otherGravityModelsEnabled()) {
      for (const auto& planet : ceph) {
        if (planet.first == "moon"  ||  planet.first == "sun") {
          continue;
        }
        double gm_planet {0.0};
        if (planet.first == "mercury") {
          gm_planet = phy_const::gm_mercury;
        } else if (planet.first == "venus") {
          gm_planet = phy_const::gm_venus;
        } else if (planet.first == "mars") {
          gm_planet = phy_const::gm_mars;
        } else if (planet.first == "jupiter") {
          gm_planet = phy_const::gm_jupiter;
        } else if (planet.first == "saturn") {
          gm_planet = phy_const::gm_saturn;
        } else if (planet.first == "uranus") {
          gm_planet = phy_const::gm_uranus;
        } else if (planet.first == "neptune") {
          gm_planet = phy_const::gm_neptune;
        } else if (planet.first == "pluto") {
          gm_planet = phy_const::gm_pluto;
        }
        std::unique_ptr<Ephemeris> planetEph =
            std::make_unique<Hermite1TcEph>(planet.first,
                                            ceph.at(planet.first),
                                            ceph.at("sun"),
                                            pCfg.getStartTime(),
                                            pCfg.getStopTime(),
                                            ecfeciSys);
        std::unique_ptr<ForceModel> planetGrav =
            std::make_unique<ThirdBodyGravity>(gm_planet,
                                               std::move(planetEph));
        deq->addForceModel(std::move(planetGrav));
      }
    }
    if (pCfg.getSrpModel() == SrpModel::spherical) {
      std::unique_ptr<Ephemeris> sunEph = std::make_unique<SunMeeus>(ecfeciSys);
      std::unique_ptr<ForceModel> srp =
          std::make_unique<SrpSpherical>(pCfg.getReflectivity(),
                                         pCfg.getAreaOverMass(),
                                         std::move(sunEph));
      deq->addForceModel(std::move(srp));
    }
      // Integrator
    std::unique_ptr<OdeSolver<JulianDate, double, 6>> sp {nullptr};
    if (pCfg.getPropagator() == Propagator::rk4) {
      sp = std::make_unique<Rk4>(std::move(deq),
                                 pCfg.getStepSize(),
                                 orbitParams.getEpoch(),
                                 xeciVec);
    } else if (pCfg.getPropagator() == Propagator::rk4s) {
      sp = std::make_unique<Rk4s>(std::move(deq),
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
    return std::make_unique<SpEphemeris>(orbitParams.getOrbitName(),
                                         pCfg.getStartTime(),
                                         pCfg.getStopTime(),
                                         ecfeciSys,
                                         std::move(sp));
  } else if (pCfg.getPropagatorType() == PropagatorType::fandg) {
    return std::make_unique<FandG>(orbitParams.getOrbitName(),
                                   orbitParams.getEpoch(),
                                   xeciVec,
                                   ecfeciSys);
  } else if (pCfg.getPropagatorType() == PropagatorType::secular_j2) {
    return std::make_unique<SecularJ2>(orbitParams.getOrbitName(),
                                       orbitParams.getEpoch(),
                                       xeciVec,
                                       ecfeciSys);
  } else if (pCfg.getPropagatorType() == PropagatorType::gpx) {
    return std::make_unique<GpX>(orbitParams.getOrbitName(),
                                 orbitParams.getEpoch(),
                                 xeciVec,
                                 ecfeciSys);
  } else if (pCfg.getPropagatorType() == PropagatorType::kepler1) {
    return std::make_unique<Kepler>(orbitParams.getOrbitName(),
                                    orbitParams.getEpoch(),
                                    xeciVec,
                                    ecfeciSys);
  } else if (pCfg.getPropagatorType() == PropagatorType::kepler1mod) {
    return std::make_unique<KeplerProp>(orbitParams.getOrbitName(),
                                        orbitParams.getEpoch(),
                                        xeciVec,
                                        ecfeciSys);
  } else if (pCfg.getPropagatorType() == PropagatorType::vinti6) {
    return std::make_unique<Vinti>(orbitParams.getOrbitName(),
                                   orbitParams.getEpoch(),
                                   xeciVec,
                                   ecfeciSys);
  } else if (pCfg.getPropagatorType() == PropagatorType::vinti_j2) {
    return std::make_unique<Vinti>(orbitParams.getOrbitName(),
                                   orbitParams.getEpoch(),
                                   xeciVec,
                                   ecfeciSys,
                                   VintiPertModel::J2_ONLY);
  } else if (pCfg.getPropagatorType() == PropagatorType::vinti6mod) {
    return std::make_unique<VintiProp>(orbitParams.getOrbitName(),
                                       orbitParams.getEpoch(),
                                       xeciVec,
                                       ecfeciSys);
  } else if (pCfg.getPropagatorType() == PropagatorType::sgp4) {
    return std::make_unique<Sgp4>(orbitParams.getOrbitName(),
                                  orbitParams.getTle(),
                                  ecfeciSys);
#ifdef GENPL
  } else if (pCfg.getPropagatorType() == PropagatorType::sec_j2) {
    return std::make_unique<SecJ2>(orbitParams.getOrbitName(),
                                   orbitParams.getEpoch(),
                                   xeciVec,
                                   ecfeciSys);
  } else if (pCfg.getPropagatorType() == PropagatorType::osc_j2) {
    return std::make_unique<OscJ2>(orbitParams.getOrbitName(),
                                   orbitParams.getEpoch(),
                                   xeciVec,
                                   ecfeciSys);
#endif
  } else {
    throw std::invalid_argument("Invalid Propagator Type");
  }
}


// Relative orbit definition
std::unique_ptr<Ephemeris>
build_orbit(const RelOrbitDef& relOrbit,
            const OrbitDef& refOrbit,
            const eom::Ephemeris& refEph,
            const std::shared_ptr<const EcfEciSys>& ecfeciSys,
            const std::unordered_map<std::string,
                                     std::vector<eom::state_vector_rec>>& ceph)
{
    // Ephemeris files rejected during input file parsing any before
    // real processing (eomx_parse_input_file) because an OrbitDef is
    // needed to define the propagation method.  Could possibly allow
    // non-propagation based offsets if the need arises.
    //
    // TLE based relative orbits are rejected here for now,
    // unfortunately after potentially time consuming processing has
    // occured.  One option is to attempt to reject earlier - or just
    // allow an offset TLE to be created.
  if (refOrbit.getCoordinateType() == CoordType::tle) {
    throw std::invalid_argument("Relative Orbit not compatible with TLE");
  }
    // Only a single relative orbit definition in RelCoordType exists so
    // no decisions to make.
    // Algorithm for computing orbital elements given radial,
    // transverse, cross-track, and transverse offset distances derived
    // based on "Spacecraft Relative Orbit Geometry Description Through
    // Orbit Element Differences" by Hanspeter Schaub.  This method of
    // defining a bounding box and offset automatically guarantees the
    // energy matching constraint.
  //if (relOrbit.getRelCoordType == rtct)
  Keplerian refOe(refEph.getStateVector(refEph.getEpoch(),
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
                    refEph.getEpoch(),
                    xarr,
                    eom::CoordType::cartesian, eom::FrameType::gcrf);
  return build_orbit(newOrbit, ecfeciSys, ceph);
}


}

