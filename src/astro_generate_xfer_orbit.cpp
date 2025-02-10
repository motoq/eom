/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>

#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Eigen/Dense>

#include <phy_const.h>
#include <astro_build.h>
#include <cal_julian_date.h>
#include <cal_duration.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_fandg.h>
#include <astro_keplerian.h>
#include <astro_orbit_def.h>
#include <astro_propagator_config.h>

#include <astro_generate.h>

namespace {
  constexpr int maxitr {100};
}

namespace eom {

/*
 * BMW Gauss problem via universal variables
 */
Eigen::Matrix<double, 6, 1>
generate_gauss_fg_xfer(const Eigen::Matrix<double, 2, 1>& r1,
                       const Eigen::Matrix<double, 2, 1>& r2,
                       const Duration& dur)
{
  // Units of DU and TU
  //   GM = 1.0

    // Short way around placeholder
  constexpr double dm {1.0};

  Eigen::Matrix<double, 6, 1> pv {Eigen::Matrix<double, 6, 1>::Zero()};

  auto r1mag = r1.norm();
  auto r2mag = r2.norm();
  auto r1dotr2 = r1.dot(r2);
  auto avar = dm*std::sqrt(r1mag*r2mag + r1dotr2);

  double zvar {0.0};
  for (int ii=0; ii<10; ++ii) {
    auto [cz, sz] = astro_fg_cands(zvar);
    auto yvar = r1mag + r2mag - avar*(1.0 - zvar*sz)/std::sqrt(cz);
    if (yvar < 0.0) {
      std::cerr << "\n\nNegative yvar";
      return pv;
    }
    auto sqrty = std::sqrt(yvar);
    auto xvar = sqrty/std::sqrt(cz);
    auto xxx = xvar*xvar*xvar;
  
      // Transfer time vs. computed time
    auto dt = dur.getTu();
    auto dtc = xxx*sz + avar*sqrty;
    auto [dcdz, dsdz] = astro_fg_dcands_dz(zvar);
    auto cinv = 1.0/cz;
    auto dtdz = xxx*(dsdz - 1.5*sz*dcdz*cinv) +
                0.125*avar*(3.0*sz*sqrty*cinv + avar/xvar);
    auto dz = (dt - dtc)/dtdz;
    zvar += dz;
    std::cerr << "\ndz  " << dz;

    if (std::abs(dz) < 1.0e-6) {
      std::cerr << "\n\nFound zvar,  dt:  " << dt << "  dtc:  " << dtc;
      break;
    }
  }
  std::cerr << "\nzvar  " << zvar;

  auto [cz, sz] = astro_fg_cands(zvar);
  auto yvar = r1mag + r2mag - avar*(1.0 - zvar*sz)/std::sqrt(cz);
  double f = {1.0 - yvar/r1mag};
  double g = {avar*std::sqrt(yvar)};

  pv.block<2,1>(0,0) = r1;
  pv.block<2,1>(3,0) = (r2 - f*r1)/g;

  return pv;
}

/*
 * Based on William E. Wiesel's "Modern Astrodynamics", 2nd Ed.,
 * Section 1.5 Boundary Value Problems, this function illustrates
 * the "shooting method".
 *
 * Current implementation uses the start orbit's state vector at the
 * beginning of the time at which to start the transfer as the initial
 * guess.  Instantaneous delta-v is assumed vs. any sort of finite burn
 * model.  Proper convergence at this time should use 1/2 the smaller
 * period of the initial or destination orbit as the maximum transfer
 * time (this is implemented as a direct transfer, not a more efficient
 * "long way around" multi-rev transfer).  Little things to add include
 * an option to prime the pump via 2-body Lambert's problem, efficient
 * multi-rev transfer, etc.
 *
 * The primary purpose is to illustrate this rather cool "shooting
 * method" described by Wiesel illustrating differential correction
 * techniques.  Even with a sloppy initial guess, proper bounding
 * keeps convergence in check and the method doesn't really care
 * about what orbit propagator is used.  The technique is incredibly
 * simple.
 */
std::pair<std::unique_ptr<Ephemeris>, int>
generate_xfer_orbit(const std::string& orbit_name,
                    const Ephemeris& startOrbit,
                    const Ephemeris& endOrbit,
                    const JulianDate& xferStartTime,
                    const Duration& xferDur,
                    const PropagatorConfig& propCfg,
                    std::shared_ptr<const EcfEciSys> ecfeciSys)
{

    // Copy propagator configuration and switch to analytic methods
    // for sun/moon vs. reloading file based ephemerides.  May add
    // in later but of little benefit at this point given instant
    // delta-v model and half rev constraint on transfer time.
  PropagatorConfig xferPropCfg = propCfg;
  xferPropCfg.disableOtherGravityModels();
  if (xferPropCfg.getSunGravityModel() == SunGravityModel::eph) {
    xferPropCfg.setSunGravityModel(SunGravityModel::meeus);
  }
  if (xferPropCfg.getMoonGravityModel() == MoonGravityModel::eph) {
    xferPropCfg.setMoonGravityModel(MoonGravityModel::meeus);
  }

    // Inverse of partials of end position vector w.r.t. start velocity
    // vector.  See Vallado's Fundamentals of Astrodynamics and Applications,
    // Orbit Determination and Estimation chapter to see the partials
    // spelled out.
  Eigen::Matrix<double, 3, 3> dv1_dr2 = (1.0/xferDur.getTu())*
                                        Eigen::Matrix<double, 3, 3>::Identity();
    // Desired end state - find velocity needed at
    // start time to reach this location at end time.
  auto xferEndTime = xferStartTime + xferDur;
  Eigen::Matrix<double, 3, 1> r2 = endOrbit.getPosition(xferEndTime,
                                                        EphemFrame::eci);

    // Initial guess is state vector at start time
  Eigen::Matrix<double, 6, 1> rv = startOrbit.getStateVector(xferStartTime,
                                                             EphemFrame::eci);


    // r1 to r2
  Eigen::Matrix<double, 3, 1> ihat {rv.block<3,1>(0,0)};
  Eigen::Matrix<double, 3, 1> khat {ihat.cross(r2)};
  Eigen::Matrix<double, 3, 1> jhat {khat.cross(ihat)};
  ihat.normalize();
  jhat.normalize();
  khat.normalize();
  Eigen::Matrix<double, 3, 3> cp;
  cp.row(0) = ihat;
  cp.row(1) = jhat;
  cp.row(2) = khat;
  Eigen::Matrix<double, 3, 1> r1p = cp*rv.block<3,1>(0,0);
  Eigen::Matrix<double, 3, 1> r2p = cp*r2;
  Eigen::Matrix<double, 6, 1> rv2d = generate_gauss_fg_xfer(r1p.block<2,1>(0,0),
                                                            r2p.block<2,1>(0,0),
                                                            xferDur);
  cp.transposeInPlace();
  rv.block<3,1>(0,0) = cp*rv2d.block<3,1>(0,0);
  rv.block<3,1>(3,0) = cp*rv2d.block<3,1>(3,0);



    // Dummy parameter
  std::unordered_map<std::string, std::vector<eom::state_vector_rec>> ceph;
    // Transfer ephemeris to determine
  std::unique_ptr<Ephemeris> xeph {nullptr};
    // Use aggressive scaling of initial correction on first go
  double bnds {.25};
  double old_miss {1.0};
  int nitr {-1};
  for (int ii=0; ii<maxitr; ++ii) {
    std::array<double, 6> x1 = {rv(0), rv(1), rv(2), rv(3), rv(4), rv(5)};
      // Update orbit def to create new xeph based on current guess
    OrbitDef orbit {orbit_name,
                    xferPropCfg,
                    xferStartTime,
                    x1,
                    CoordType::cartesian,
                    FrameType::gcrf};
    xeph = build_orbit(orbit, ecfeciSys, ceph);
    Eigen::Matrix<double, 3, 1> r2x = xeph->getPosition(xferEndTime,
                                                        EphemFrame::eci);
      // Update
    Eigen::Matrix<double, 3, 1> dr2 = r2 - r2x;
    double miss {dr2.norm()};
    if (miss < 1.0*phy_const::du_per_m) {
      nitr = ii + 1;
      break;
    } else if (miss > old_miss) {
      bnds /= 2.0;
    } else {
      old_miss = miss;
      bnds *= 1.5;
      //bnds = bnds > 0.75 ? 0.75 : bnds;
      bnds = bnds > 1.0 ? 1.0 : bnds;
    }
    Eigen::Matrix<double, 3, 1> dv1 = bnds*dv1_dr2*dr2;
    rv.block<3,1>(3,0) += dv1;
  }

  return std::make_pair(std::move(xeph), nitr);
}


}
