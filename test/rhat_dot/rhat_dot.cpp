/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <vector>

#include <Eigen/Dense>

#include <phy_const.h>
#include <utl_units.h>
#include <mth_derivative.h>
#include <mth_unit_vector.h>
#include <cal_duration.h>
#include <astro_ephemeris.h>
#include <cal_greg_date.h>
#include <cal_julian_date.h>
#include <astro_build.h>
#include <astro_ecfeci_sys.h>
#include <astro_gravity_jn.h>
#include <astro_orbit_def.h>
#include <astro_propagator_config.h>

//
// Uncertainty in position vector when testing covariance transformation
// (diagonal elements of the square root of the position covariance)
// Units of ER.
//
namespace {
  constexpr double k952D {2.448};           // 95% 2D Mahalanobis distance
  constexpr double k953D {2.796};           // 95% 3D Mahalanobis distance

  constexpr double sigma_x {0.03};
  constexpr double sigma_y {0.07};
  constexpr double sigma_z {0.05};
}

/**
 * Tests the algorithms for taking the 1st and 2nd derivatives of
 * a unit vector as implemented in mth_unit_vector.h.  A HEO orbit is
 * created and the earth fixed position, velocity, and acceleration
 * vectors are used.  The primary weakness to this test is the majority
 * of the acceleration is centripetal, therefore not rigorously testing
 * one of the components of the 2nd derivative.
 *
 * Another test is done to validate the mapping from the position
 * vector to its unit vector via partial derivatives.  Perturbations are
 * added to the position vector.  Expected containment of this perturbed
 * vector is verified via Mahalanobis distance as a sanity check.  Next,
 * the unit vector is created, the covariance transformed, and the
 * Mahalanobis distance once again used to check containment.
 *
 * @author  Kurt Motekew
 * @date    2025/01/03
 */
int main()
{
  using namespace eom;
  using namespace utl_units;

    // Random number generator
  std::random_device r;
  std::seed_seq seed {r()};
  std::mt19937 sgen(seed);
  std::normal_distribution<> randn(0.0, 1.0);
    // Position covariance - units of DU (~ER)
  Eigen::Matrix<double, 3, 3> Cov = Eigen::Matrix<double, 3, 3>::Zero();
  Cov(0,0) = sigma_x*sigma_x;
  Cov(1,1) = sigma_y*sigma_y;
  Cov(2,2) = sigma_z*sigma_z;
  Eigen::Matrix<double, 3, 3> W = Cov.inverse();

    // Satellite state vector
  GregDate gd {2021, 11, 12};
  JulianDate jd1 {gd, 17, 0, 0.0};
  std::array<double, 6> oe = {4.1632, 0.741,
                              63.4_deg, 345.0_deg, 270.0_deg, 0.0_deg};

    // Create ECFECI transformation service
  Duration dur {1.0, 1.0_day};
  auto jd2 = jd1 + dur;
  Duration fi_dt {0.0, 1.0};
  auto f2i = std::make_shared<EcfEciSys>(jd1, jd2, fi_dt, nullptr);

    // J2 gravity model plus sun/moon
  PropagatorConfig pcfg {PropagatorType::sp};
  pcfg.setStartStopTime(jd1, jd2);
  pcfg.setGravityModel(GravityModel::jn);
  pcfg.setDegreeOrder(2, 0);
  pcfg.setSunGravityModel(SunGravityModel::meeus);
  pcfg.setMoonGravityModel(MoonGravityModel::meeus);
  pcfg.setPropagator(Propagator::adams4);
    // Tight integration step size for perigee since using simple integrator
  pcfg.setStepSize({0.25, 1.0_min});

  OrbitDef odef {"heo_sat", pcfg, jd1, oe, 
                 CoordType::keplerian, FrameType::gcrf};
    // Orbit builder utility requires source of external celestial
    // ephemerides even if not used
  std::unordered_map<std::string,
                     std::vector<state_vector_rec>> ceph;
  std::unique_ptr<Ephemeris> eph = build_orbit(odef, f2i, ceph);
    // Local gravity model since ephemeris services only
    // provide position and velocity
  GravityJn j2Grav {2};

    // Output step size
  Duration dt {1.0, 1.0_min};
    // Make sure ephemeris covers time span that will be needed
    // to generate derivatives
  auto jd = jd1 + dt.getDays();
  auto jdStop = jd2 + -1.0*dt.getDays();
    // Numeric differentiation step size
  Duration dx {1.0, 0.1_sec};
  auto dx_days = dx.getDays();
  auto dx_tu = dx.getTu();
    // Accumulated error in rhat_dot and rhat_ddot
  double max_rhat_dot_err {0.0};
  double max_rhat_ddot_err {0.0};
  double max_rhat_dot_func_err {0.0};
    // Used to track containment from mapping position error
    // to unit vector
  int n_contained {0};
  int nhat_contained {0};
  int npts {0};
  while (jd < jdStop) {
    //
    // First check derivatives
    //
      // ECEF position, earth fixed velocity and acceleration
    Eigen::Matrix<double, 6, 1> pvf = eph->getStateVector(jd, EphemFrame::ecf);
    Eigen::Matrix<double, 3, 1> r_s_o_f = pvf.block<3,1>(0,0);
    Eigen::Matrix<double, 3, 1> v_s_f_f = pvf.block<3,1>(3,0);
    Eigen::Matrix<double, 3, 1> a_s_i_f = j2Grav.getAcceleration(r_s_o_f);
    Eigen::Matrix<double, 3, 1> a_s_f_f = f2i->gravity2ecf(jd, r_s_o_f,
                                                               v_s_f_f,
                                                               a_s_i_f);
      // derivatives via unit vector class
    UnitVector<double, 3> uv {r_s_o_f, v_s_f_f, a_s_f_f};
    Eigen::Matrix<double, 3, 1> rhat_dot = uv.getNormalizedDot();
    Eigen::Matrix<double, 3, 1> rhat_ddot = uv.getNormalizedDDot();

      // Numerical differentiation for comparison/validation
    const Eigen::Matrix<double, 3, 1> rhat = r_s_o_f.normalized();
    Eigen::Matrix<double, 3, 1> rhatf =
        eph->getPosition(jd + +1.0*dx_days, EphemFrame::ecf);
    rhatf.normalize();
    Eigen::Matrix<double, 3, 1> rhatb =
        eph->getPosition(jd + -1.0*dx_days, EphemFrame::ecf);
    rhatb.normalize();
      //
    Eigen::Matrix<double, 3, 1> rhat_dot_num =
        derivative::first(dx_tu, rhatb, rhatf);
    Eigen::Matrix<double, 3, 1> rhat_ddot_num =
        derivative::second(dx_tu, rhatb, rhat, rhatf);

      // Also validate rhat_dot stand-alone function
    auto drdhd = -1.0*(rhat_dot - rhat_dot_num).norm();
    auto drdhdd = -1.0*(rhat_ddot - rhat_ddot_num).norm();
    auto drdhfd = -1.0*(rhat_dot - unit_vector_dot(r_s_o_f, v_s_f_f)).norm();
    if (drdhd < max_rhat_dot_err) {
      max_rhat_dot_err = drdhd;
    }
    if (drdhdd < max_rhat_ddot_err) {
      max_rhat_ddot_err = drdhdd;
    }
    if (drdhfd < max_rhat_dot_func_err) {
      max_rhat_dot_func_err = drdhfd;
    }

    //
    // Test covariance mapping
    //

      // Add noise to position vector - check 95% containment
    Eigen::Matrix<double, 3, 1> pert = {sigma_x*randn(sgen),
                                        sigma_y*randn(sgen),
                                        sigma_z*randn(sgen)};
    Eigen::Matrix<double, 3, 1> r_pert = r_s_o_f + pert;
    Eigen::Matrix<double, 3, 1> dr = r_s_o_f - r_pert;
    if (std::sqrt((dr.transpose()*W*dr)(0,0)) <= k953D) {
      n_contained++;
    }
      // Transform Cov to unit vector covariance
      // Now rank deficient
    Eigen::Matrix<double, 3, 3> A = unit_vector_partials(r_s_o_f);
    Eigen::Matrix<double, 3, 3> Cov_hat_rd = A*Cov*A.transpose();
      // Becomes rank deficient due to projection, so can't invert
      // for Mahalanobis in this reference frame/dimension
      // Transform to plane normal to position vector, and reduce to
      // 2D
    Eigen::Matrix<double, 3, 1> zhat = {0.0, 0.0, 1.0};
    Eigen::Matrix<double, 3, 1> khat {rhat};
    Eigen::Matrix<double, 3, 1> ihat {zhat.cross(khat)};
    Eigen::Matrix<double, 3, 1> jhat {khat.cross(ihat)};
    ihat.normalize();
    jhat.normalize();
    khat.normalize();
    Eigen::Matrix<double, 3, 3> Q;
    Q.row(0) = ihat;
    Q.row(1) = jhat;
    Q.row(2) = khat;
      // Upper full rank 2x2 is in plane orthogonal to rhat
      // after transformation - extract that for W
    Eigen::Matrix<double, 3, 3> Cov_hat = Q*Cov_hat_rd*Q.transpose();
    Eigen::Matrix<double, 2, 2> W_hat = Cov_hat.block<2,2>(0,0).inverse();
    Eigen::Matrix<double, 3, 1> rhat_pert = r_pert.normalized();
    Eigen::Matrix<double, 2, 1> drhat = (Q*(rhat - rhat_pert)).block<2,1>(0,0);
    if (std::sqrt((drhat.transpose()*W_hat*drhat)(0,0)) <= k952D) {
      nhat_contained++;
    }

    npts++;
    jd += dt;
  }

  std::cout << "\n  " << npts  << " test points";
  std::cout << "\n--- Test derivative functions ---";
  std::cout << "\nMax rhat_dot error:   " << std::abs(max_rhat_dot_err);
  std::cout << "\nMax rhat_ddot error:  " << std::abs(max_rhat_ddot_err);
  std::cout << "\nMax function diff:    " << std::abs(max_rhat_dot_func_err);

  std::cout << "\n--- Test transformation function ---";
  std::cout << '\n' << (100.0*n_contained)/npts <<
               "\% containment vs expected 95\% for perturbed vector";
  std::cout << '\n' << (100.0*nhat_contained)/npts <<
               "\% containment vs expected 95\% for perturbed unit vector";
  std::cout << '\n';
}
