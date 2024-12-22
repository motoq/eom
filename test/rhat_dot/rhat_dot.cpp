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

int main()
{

  using namespace eom;
  using namespace utl_units;

    // Initial state
  GregDate gd {2021, 11, 12};
  JulianDate jd1 {gd, 17, 0, 0.0};
  std::array<double, 6> oe = {4.1632, 0.741,
                              63.4_deg, 345.0_deg, 270.0_deg, 0.0_deg};

    // Create ECFECI service
  Duration dur {1.0, 1.0_day};
  auto jd2 = jd1 + dur;
  Duration fi_dt {0.0, 1.0};
  auto f2i = std::make_shared<EcfEciSys>(jd1, jd2, fi_dt, nullptr);

  PropagatorConfig pcfg {PropagatorType::sp};
  pcfg.setStartStopTime(jd1, jd2);
  pcfg.setGravityModel(GravityModel::jn);
  pcfg.setDegreeOrder(2, 0);
  pcfg.setSunGravityModel(SunGravityModel::meeus);
  pcfg.setMoonGravityModel(MoonGravityModel::meeus);
  pcfg.setPropagator(Propagator::adams4);
  pcfg.setStepSize({1.0, 1.0_min});

  OrbitDef odef {"heo_sat", pcfg, jd1, oe, 
                 CoordType::keplerian, FrameType::gcrf};
  std::unordered_map<std::string,
                     std::vector<state_vector_rec>> ceph;

  std::unique_ptr<Ephemeris> eph = build_orbit(odef, f2i, ceph);
  GravityJn j2Grav {2};

  Duration dt {5.0, 1.0_min};
    // Make sure ephemeris covers time span that will be needed
    // to generate derivatives
  auto jd = jd1 + dt.getDays();
  auto jdStop = jd2 + -1.0*dt.getDays();
  Duration dx {1.0, 0.1_sec};
  auto dx_days = dx.getDays();
  auto dx_tu = dx.getTu();
  double max_rhat_dot_err {0.0};
  double max_rhat_ddot_err {0.0};
  double max_rhat_dot_func_err {0.0};
  while (jd < jdStop) {
    Eigen::Matrix<double, 6, 1> pvf = eph->getStateVector(jd, EphemFrame::ecf);
    Eigen::Matrix<double, 3, 1> r_s_o_f = pvf.block<3,1>(0,0);
    Eigen::Matrix<double, 3, 1> v_s_f_f = pvf.block<3,1>(3,0);
    Eigen::Matrix<double, 3, 1> a_s_i_f = j2Grav.getAcceleration(r_s_o_f);
    Eigen::Matrix<double, 3, 1> a_s_f_f = f2i->gravity2ecf(jd, r_s_o_f,
                                                               v_s_f_f,
                                                               a_s_i_f);
    UnitVector<double, 3> uv {r_s_o_f, v_s_f_f, a_s_f_f};
    Eigen::Matrix<double, 3, 1> rhat_dot = uv.getNormalizedDot();
    Eigen::Matrix<double, 3, 1> rhat_ddot = uv.getNormalizedDDot();

    Eigen::Matrix<double, 3, 1> rhat = r_s_o_f.normalized();
    Eigen::Matrix<double, 3, 1> rhatf =
        eph->getPosition(jd + +1.0*dx_days, EphemFrame::ecf);
    rhatf.normalize();
    Eigen::Matrix<double, 3, 1> rhatb =
        eph->getPosition(jd + -1.0*dx_days, EphemFrame::ecf);
    rhatb.normalize();

    Eigen::Matrix<double, 3, 1> rhat_dot_num =
        derivative::first(dx_tu, rhatb, rhatf);
    Eigen::Matrix<double, 3, 1> rhat_ddot_num =
        derivative::second(dx_tu, rhatb, rhat, rhatf);

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

    jd += dt;
  }
  std::cout << "\nMax rhat_dot error:   " << std::abs(max_rhat_dot_err);
  std::cout << "\nMax rhat_ddot error:  " << std::abs(max_rhat_ddot_err);
  std::cout << "\nMax function diff:    " << std::abs(max_rhat_dot_func_err);

  std::cout << '\n';
  
}
