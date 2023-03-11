/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_sp3_orbit.h>

#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <mth_hermite2.h>
#include <astro_ephemeris.h>
#include <astro_gravity_jn.h>
#include <mth_index_mapper.h>

namespace eom {

Sp3Orbit::Sp3Orbit(const std::string& name,
                   const std::vector<state_vector_rec>& sp3_records,
                   const JulianDate& jdStart,
                   const JulianDate& jdStop,
                   const std::shared_ptr<const EcfEciSys>& ecfeciSys)
{
  m_name = name;
  m_ecfeciSys = ecfeciSys;

  if (sp3_records.size() > 1) {
    m_jdEpoch = sp3_records.front().t;
    m_jdStart = m_jdEpoch;
    m_jdStop = sp3_records.back().t;
  } else {
    throw std::runtime_error(
        "Sp3Orbit::Sp3Orbit() Less than two eph records: " + m_name);
  }
  if (jdStart < m_jdStart) {
    throw std::runtime_error(
        "Sp3Orbit::Sp3Orbit() Ephemeris begins too late: " + m_name);
  }
  if (m_jdStop < jdStop) {
    throw std::runtime_error(
        "Sp3Orbit::Sp3Orbit() Ephemeris ends too early: " + m_name);
  }

  std::vector<std::pair<JulianDate, JulianDate>> times;
  GravityJn grv(4);
    // Generate and store Hermite interpolation objects
  for (unsigned long ii=1UL; ii<sp3_records.size(); ++ii) {
    const state_vector_rec& r1 = sp3_records[ii-1UL];
    const state_vector_rec& r2 = sp3_records[ii];
    Eigen::Matrix<double, 3, 1> a1 = grv.getAcceleration(r1.p);
    a1 = m_ecfeciSys->gravity2ecf(r1.t, r1.p, r1.v, a1);
    Eigen::Matrix<double, 3, 1> a2 = grv.getAcceleration(r2.p);
    a2 = m_ecfeciSys->gravity2ecf(r2.t, r2.p, r2.v, a2);
    double dt_tu {phy_const::tu_per_day*(r2.t - r1.t)};
    Hermite2<double, 3> hItp(dt_tu,
                             r1.p, r1.v, a1,
                             r2.p, r2.v, a2,
                             phy_const::epsdt);
    m_eph_interpolators.emplace_back(r1.t, r2.t, hItp);
    times.emplace_back(r1.t, r2.t);
  }
  m_ndxr = std::make_unique<IndexMapper<JulianDate>>(times);
}


Eigen::Matrix<double, 6, 1> Sp3Orbit::getStateVector(const JulianDate& jd,
                                                        EphemFrame frame) const
{
  unsigned long ndx {};
  try {
    ndx = m_ndxr->getIndex(jd);
  } catch (const std::out_of_range& ia) {
    throw std::out_of_range("Sp3Orbit::getStateVector() - bad time");
  }
  const auto& irec = m_eph_interpolators[ndx];
  double dt_tu {phy_const::tu_per_day*(jd - irec.jd1)};
  Eigen::Matrix<double, 6, 1> xecf;
  xecf.block<3,1>(0,0) = irec.hItp.getPosition(dt_tu);
  xecf.block<3,1>(3,0) = irec.hItp.getVelocity(dt_tu);

  if (frame == EphemFrame::eci) {
    return m_ecfeciSys->ecf2eci(jd, xecf.block<3,1>(0,0), xecf.block<3,1>(3,0));
  }

  return xecf;
}


Eigen::Matrix<double, 3, 1> Sp3Orbit::getPosition(const JulianDate& jd,
                                                  EphemFrame frame) const
{
  unsigned long ndx {};
  try {
    ndx = m_ndxr->getIndex(jd);
  } catch (const std::out_of_range& ia) {
    throw std::out_of_range("Sp3Orbit::getPosition - bad time");
  }
  const auto& irec = m_eph_interpolators[ndx];
  double dt_tu {phy_const::tu_per_day*(jd - irec.jd1)};
  Eigen::Matrix<double, 3, 1> xecf = irec.hItp.getPosition(dt_tu);

  if (frame == EphemFrame::eci) {
    return m_ecfeciSys->ecf2eci(jd, xecf);
  }

  return xecf;
}


}
