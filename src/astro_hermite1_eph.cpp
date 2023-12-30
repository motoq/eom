/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_hermite1_eph.h>

#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <utility>
#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <mth_hermite1.h>
#include <astro_ephemeris.h>
#include <mth_index_mapper.h>

namespace eom {

Hermite1Eph::Hermite1Eph(const std::string& name,
                         const std::vector<state_vector_rec>& pv_records,
                         const JulianDate& jdStart,
                         const JulianDate& jdStop,
                         std::shared_ptr<const EcfEciSys> ecfeciSys)
{
  m_name = name;
  m_ecfeciSys = std::move(ecfeciSys);

  if (pv_records.size() > 1) {
    m_jdEpoch = pv_records.front().t;
    m_jdStart = m_jdEpoch;
    m_jdStop = pv_records.back().t;
  } else {
    throw std::runtime_error(
        "Hermite1Eph::Hermite1Eph() Less than two eph records: " + m_name);
  }
  if (jdStart < m_jdStart) {
    throw std::runtime_error(
        "Hermite1Eph::Hermite1Eph() Ephemeris begins too late: " + m_name);
  }
  if (m_jdStop < jdStop) {
    throw std::runtime_error(
        "Hermite1Eph::Hermite1Eph() Ephemeris ends too early: " + m_name);
  }

  std::vector<std::pair<JulianDate, JulianDate>> times;
    // Generate and store Hermite interpolation objects
  for (unsigned long ii=1UL; ii<pv_records.size(); ++ii) {
    const state_vector_rec& r1 = pv_records[ii-1UL];
    const state_vector_rec& r2 = pv_records[ii];
    double dt_tu {phy_const::tu_per_day*(r2.t - r1.t)};
    Hermite1<double, 3> hItp(dt_tu,
                             r1.p, r1.v,
                             r2.p, r2.v,
                             phy_const::epsdt);
    m_eph_interpolators.emplace_back(r1.t, r2.t, hItp);
    times.emplace_back(r1.t, r2.t);
  }
  m_ndxr = std::make_unique<IndexMapper<JulianDate>>(times);
}


Eigen::Matrix<double, 6, 1> Hermite1Eph::getStateVector(const JulianDate& jd,
                                                        EphemFrame frame) const
{
  unsigned long ndx {};
  try {
    ndx = m_ndxr->getIndex(jd);
  } catch (const std::out_of_range& ia) {
    throw std::out_of_range("Hermite1Eph::getStateVector() - bad time");
  }
  const auto& irec = m_eph_interpolators[ndx];
  double dt_tu {phy_const::tu_per_day*(jd - irec.jd1)};
  Eigen::Matrix<double, 6, 1> xeci;
  xeci.block<3,1>(0,0) = irec.hItp.getPosition(dt_tu);
  xeci.block<3,1>(3,0) = irec.hItp.getVelocity(dt_tu);

  if (frame == EphemFrame::ecf) {
    return m_ecfeciSys->eci2ecf(jd, xeci.block<3,1>(0,0), xeci.block<3,1>(3,0));
  }

  return xeci;
}


Eigen::Matrix<double, 3, 1> Hermite1Eph::getPosition(const JulianDate& jd,
                                                     EphemFrame frame) const
{
  unsigned long ndx {};
  try {
    ndx = m_ndxr->getIndex(jd);
  } catch (const std::out_of_range& ia) {
    throw std::out_of_range("Hermite1Eph::getPosition - bad time");
  }
  const auto& irec = m_eph_interpolators[ndx];
  double dt_tu {phy_const::tu_per_day*(jd - irec.jd1)};
  Eigen::Matrix<double, 3, 1> xeci = irec.hItp.getPosition(dt_tu);

  if (frame == EphemFrame::ecf) {
    return m_ecfeciSys->eci2ecf(jd, xeci);
  }

  return xeci;
}


}
