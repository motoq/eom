/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_sp3_ephem.h>

#include <string>
#include <array>
#include <vector>
#include <memory>
#include <utility>
#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_granule.h>
#include <astro_ephemeris.h>
#include <mth_index_mapper.h>

namespace eom {

Sp3Ephem::Sp3Ephem(const std::string& name,
                   const std::vector<state_vector_rec>& sp3_records,
                   const JulianDate& jdStart,
                   const JulianDate& jdStop,
                   const std::shared_ptr<const EcfEciSys>& ecfeciSys)
{
  m_name = name;
  m_ecfeciSys = ecfeciSys;

  unsigned long npts {static_cast<unsigned long>(sp3::np)};
  if (sp3_records.size() >= npts) {
    m_jdEpoch = sp3_records.front().t;
    m_jdStart = m_jdEpoch;
    m_jdStop = sp3_records.back().t;
  } else {
    throw std::runtime_error(
        "Sp3Ephem::Sp3Ephem() insufficient eph records: " + m_name);
  }
  if (jdStart < m_jdStart) {
    throw std::runtime_error(
        "Sp3Ephem::Sp3Ephem() Ephemeris begins too late: " + m_name);
  }
  if (m_jdStop < jdStop) {
    throw std::runtime_error(
        "Sp3Ephem::Sp3Ephem() Ephemeris ends too early: " + m_name);
  }

  unsigned long nrec = sp3_records.size()/static_cast<unsigned long>(npts-1UL);
  nrec--;
  std::vector<std::pair<JulianDate, JulianDate>> times;
    // Generate and store Chebyshev granules - separate position and
    // velocity coefficients
  std::array<JulianDate, sp3::np> jds;
  Eigen::Matrix<double, 3, sp3::np> pvecs;
  Eigen::Matrix<double, 3, sp3::np> vvecs;
  for (unsigned long ii=0; ii<nrec; ++ii) {
    unsigned long ndx {ii*(npts-1UL)};
    for (int jj=0; jj<sp3::np; ++jj) {
      const state_vector_rec& erec = sp3_records[ndx];
      jds[jj] = erec.t;
      pvecs.block(0,jj,3,1) = erec.p;
      vvecs.block(0,jj,3,1) = erec.v;
      ndx++;
    }
    Granule<sp3::order, sp3::np> tItp(jds, pvecs, vvecs);
    m_eph_interpolators.emplace_back(jds[0], jds[sp3::np-1], tItp);
    times.emplace_back(jds[0], jds[sp3::np-1]);
  }
  m_ndxr = std::make_unique<IndexMapper<JulianDate>>(times);
}


Eigen::Matrix<double, 6, 1> Sp3Ephem::getStateVector(const JulianDate& jd,
                                                     EphemFrame frame) const
{
  unsigned long ndx {};
  try {
    ndx = m_ndxr->getIndex(jd);
  } catch (const std::out_of_range& ia) {
    throw std::out_of_range("Sp3Ephem::getStateVector() - bad time");
  }
  const auto& irec = m_eph_interpolators[ndx];
  Eigen::Matrix<double, 6, 1> xecf;
  xecf.block<3,1>(0,0) = irec.tItp.getPosition(jd);
  xecf.block<3,1>(3,0) = irec.tItp.getVelocity(jd);

  if (frame == EphemFrame::eci) {
    return m_ecfeciSys->ecf2eci(jd, xecf.block<3,1>(0,0), xecf.block<3,1>(3,0));
  }

  return xecf;
}


Eigen::Matrix<double, 3, 1> Sp3Ephem::getPosition(const JulianDate& jd,
                                                  EphemFrame frame) const
{
  unsigned long ndx {};
  try {
    ndx = m_ndxr->getIndex(jd);
  } catch (const std::out_of_range& ia) {
    throw std::out_of_range("Sp3Ephem::getPosition - bad time");
  }
  const auto& irec = m_eph_interpolators[ndx];
  Eigen::Matrix<double, 3, 1> xecf = irec.tItp.getPosition(jd);

  if (frame == EphemFrame::eci) {
    return m_ecfeciSys->ecf2eci(jd, xecf);
  }

  return xecf;
}


}
