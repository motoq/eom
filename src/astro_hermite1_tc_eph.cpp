/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_hermite1_tc_eph.h>

#include <string>
#include <memory>
#include <vector>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_hermite1_eph.h>

namespace eom {

Hermite1TcEph::Hermite1TcEph(const std::string& name,
                             const std::vector<state_vector_rec>& target_recs,
                             const std::vector<state_vector_rec>& center_recs,
                             const JulianDate& jdStart,
                             const JulianDate& jdStop,
                             std::shared_ptr<const EcfEciSys> ecfeciSys)
{
    // Start with generating ephemerides
  m_targetEph = std::make_unique<Hermite1Eph>("target",
                                              target_recs,
                                              jdStart,
                                              jdStop,
                                              ecfeciSys);
  m_centerEph = std::make_unique<Hermite1Eph>("center",
                                              center_recs,
                                              jdStart,
                                              jdStop,
                                              ecfeciSys);
    // Finish by setting class variables
  m_name = name;
  m_jdStart = m_targetEph->getBeginTime();
  m_jdStop = m_targetEph->getEndTime();
  m_jdEpoch = m_targetEph->getEpoch();
  m_ecfeciSys = ecfeciSys;
}


Eigen::Matrix<double, 6, 1>
Hermite1TcEph::getStateVector(const JulianDate& jd, EphemFrame frame) const
{
  Eigen::Matrix<double, 6, 1> pv;
  if (frame == EphemFrame::eci) {
    pv = m_centerEph->getStateVector(jd, EphemFrame::eci) +
         m_targetEph->getStateVector(jd, EphemFrame::eci);
  } else {
    pv = m_centerEph->getStateVector(jd, EphemFrame::ecf) +
         m_targetEph->getStateVector(jd, EphemFrame::ecf);
  }

  return pv;
}


Eigen::Matrix<double, 3, 1> Hermite1TcEph::getPosition(const JulianDate& jd,
                                                       EphemFrame frame) const
{
  Eigen::Matrix<double, 3, 1> pos;
  if (frame == EphemFrame::eci) {
    pos = m_centerEph->getPosition(jd, EphemFrame::eci) +
          m_targetEph->getPosition(jd, EphemFrame::eci);
  } else {
    pos = m_centerEph->getPosition(jd, EphemFrame::ecf) +
          m_targetEph->getPosition(jd, EphemFrame::ecf);
  }

  return pos;
}


}
