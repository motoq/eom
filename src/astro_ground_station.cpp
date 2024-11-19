/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_ground_station.h>

#include <memory>
#include <utility>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_ground_point.h>

namespace eom {

GroundStation::GroundStation(const JulianDate& jd0,
                             const Eigen::Matrix<double, 3, 1>& pos0f,
                             const Eigen::Matrix<double, 3, 1>& vel0f,
                             std::shared_ptr<const EcfEciSys> ecfeciSys,
                             const std::string& name) :
                             m_jd0 {jd0}, m_pos0f {pos0f}, m_vel0f {vel0f},
                             m_ecfeci {std::move(ecfeciSys)}, m_name {name}
{
  m_jdStart = m_ecfeci->getBeginTime();
  m_jdStop = m_ecfeci->getEndTime();
}


Eigen::Matrix<double, 6, 1>
GroundStation::getStateVector(const JulianDate& jd, EphemFrame frame) const 
{
  Eigen::Matrix<double, 6, 1> xecf;
  xecf.block<3,1>(0,0) = this->getPosition(jd, EphemFrame::ecf);
  xecf.block<3,1>(3,0) = m_vel0f;

  if (frame == EphemFrame::eci) {
    return m_ecfeci->ecf2eci(jd, xecf.block<3,1>(0,0), xecf.block<3,1>(3,0));
  }
  return xecf;
}


Eigen::Matrix<double, 3, 1> GroundStation::getPosition(const JulianDate& jd,
                                                       EphemFrame frame) const
{
  if (frame == EphemFrame::eci  &&
      (jd < m_jdStart  ||  m_jdStop < jd)) {
    throw std::out_of_range("GroundStation::getPosition() - bad time: " +
                            jd.to_string());
  }

  auto dt = phy_const::tu_per_day*(jd - m_jd0);
  Eigen::Matrix<double, 3, 1> posf = m_pos0f + dt*m_vel0f;

  if (frame == EphemFrame::eci) {
    return m_ecfeci->ecf2eci(jd, posf);
  }
  return posf;
}


GroundPoint GroundStation::getGroundPoint(const JulianDate& jd) const
{
  return GroundPoint(this->getPosition(jd, EphemFrame::ecf), m_name);
}


}
