/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_sgp4.h>

#include <array>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <Eigen/Dense>

#include <SGP4.h>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <astro_tle.h>

namespace eom {

Sgp4::Sgp4(const std::string& orbit_name,
           const Tle& tle,
           std::shared_ptr<const EcfEciSys> ecfeciSys,
           Sgp4EarthModel earth_model,
           Sgp4Mode mode)
{
  m_name = orbit_name;
    // Parse TLE
  m_tle = tle;
  m_jd0 = m_tle.getEpoch();
  m_ecfeci = std::move(ecfeciSys);

    // Model/propagator options
  gravconsttype gctype {gravconsttype::wgs72};
  if (earth_model == Sgp4EarthModel::wgs84) {
    gctype = gravconsttype::wgs84;
  }
  char opsmode {'a'};
  if (mode == Sgp4Mode::improved) {
    opsmode = 'i';
  }

    // radians/minute per rev/day = 2pi/1440
  double rpm_per_rpd {utl_const::tpi*utl_const::day_per_min};
    // Initialize SATREC for orbit propagation
  SGP4Funcs::sgp4init(gctype,
                      opsmode,
                      m_tle.getSatN().data(),
                      m_tle.getTleSgpEpoch(),
                      m_tle.getBstar(),
                      rpm_per_rpd*utl_const::day_per_min*
                          m_tle.getMeanMotionRate(),
                      rpm_per_rpd*utl_const::day_per_min*utl_const::day_per_min*
                          m_tle.getMeanMotionSecondRate(),
                      m_tle.getEccentricity(),
                      utl_const::rad_per_deg*m_tle.getArgumentOfPerigee(),
                      utl_const::rad_per_deg*m_tle.getInclination(),
                      utl_const::rad_per_deg*m_tle.getMeanAnomaly(),
                      rpm_per_rpd*m_tle.getMeanMotion(),
                      utl_const::rad_per_deg*m_tle.getRaan(),
                      m_satrec);

    // Check error message with zero forward propagation
  double delta_t {0.0};
  std::array<double, 3> pos;
  std::array<double, 3> vel;
  if (!SGP4Funcs::sgp4(m_satrec, delta_t, pos.data(), vel.data())) {
    throw std::invalid_argument("Sgp4::Sgp4(): Invalid TLE: " +
                                m_tle.getLineOne() + " " +
                                m_tle.getLineTwo());
  }
}


std::string Sgp4::getName() const
{
  return m_name;
}


Eigen::Matrix<double, 6, 1> Sgp4::getStateVector(const JulianDate& jd,
                                                  EphemFrame frame) const
{
  double delta_t {cal_const::min_per_day*(jd - m_jd0)};
  std::array<double, 3> pos;
  std::array<double, 3> vel;
  elsetrec satrec = m_satrec;
  SGP4Funcs::sgp4(satrec, delta_t, pos.data(), vel.data());

  Eigen::Matrix<double, 6, 1> xteme;
  xteme(0) = phy_const::du_per_km*pos[0];
  xteme(1) = phy_const::du_per_km*pos[1];
  xteme(2) = phy_const::du_per_km*pos[2];
  xteme(3) = phy_const::du_per_km*vel[0]*phy_const::sec_per_tu;
  xteme(4) = phy_const::du_per_km*vel[1]*phy_const::sec_per_tu;
  xteme(5) = phy_const::du_per_km*vel[2]*phy_const::sec_per_tu;
  Eigen::Matrix<double, 6, 1> xecf = m_ecfeci->teme2ecf(jd,
                                                        xteme.block<3,1>(0,0),
                                                        xteme.block<3,1>(3,0));

  if (frame == EphemFrame::eci) {
    return m_ecfeci->ecf2eci(jd, xecf.block<3,1>(0,0), xecf.block<3,1>(3,0));
  }
  return xecf;
}


Eigen::Matrix<double, 3, 1> Sgp4::getPosition(const JulianDate& jd,
                                               EphemFrame frame) const
{
  double delta_t {cal_const::min_per_day*(jd - m_jd0)};
  std::array<double, 3> pos;
  std::array<double, 3> vel;
  elsetrec satrec = m_satrec;
  SGP4Funcs::sgp4(satrec, delta_t, pos.data(), vel.data());

  Eigen::Matrix<double, 3, 1> xteme;
  xteme(0) = phy_const::du_per_km*pos[0];
  xteme(1) = phy_const::du_per_km*pos[1];
  xteme(2) = phy_const::du_per_km*pos[2];

  Eigen::Matrix<double, 3, 1> xecf = m_ecfeci->teme2ecf(jd, xteme);

  if (frame == EphemFrame::eci) {
    return m_ecfeci->ecf2eci(jd, xecf);
  }
  return xecf;
}


}
