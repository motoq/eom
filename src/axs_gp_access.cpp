/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_gp_access.h>

#include <string>
#include <utility>
#include <memory>
#include <cmath>
#include <stdexcept>

#include <Eigen/Dense>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ground_point.h>
#include <astro_ephemeris.h>
#include <astro_keplerian.h>
#include <axs_gp_constraints.h>

namespace {
    // One second increment after end of access window to start
    // searching for the next access window
  constexpr double jd_inc {1.0*utl_const::day_per_sec};
    // Search increment
  constexpr double dt_days {20.0*utl_const::day_per_sec};
    // Factor to cut max search time by
  constexpr double dt_max_sf {0.25};
    // Access convergence
  constexpr double tol_dt_day {0.1*utl_const::day_per_sec};
  constexpr int max_itr {42};
}

namespace eom {

GpAccess::GpAccess(const JulianDate& jdStart,
                   const JulianDate& jdStop,
                   const GroundPoint& gp,
                   const GpConstraints& xcs,
                   std::shared_ptr<const Ephemeris> eph) : m_jdStart {jdStart},
                                                           m_jdStop {jdStop},
                                                           m_gp {gp},
                                                           m_xcs {xcs},
                                                           m_eph {
                                                             std::move(eph)
                                                           }
{
  m_jd = jdStart;

  try {
    Keplerian oe(m_eph->getStateVector(m_jd, EphemFrame::eci));
    double max_vel {phy_const::earth_equatorial_speed() + oe.getPerigeeSpeed()};
    m_theta_dot = max_vel/oe.getPerigeeRadius();

    double alpha {utl_const::pio2 + m_xcs.getMinEl()};
    double phi {std::sin(alpha)*
                std::asin(phy_const::earth_smaj/oe.getPerigeeRadius())};
    double theta {utl_const::pi - (alpha + phi)};
    m_max_dt_days = dt_max_sf*phy_const::day_per_tu*(theta/m_theta_dot);
  } catch (const std::invalid_argument& ia) {
    throw std::invalid_argument("Non-orbital Ephemeris Sent to GpAccess: " +
                                m_eph->getName());
  }
}


bool GpAccess::findNextAccess()
{
    // Don't search if current time is past stop time
  if (m_jdStop <= m_jd) {
    return false;
  }

  axs_interval riseset;
  bool in_interval {false};

    // Either already in an access window, or need to locate start time
  if (is_visible(m_jd)) {
    in_interval = true;
    riseset.rise = m_jd;
    // Set elevation too
  } else {
    in_interval = findRise(riseset);
  }

  if (in_interval) {
    findSet(riseset);
    setRiseSetStatus(riseset);
    m_intervals.push_back(riseset);
    return true;
  }

  return false;
}


void GpAccess::findAllAccesses()
{
    // Locate next interval and store
  while (findNextAccess())
  {
      // Set time to search for next access interval
    m_jd = m_intervals.back().set + jd_inc;
  }
}


std::string GpAccess::getGpName() const
{
  return m_gp.getName();
}


std::string GpAccess::getOrbitName() const
{
  return (*m_eph).getName();
}


bool GpAccess::is_visible(const JulianDate& jd)
{
    // Ensure ephemeris past availability is not requested
  if (m_jdStop < m_jd  ||  m_jd < m_jdStart) {
    return false;
  }

  Eigen::Matrix<double, 3, 1> pos = m_eph->getPosition(jd, EphemFrame::ecf);
  if (m_gp.getSinElevation(pos) >= m_xcs.getSineMinEl()) {
    return true;
  }

  return false;
}


bool GpAccess::findRise(axs_interval& axs)
{
  bool found_rise {false};
  while (!(found_rise = is_visible(m_jd))) {
    m_jd += dt_days;
      // If past stop time, access was not found - done
    if (m_jdStop <= m_jd) {
      m_jd = m_jdStop;
      break;
    }
  }

    // Bisection method to refine time:  jd2 < jd1, jd3 is average
  if (found_rise) {
    auto jd1 = m_jd;                        // In access
    auto jd2 = jd1 + -dt_days;              // Behind access
    auto jd3 = jd2 + 0.5*(jd1 - jd2);
    for (int ii=0; (jd1 - jd2)>tol_dt_day  &&  ii<max_itr; ++ii) {
      if (is_visible(jd3)) {
        jd1 = jd3;
      } else {
        jd2 = jd3;
      }
      jd3 = jd2 + 0.5*(jd1 - jd2);
    }
    m_jd = jd3;
  }

  if (found_rise) {
    axs.rise = m_jd;
    m_jd += dt_days;
    return true;
  } else {
    return false;
  }
}


void GpAccess::findSet(axs_interval& axs)
{
  bool found_set {true};
  while (is_visible(m_jd)) {
    m_jd += dt_days;
    if (m_jdStop <= m_jd) {
      m_jd = m_jdStop;
      found_set = false;
      break;
    }
  }

    // Bisection method to refine time:  jd2 < jd1, jd3 is average
  if (found_set) {
    auto jd1 = m_jd;                        // In access
    auto jd2 = jd1 + -dt_days;              // Behind access
    auto jd3 = jd2 + 0.5*(jd1 - jd2);
    for (int ii=0; (jd1 - jd2)>tol_dt_day  &&  ii<max_itr; ++ii) {
      if (!is_visible(jd3)) {
        jd1 = jd3;
      } else {
        jd2 = jd3;
      }
      jd3 = jd2 + 0.5*(jd1 - jd2);
    }
    m_jd = jd3;
  }

  axs.set = m_jd;
  if (found_set) {
    m_jd += dt_days;
  }
}

void GpAccess::setRiseSetStatus(axs_interval& axs)
{
  axs.sinel_rise = m_gp.getSinElevation(m_eph->getPosition(axs.rise,
                                                           EphemFrame::ecf));
  axs.sinel_set = m_gp.getSinElevation(m_eph->getPosition(axs.set,
                                                          EphemFrame::ecf));
}

}
