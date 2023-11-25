/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_gp_access.h>

#include <string>
#include <memory>

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

  Keplerian oe(m_eph->getStateVector(m_jd, EphemFrame::eci));
  m_max_vel = phy_const::earth_equatorial_speed() + oe.getPerigeeSpeed();
}


bool GpAccess::findNextAccess()
{
    // Don't search if current time is past stop time
  if (m_jdStop <= m_jd) {
    return false;
  }

  axs_interval axs;

    // Updates m_jd, axs
  auto findRise = [&]() {
    double dt_days = utl_const::day_per_sec;
    bool found_rise {false};
    while (!(found_rise = is_visible(m_jd))) {
      m_jd += dt_days;
        // If past stop time, access was not found - done
      if (m_jdStop <= m_jd) {
        m_jd = m_jdStop;
        break;
      }
    }

    /*
    Eigen::Matrix<double, 3, 1> pos = m_eph->getPosition(m_jd, EphemFrame::ecf);
    auto dist = (pos - m_gp.getCartesian()).norm();

    m_jd += phy_const::day_per_tu*(dist/m_max_vel);
    axs.rise = m_jd;
    // Set elevation too
    */

    if (found_rise) {
      axs.rise = m_jd;
      return true;
    } else {
      return false;
    }
  };

    // Assumes inside a rise, looking for a set
    // Always a set if currently within an access window
    // Updates m_jd, axs
  auto findSet = [&]() {
    double dt_days = utl_const::day_per_sec;
    while (is_visible(m_jd)) {
      m_jd += dt_days;
      if (m_jdStop <= m_jd) {
        m_jd = m_jdStop;
        break;
      }
    }
    axs.set = m_jd;
  };

  bool in_interval {false};

    // Either already in an access window, or need to locate start time
  if (is_visible(m_jd)) {
    in_interval = true;
    axs.rise = m_jd;
    // Set elevation too
  } else {
    in_interval = findRise();
  }

  if (in_interval) {
    findSet();
    m_intervals.push_back(axs);
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


std::string GpAccess::getGpName() const noexcept
{
  return m_gp.getName();
}


std::string GpAccess::getOrbitName() const
{
  return (*m_eph).getName();
}


bool GpAccess::is_visible(const JulianDate& jd)
{
  Eigen::Matrix<double, 3, 1> pos = m_eph->getPosition(jd, EphemFrame::ecf);

  if (m_gp.getSinElevation(pos) >= m_xcs.getSineMinEl()) {
    return true;
  }

  return false;
}


}
