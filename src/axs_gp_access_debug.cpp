/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_gp_access_debug.h>

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
  constexpr double dt_days {8.0*utl_const::day_per_sec};
    // Access convergence
  constexpr double tol_dt_day {0.1*utl_const::day_per_sec};
  constexpr int max_itr {42};
}

namespace eom {

GpAccessDebug::
GpAccessDebug(const JulianDate& jdStart,
              const JulianDate& jdStop,
              const GroundPoint& gp,
              const GpConstraints& xcs,
              std::shared_ptr<const Ephemeris> eph) : m_jdStart {jdStart},
                                                      m_jdStop {jdStop},
                                                      m_gp {gp},
                                                      m_xcs {xcs},
                                                      m_eph { std::move(eph) }
{
  m_jd = jdStart;

  try {
    Keplerian oe(m_eph->getStateVector(m_jd, EphemFrame::eci));
  } catch (const std::invalid_argument& ia) {
    throw std::invalid_argument("Non-orbital Ephemeris in GpAccessDebug: " +
                                m_eph->getName());
  }
}


bool GpAccessDebug::findNextAccess()
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


bool GpAccessDebug::findAllAccesses()
{
  bool found_interval {false};

    // Locate next interval and store
  while (findNextAccess())
  {
    found_interval = true;
      // Set time to search for next access interval
    m_jd = m_intervals.back().set + jd_inc;
  }

  return found_interval;
}


std::string GpAccessDebug::getGpName() const
{
  return m_gp.getName();
}


std::string GpAccessDebug::getOrbitName() const
{
  return (*m_eph).getName();
}


bool GpAccessDebug::is_visible(const JulianDate& jd)
{
    // Ensure ephemeris past availability is not requested
  if (m_jdStop < m_jd  ||  m_jd < m_jdStart) {
    return false;
  }

  Eigen::Matrix<double, 3, 1> pos = m_eph->getPosition(jd, EphemFrame::ecf);

  return m_xcs.isVisible(jd, m_gp, pos);
}


bool GpAccessDebug::findRise(axs_interval& axs)
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


void GpAccessDebug::findSet(axs_interval& axs)
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

void GpAccessDebug::setRiseSetStatus(axs_interval& axs)
{
  axs.rasel_rise = m_gp.getRngAzSinEl(m_eph->getPosition(axs.rise,
                                                         EphemFrame::ecf));
  axs.rasel_set = m_gp.getRngAzSinEl(m_eph->getPosition(axs.set,
                                                        EphemFrame::ecf));
}

}
