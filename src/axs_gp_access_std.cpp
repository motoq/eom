/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_gp_access_std.h>

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
    // Eccentricity switchpoints for interpolation
  constexpr double ecc_t {0.07};
  constexpr double ecc_t2 {0.3};
    // Access convergence
  constexpr double tol_dt_day {0.1*utl_const::day_per_sec};
  constexpr int max_itr {42};
}

namespace eom {

GpAccessStd::GpAccessStd(const JulianDate& jdStart,
                         const JulianDate& jdStop,
                         const GroundPoint& gp,
                         const GpConstraints& xcs,
                         std::shared_ptr<const Ephemeris> eph) :
                             m_jdStart {jdStart},
                             m_jdStop {jdStop},
                             m_gp {gp},
                             m_xcs {xcs},
                             m_eph {std::move(eph)}
{
  m_jd = jdStart;

  try {
    Keplerian oe(m_eph->getStateVector(m_jd, EphemFrame::eci));
    m_rp = oe.getPerigeeRadius();
    m_ra = oe.getApogeeRadius();
    double ecc {oe.getEccentricity()};
    if (ecc > ecc_t2) {
      m_exp_dt = true;
    } else if (ecc > ecc_t) {
      m_linear_dt = true;
    }
      // Perigee
    double theta_dot_p {oe.getPerigeeSpeed()/m_rp};
    m_dt_days_p = phy_const::day_per_tu*search_stepsize(theta_dot_p);
      // Apogee
    double theta_dot_a {oe.getApogeeSpeed()/m_ra};
    m_dt_days_a = phy_const::day_per_tu*search_stepsize(theta_dot_a);
    if (m_exp_dt) {
        // Semilatus Rectum, get flight path angle for v = 90
      double fpa_slr {std::atan2(ecc/sqrt(1.0 + ecc*ecc),
                                 1.0/sqrt(1.0 + ecc*ecc))};
      double a {0.5*(m_ra + m_rp)};
      double r_slr {oe.getSemilatusRectum()};
      double v_slr {std::sqrt(phy_const::gm*(2.0/r_slr - 1.0/a))};
      double v_slr_t {v_slr*std::cos(fpa_slr)};
      double theta_dot_slr {v_slr_t/r_slr};
      m_dt_days_slr = phy_const::day_per_tu*search_stepsize(theta_dot_slr);
      Eigen::Matrix<double, 3, 1> y = {m_dt_days_p, m_dt_days_a, m_dt_days_slr};
      Eigen::Matrix<double, 3, 3> ap;
      ap(0,0) = 1;
      ap(1,0) = 1;
      ap(2,0) = 1;
      ap(0,1) = m_rp;
      ap(1,1) = r_slr;
      ap(2,1) = m_ra;
      ap(0,2) = std::pow(10, m_rp);
      ap(1,2) = std::pow(10, r_slr);
      ap(2,2) = std::pow(10, m_ra);
      Eigen::ColPivHouseholderQR<Eigen::Matrix<double, 3, 3>> qr(ap);
      Eigen::Matrix<double, 3, 1> phat = qr.solve(y);
      m_exp_a0 = phat(0);
      m_exp_a1 = phat(1);
      m_exp_a10 = phat(2);
    }
  } catch (const std::invalid_argument& ia) {
    throw std::invalid_argument("Non-orbital Ephemeris Sent to GpAccessStd: " +
                                m_eph->getName());
  }
}


bool GpAccessStd::findNextAccess()
{
    // Don't search if current time is past stop time
  if (m_jdStop <= m_jd) {
    return false;
  }

  axs_interval riseset;
  bool in_interval {false};

    // Either already in an access window, or need to locate start time
  if (is_visible(m_jd, nullptr)) {
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


void GpAccessStd::findAllAccesses()
{
    // Locate next interval and store
  while (findNextAccess())
  {
      // Set time to search for next access interval
    m_jd = m_intervals.back().set + jd_inc;
  }
}


std::string GpAccessStd::getGpName() const
{
  return m_gp.getName();
}


std::string GpAccessStd::getOrbitName() const
{
  return (*m_eph).getName();
}


bool GpAccessStd::is_visible(const JulianDate& jd, double* new_dt_days) const
{
    // Ensure ephemeris past availability is not requested
  if (m_jdStop < m_jd  ||  m_jd < m_jdStart) {
    return false;
  }

    // Update time step size for access boundary search
  Eigen::Matrix<double, 3, 1> pos = m_eph->getPosition(jd, EphemFrame::ecf);
  if (m_linear_dt  &&  new_dt_days != nullptr) {
    double r {pos.norm()};
    double factor {(r - m_rp)/(m_ra - m_rp)};
      // Could be negative if orbit degrades since epoch
    if (factor > 0) {
      *new_dt_days = m_dt_days_p + factor*(m_dt_days_a - m_dt_days_p);
    }
  } else if (m_exp_dt  &&  new_dt_days != nullptr) {
    double r {pos.norm()};
    *new_dt_days = std::min(std::max(m_exp_a0 +
                                     m_exp_a1*r +
                                     m_exp_a10*r*std::pow(10, r), lb), ub);
  }

  if (m_xcs.onlyUseMinEl()  &&
      m_gp.getSinElevation(pos) >= m_xcs.getSineMinEl()) {
    return true;
  } else {
    return m_xcs.isVisible(jd, m_gp, pos);
  }

  return false;
}


bool GpAccessStd::findRise(axs_interval& axs)
{
  double dt_days {m_dt_days_p};
  bool found_rise {false};
  while (!(found_rise = is_visible(m_jd, &dt_days))) {
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
      if (is_visible(jd3, nullptr)) {
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


void GpAccessStd::findSet(axs_interval& axs)
{
  double dt_days {m_dt_days_p};
  bool found_set {true};
  while (is_visible(m_jd, &dt_days)) {
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
      if (!is_visible(jd3, nullptr)) {
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

void GpAccessStd::setRiseSetStatus(axs_interval& axs)
{
  axs.sinel_rise = m_gp.getSinElevation(m_eph->getPosition(axs.rise,
                                                           EphemFrame::ecf));
  axs.sinel_set = m_gp.getSinElevation(m_eph->getPosition(axs.set,
                                                          EphemFrame::ecf));
}

}
