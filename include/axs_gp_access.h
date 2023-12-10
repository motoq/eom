/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef AXS_GP_ACCESS_H
#define AXS_GP_ACCESS_H

#include <string>
#include <memory>
#include <vector>
#include <algorithm>

#include <utl_const.h>
#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ground_point.h>
#include <astro_ephemeris.h>
#include <axs_gp_constraints.h>
#include <axs_interval.h>
#include <axs_gp_visibility.h>

namespace {
  constexpr double k {0.99547*utl_const::rad_per_deg};
  constexpr double c {-0.1481*phy_const::tu_per_min};
  constexpr double lb {8.0*phy_const::tu_per_sec};
  constexpr double ub {2.0*phy_const::tu_per_min};
  /**
   * Given the angular velocity of the satellite w.r.t. the center of
   * the earth (rate of change in true anomaly), compute the time
   * increment used when searching for access interval bounds.
   *
   * @param  theta_dot  True anomaly rate, rad/TU
   *
   * @return  Time increment, TU
   */
  constexpr double search_stepsize(double theta_dot) {
    return std::min(std::max(k/theta_dot + c, lb), ub);
  }
}

namespace eom {


/**
 * Performs access analysis between a ground point and ephemeris
 * resource.  The ephemeris resource is assumed to be a valid orbit.
 *
 * @author  Kurt Motekew
 * @date    20231027
 */
class GpAccess : public GpVisibility {
public:
  /**
   * Initialize with ground point, ephemeris, and static constraints
   *
   * @param  gp   Ground point definition
   * @param  eph  Ephemeris source
   * @param  xcs  Access constraints
   *
   * @throws  invalid_argument if not orbital ephemeris
   */
  GpAccess(const JulianDate& jdStart,
           const JulianDate& jdStop,
           const GroundPoint& gp,
           const GpConstraints& xcs,
           std::shared_ptr<const Ephemeris> eph);
  ~GpAccess() = default;
  GpAccess(const GpAccess&) = default;
  GpAccess& operator=(const GpAccess&) = default;
  GpAccess(GpAccess&&) = default;
  GpAccess& operator=(GpAccess&&) = default;

  /**
   * Computes access analysis over entire duration.
   */
  bool findNextAccess() override;

  /**
   * Computes access analysis over entire duration.
   */
  void findAllAccesses() override;

  /**
   * @return  Name (string identifier) associated with ground point
   */
  std::string getGpName() const override;

  /**
   * @return  Name (string identifier) associated with orbit
   */
  std::string getOrbitName() const override;

  /**
   * Provides constant iterator access to access interval structures
   *
   * @return  Beginning of range iterator for axs_interval structure
   */
  std::vector<axs_interval>::const_iterator cbegin() const override {
    return m_intervals.cbegin();
  }

  /**
   * Provides constant iterator access to access interval structures.
   *
   * @return  End of range iterator for axs_interval structure
   */
  std::vector<axs_interval>::const_iterator cend() const override {
    return m_intervals.cend();
  }

  /**
   * Provides constant iterator access to access interval structures.
   * Note constant implementation.
   *
   * @return  Beginning of range iterator for axs_interval structure
   */
  std::vector<axs_interval>::const_iterator begin() const override {
    return m_intervals.cbegin();
  }

  /**
   * Provides constant iterator access to access interval structures.
   * Note constant implementation.
   *
   * @return  End of range iterator for axs_interval structure
   */
  std::vector<axs_interval>::const_iterator end() const override {
    return m_intervals.cend();
  }

private:
  /*
   * Given the time of interest, evaluates if access is satisfied based
   * on stored constraints.  Returns false if requested time is outside
   * the open interval defined by m_jdStart and m_jdStop.
   *
   * @param  jd           Time to evaluate if access constraints are met
   * @param  new_dt_days  Time increment to use for next time increment
   *                      in search.  If nullptr, then ignored.
   */
  bool is_visible(const JulianDate& jd, double* new_dt_days) const; 

  /*
   * Locate the start of an access window based on the assumption that
   * m_jd is currently before an access window.  m_jd will be updated.
   *
   * @param  axs  Access interval to be updated with rise time related
   *              values if access constraints are satisfied.
   *
   * @return  If the start of an access window was found before
   *          m_jdStop, return true.
   *
   * Requirement:  Ensure currently outside of an access window before
   *               calling (is_visible(m_jd) == false;).
   */
  bool findRise(axs_interval& axs);

  /*
   * Locate the end of an access window based on the assumption that
   * m_jd is currently within an access window.  If m_jd exceeds
   * m_jdStop, then the set time will be set to m_jdStop.  A set time
   * will always exist based on the above assumption and process.
   *
   * @param  axs  Access interval to be updated with set time related
   *              values.
   *
   * Requirement:  Ensure currently inside of an access window before
   *               calling (is_visible(m_jd) == true;).
   */
  void findSet(axs_interval& axs);

  /*
   * Set geometry constraints in axs_interval for located rise and set
   * times.
   *
   * @param  axs  Access interval to be updated with rise and set time
   *              geometry values.
   *
   * Requirement:  Valid rise and set times have been set in axs
   */
  void setRiseSetStatus(axs_interval& axs);

  JulianDate m_jdStart;
  JulianDate m_jdStop;
  GroundPoint m_gp;
  GpConstraints m_xcs;
  std::shared_ptr<const Ephemeris> m_eph;

  JulianDate m_jd;
  double m_dt_days_p {20.0*utl_const::day_per_sec};
  double m_dt_days_a {20.0*utl_const::day_per_sec};
  double m_rp {1.0};
  double m_ra {1.0};
  double m_ecc {0.0};

  std::vector<axs_interval> m_intervals;
};


}

#endif
