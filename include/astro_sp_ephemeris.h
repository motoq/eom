/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_SP_EPHEMERIS_H
#define ASTRO_SP_EPHEMERIS_H

#include <string>
#include <vector>
#include <memory>
#include <utility>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <mth_hermite.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <astro_propagator_config.h>

namespace eom {

struct eph_record {
  JulianDate t;
  Eigen::Matrix<double, 3, 1> p;
  Eigen::Matrix<double, 3, 1> v;
  Eigen::Matrix<double, 3, 1> a;

  eph_record(const JulianDate jd,
             const Eigen::Matrix<double, 3, 1>& x,
             const Eigen::Matrix<double, 3, 1>& dx,
             const Eigen::Matrix<double, 3, 1>& ddx) : t(jd),
                                                       p(x),
                                                       v(dx),
                                                       a(ddx)
  {
  }
};


/**
 * SP Ephemeris
 */
class SpEphemeris : public Ephemeris {
public:
  ~SpEphemeris() = default;
  SpEphemeris(const SpEphemeris&) = default;             // copy constructor
  SpEphemeris& operator=(const SpEphemeris&) = default;  // copy assignment
  SpEphemeris(SpEphemeris&&) = default;                  // move constructor
  SpEphemeris& operator=(SpEphemeris&&) = default;       // move assignment

  SpEphemeris(const std::string& name,
              const JulianDate& epoch,
              const Eigen::Matrix<double, 6, 1>& xeci,
              const std::shared_ptr<const EcfEciSys>& ecfeciSys,
              const PropagatorConfig& propCfg);

  /**
   * @return  Empty string
   */
  std::string getName() const override
  {
    return "";
  }

  /**
   * @return  Default Julian Date
   */
  JulianDate getEpoch() const override
  {
    return m_jdStart;
  }

  /**
   * Center of the earth
   *
   * @param  jd     Not used
   * @param  frame  Not used
   *
   * @return  Zero vector
   */
  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate&,
                                             EphemFrame) const override
  {
    return nullState;
  }

private:
  std::string m_name {""};
  JulianDate m_jdStart;
  JulianDate m_jdStop;
  Eigen::Matrix<double, 6, 1> nullState;
  std::shared_ptr<const EcfEciSys> m_ecfeciSys {nullptr};

  std::vector<std::pair<JulianDate, JulianDate>> m_eph_times;
  std::vector<Hermite<double, 3>> m_eph;
};


}

#endif
