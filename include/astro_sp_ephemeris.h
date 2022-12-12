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

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <astro_propagator_config.h>


namespace eom {

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

};


}

#endif
