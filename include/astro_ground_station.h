/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_GROUND_STATION_H
#define ASTRO_GROUND_STATION_H

#include <memory>
#include <string>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_ground_point.h>


namespace eom {

/**
 * Defines a ground station based on an initial position and velocity at
 * a given epoch.  This definition is created specifically for SLR
 * sites in mind, but should be compatible with any ground site where
 * crustial dynamics should be accounted for.  The class is meant for
 * objects moving at rates on the geological time scale such that most
 * users will simply request the ITRF relative/fixed GroundPoint once
 * over the course of a simulation vs.  getPosition(utc) repeatedly.
 *
 * @author  Kurt Motekew
 * @date    2024/09/26
 */
class GroundStation : public Ephemeris {
public:
  /**
   * Initialize with initial state, ECF/ECI services, and optional ID
   *
   * @param  jd0        UTC time of initial state
   * @param  pos0f      ECF Position, DU
   * @param  vel0f      ECF Velocity, DU/TU
   * @param  ecfeciSys  ECF/ECI conversion resource
   * @param  name       Optional unique ID if needed.
   */
  GroundStation(const JulianDate& jd0,
                const Eigen::Matrix<double, 3, 1>& pos0f,
                const Eigen::Matrix<double, 3, 1>& vel0f,
                std::shared_ptr<const EcfEciSys> ecfeciSys,
                const std::string& name = "GroundStation");

  /**
   * @return  Identifier, with defalut value GroundStation unless set
   *          during construction.
   */
  std::string getName() const override
  {
    return m_name;
  }

  /**
   * @return  Orbit epoch
   */
  JulianDate getEpoch() const override
  {
    return m_jdStart;
  }

  /**
   * @return  Earliest time for which ephemeris can be retrieved
   */
  JulianDate getBeginTime() const override
  {
    return m_jdStart;
  }

  /**
   * @return  Latest time for which ephemeris can be retrieved
   */
  JulianDate getEndTime() const override
  {
    return m_jdStop;
  }

  /**
   * Compute position and velocity state vector of the ground site
   *
   * @param  jd     Time of desired state vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian state vector at requested time in the requested
   *          reference frame, DU and DU/TU
   *
   * @throws  out_of_range if the requested time is out of range.  This
   *          would be due to requesting the state vector in the ECI
   *          frame at a time for which ECF/ECI transformation
   *          data is not available vs. any limitations imposed on
   *          modeling the ground site state.
   */
  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate& jd,
                                             EphemFrame frame) const override;

  /**
   * Compute position of the ground station
   *
   * @param  jd     Time of desired position vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian position vector, DU
   *
   * @throws  out_of_range if the requested time is out of range.  This
   *          would be due to requesting the position vector in the ECI
   *          frame at a time for which ECF/ECI transformation
   *          data is not available vs. any limitations imposed on
   *          modeling the ground site state.
   */
  Eigen::Matrix<double, 3, 1> getPosition(const JulianDate& jd,
                                          EphemFrame frame) const override;

  /**
   * Compute position and velocity state vector of the ground site.
   * Note, no exceptions are thrown as internally, the location is
   * computed in earth fixed Cartesian (therefore not requiring any
   * ECF/ECI transformations).
   *
   * @param  jd  Time of desired state vector, UTC
   *
   * @return  Geodetic position
   */
  GroundPoint getGroundPoint(const JulianDate& jd) const;

private:
  JulianDate m_jd0;
  Eigen::Matrix<double, 3, 1> m_pos0f;
  Eigen::Matrix<double, 3, 1> m_vel0f;

  std::shared_ptr<const EcfEciSys> m_ecfeci {nullptr};
  std::string m_name {};

  JulianDate m_jdStart;
  JulianDate m_jdStop;
};


}

#endif
