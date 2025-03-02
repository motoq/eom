/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_HERMITE1_TC_EPH_H
#define ASTRO_HERMITE1_TC_EPH_H

#include <string>
#include <vector>
#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <mth_hermite1.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>

namespace eom {

/**
 * Hermite interpolation of a central body w.r.t. the earth and a target
 * body w.r.t. the central body is employed.  Ttwo position and velocity
 * pairs are employed to generate state vector from sets of inertial
 * ephemerides given a time.  See Hermite1Eph for more info.
 *
 * @author  Kurt Motekew  2023/12/31
 */
class Hermite1TcEph : public Ephemeris {
public:
  /**
   * Initialize central body and target ephemerides.
   *
   * @param  name         Unique identifier for the resulting ephemeris.
   * @param  target_recs  Position and velocity records to be used for
   *                      the target ephemeris.  Position is w.r.t.
   *                      the center_recs ephemeris.  Inertial, DU and
   *                      DU/TU
   * @param  center_recs  Position and velocity records to be used for
   *                      the center ephemeris.  ECI, DU, and DU/TU
   * @param  jdStart      Start time for which ephemeris must be available
   * @param  jdStop       End time for which ephemeris must be available
   * @param  ecfeciSys    ECF/ECI conversion resource
   *
   * @throws  runtime_error for initialization error
   */
  Hermite1TcEph(const std::string& name,
                const std::vector<state_vector_rec>& target_recs,
                const std::vector<state_vector_rec>& center_recs,
                const JulianDate& jdStart,
                const JulianDate& jdStop,
                std::shared_ptr<const EcfEciSys> ecfeciSys);

  /**
   * @return  Unique ephemeris identifier
   */
  std::string getName() const override
  {
    return m_name;
  }

  /**
   * @return  Get initializing state vector time
   */
  JulianDate getEpoch() const override
  {
    return m_jdEpoch;
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
   * Interpolate state vector from stored ephemeris for given time
   *
   * @param  jd     Time of desired state vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian state vector at requested time in the requested
   *          reference frame, DU and DU/TU
   *
   * @throws  out_of_range if the requested time is out of range
   */
  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate&,
                                             EphemFrame frame) const override;

  /**
   * @param  jd     Time of desired position vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian position vector, DU
   *
   * @throws  out_of_range if the requested time is out of range
   */
  Eigen::Matrix<double, 3, 1> getPosition(const JulianDate& jd,
                                          EphemFrame frame) const override;

private:
  std::string m_name;
  JulianDate m_jdStart;
  JulianDate m_jdStop;
  JulianDate m_jdEpoch;
  std::shared_ptr<const EcfEciSys> m_ecfeciSys {nullptr};

  std::unique_ptr<Ephemeris> m_targetEph {nullptr};
  std::unique_ptr<Ephemeris> m_centerEph {nullptr};
};


}

#endif
