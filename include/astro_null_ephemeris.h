/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_NULL_EPHEMERIS_H
#define ASTRO_NULL_EPHEMERIS_H

#include <string>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>


namespace eom {

/**
 * Null Ephemeris object
 */
class NullEphemeris : public Ephemeris {
public:
  NullEphemeris() = default;

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
    return jd0;
  }

  /**
   * @return  Earliest time for which ephemeris can be retrieved
   */
  JulianDate getBeginTime() const override
  {
    return jd0;
  }

  /**
   * @return  Latest time for which ephemeris can be retrieved
   */
  JulianDate getEndTime() const override
  {
    return jd0;
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

  /**
   * @param  jd     Not used
   * @param  frame  Not used
   *
   * @return  Zero vector
   */
  Eigen::Matrix<double, 3, 1> getPosition(const JulianDate&,
                                          EphemFrame) const override
  {
    return nullState.block<3,1>(0,0);
  }

private:
  JulianDate jd0;
  Eigen::Matrix<double, 6, 1> nullState = Eigen::Matrix<double, 6, 1>::Zero();
};


}

#endif
