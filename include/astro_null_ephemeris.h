/*
 * Copyright 2021 Kurt Motekew
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
 * Bad Orbit
 */
class NullEphemeris : public Ephemeris {
public:
  NullEphemeris() = default;
  ~NullEphemeris() = default;
  NullEphemeris(const NullEphemeris&) = default;             // copy constructor
  NullEphemeris& operator=(const NullEphemeris&) = default;  // copy assignment
  NullEphemeris(NullEphemeris&&) = default;                  // move constructor
  NullEphemeris& operator=(NullEphemeris&&) = default;       // move assignment

  /**
   * @return  Unique ephemeris identifier
   */
  std::string getName() const override
  {
    return "";
  }

  /**
   * @return  Orbit epoch
   */
  JulianDate getEpoch() const override
  {
    return jd0;
  }

  /**
   * Compute state vector given a time
   *
   * @param  jd     Time of desired state vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian state vector at requested time in the requested
   *          reference frame, DU and DU/TU
   */
  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate&,
                                             EphemFrame) const override
  {
    return nullState;
  }

private:
  JulianDate jd0;
  Eigen::Matrix<double, 6, 1> nullState = Eigen::Matrix<double, 6, 1>::Zero();
};


}

#endif
