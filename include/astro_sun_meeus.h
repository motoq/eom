/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_SUN_MEEUS_H
#define ASTRO_SUN_MEEUS_H

#include <string>
#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>


namespace eom {

/**
 * Computes solar coordinates based on Meeus' analytic model that is
 * accurate to approximately 0.01 degrees.  The coordinates are computed
 * relative to a J2000 reference and not transformed to GCRF given the
 * ~20 MAS difference is well in the noise.  Position magnitude is not
 * great, but this is a reasonable approximation for force model
 * perturbations.  Velocity is included to satisfy the Ephemeris
 * interface, computed via simple differencing.
 *
 * Up to 1' of error compared to precision ephemerides during the month
 * of Feb 2023 have been observed with this implementation (note the sun
 * appears to span an arc on the order of 30').  The conversion to the
 * J2000 equinox is valid through 2100.
 *
 * Meeus, Jean, "Astronomical Algorithms", 2nd Ed., Willmann-Bell, Inc.,
 * 1998.  Chapter 25 algorithm converted to Cartesian coordinates.
 *
 * @author  Kurt Motekew
 * @date    2023/01/30
 */
class SunMeeus : public Ephemeris {
public:
  ~SunMeeus() = default;
  SunMeeus(const SunMeeus&) = default;             // copy constructor
  SunMeeus& operator=(const SunMeeus&) = default;  // copy assignment
  SunMeeus(SunMeeus&&) = default;                  // move constructor
  SunMeeus& operator=(SunMeeus&&) = default;       // move assignment

  /**
   * Initialize with ECF/ECI service and optional ID
   *
   * @param  ecfeciSys  ECF/ECI conversion resource
   * @param  name       Optional unique ID if needed.
   */
  SunMeeus(std::shared_ptr<const EcfEciSys> ecfeciSys,
           const std::string& name = "SunMeeus");

  /**
   * @return  Identifier, with defalut value SunMeeus unless set
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
   * Compute position and velocity state vector of the sun
   *
   * @param  jd     Time of desired state vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian state vector at requested time in the requested
   *          reference frame, DU and DU/TU
   *
   * @throws  out_of_range if the requested time is out of range.  This
   *          would be due to a time for which ECF/ECI transformation
   *          data is not available.
   */
  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate& jd,
                                             EphemFrame frame) const override;

  /**
   * Compute position of the sun 
   *
   * @param  jd     Time of desired position vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian position vector, DU
   *
   * @throws  out_of_range if the requested time is out of range.  This
   *          would be due to a time for which ECF/ECI transformation
   *          data is not available.
   */
  Eigen::Matrix<double, 3, 1> getPosition(const JulianDate& jd,
                                          EphemFrame frame) const override;

private:
  std::string m_name {};
  std::shared_ptr<const EcfEciSys> m_ecfeci {nullptr};
  JulianDate m_jdStart;
  JulianDate m_jdStop;
};


}

#endif
