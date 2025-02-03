/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_FANDG_H
#define ASTRO_FANDG_H

#include <string>
#include <array>
#include <memory>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>


namespace eom {

/**
 * Implementation f and g method of Kepler's problem as
 * presented in BMW, limited to elliptical orbits.
 *
 * Fundamentals of Astrodynamics
 * Bate, Mueller, & White
 */
class FandG : public Ephemeris {
public:
  ~FandG() = default;
  FandG(const FandG&) = default;             // copy constructor
  FandG& operator=(const FandG&) = default;  // copy assignment
  FandG(FandG&&) = default;                  // move constructor
  FandG& operator=(FandG&&) = default;       // move assignment

  /**
   * Initialize Keppler1
   *
   * @param  orbit_name  Name (string identifieer) associated with orbit
   * @param  epoch       Orbit state vector epoch, UTC
   * @param  xeci        Cartesian GCRF state vector, DU and DU/TU
   * @param  ecfeciSys   ECF/ECI conversion resource
   */
  FandG(const std::string& orbit_name,
        const JulianDate& epoch,
        const Eigen::Matrix<double, 6, 1>& xeci,
        std::shared_ptr<const EcfEciSys> ecfeciSys);

  /**
   * @return  Unique ephemeris identifier
   */
  std::string getName() const override;

  /**
   * @return  Orbit epoch
   */
  JulianDate getEpoch() const override;

  /**
   * @return  Earliest time for which ephemeris can be retrieved
   */
  JulianDate getBeginTime() const override;

  /**
   * @return  Latest time for which ephemeris can be retrieved
   */
  JulianDate getEndTime() const override;

  /**
   * Compute state vector given a time
   *
   * @param  jd     Time of desired state vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian state vector at requested time in the requested
   *          reference frame, DU and DU/TU
   *
   * @throws  out_of_range if the requested time is out of range.  This
   *          would be due to a time for which ECF/ECI transformation
   *          data is not available vs. a limitation of the analytic
   *          propagator itself.
   */
  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate& jd,
                                             EphemFrame frame) const override;

  /**
   * @param  jd     Time of desired position vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian position vector, DU
   *
   * @throws  out_of_range if the requested time is out of range.  This
   *          would be due to a time for which ECF/ECI transformation
   *          data is not available vs. a limitation of the analytic
   *          propagator itself.
   */
  Eigen::Matrix<double, 3, 1> getPosition(const JulianDate& jd,
                                          EphemFrame frame) const override;

private:
  std::string m_name;
  JulianDate m_jd0;
  std::shared_ptr<const EcfEciSys> m_ecfeci {nullptr};
  Eigen::Matrix<double, 2, 1> m_r0;
  Eigen::Matrix<double, 2, 1> m_v0;

  double m_r0_mag {0.0};
    // Semimajor axis
  double m_a {0.0};
    // Perifocal to ECI reference frame transformation;
  Eigen::Matrix<double, 3, 3> m_c_ip {Eigen::Matrix<double, 3, 3>::Identity()};
  
};


}

#endif
