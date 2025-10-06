/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_SECULAR_J2_H
#define ASTRO_SECULAR_J2_H

#include <array>
#include <memory>
#include <string>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>

namespace eom {

/*
 * J2 propagator incorporating secular effects on mean motion, RAAN, and
 * argument of perigee.  Ideally, this propagator option should be
 * initialized with a mean element set.  See implementation for details.
 */
class SecularJ2 : public Ephemeris {
public:
  /**
   * Initialize SecularJ2
   *
   * @param  orbit_name  Name (string identifieer) associated with orbit
   * @param  epoch       Orbit state vector epoch, UTC
   * @param  xeci        Cartesian GCRF state vector, DU and DU/TU
   *                     Most appropriate if the state vector was derived
   *                     (converted from) a mean element set
   * @param  ecfeciSys   ECF/ECI conversion resource
   */
  SecularJ2(const std::string& orbit_name,
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
    // Orbit definition epoch
  JulianDate m_jd0;
  std::shared_ptr<const EcfEciSys> m_ecfeci {nullptr};
    // Orbital elements at epoch (true anomaly fast variable)
  std::array<double, 6> m_oe0;
    // Mean anomaly at epoch
  double m_m0 {0.0};

    // Secular variations, mean motion, RAAN_dot, arg perigee rate
  double m_nbar {1.0};
  double m_odot {0.0};
  double m_wdot {0.0};
};


}

#endif
