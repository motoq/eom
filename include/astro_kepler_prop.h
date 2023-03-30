/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_KEPLER_PROP_H
#define ASTRO_KEPLER_PROP_H

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
 * Experiement in progress...
 *
 * Driver class to the Gim J. Der & Herbert B. Reynolds Kepler 2-body
 * analytic propagator.  See the Vinti6.c source code.
 *
 * Vinti6.c and Vinti.h source:
 * Supplemental Material
 * https://arc.aiaa.org/doi/suppl/10.2514/4.866487
 *
 * Orbital and Celestial Mechanics
 * Nino L. Bonavito, Gim J. Der, and John P. Vinti
 * AIAA, 1998
 */
class KeplerProp : public Ephemeris {
public:
  ~KeplerProp() = default;
  KeplerProp(const KeplerProp&) = default;             // copy constructor
  KeplerProp& operator=(const KeplerProp&) = default;  // copy assignment
  KeplerProp(KeplerProp&&) = default;                  // move constructor
  KeplerProp& operator=(KeplerProp&&) = default;       // move assignment

  /**
   * Initialize Keppler1
   *
   * @param  orbit_name  Name (string identifieer) associated with orbit
   * @param  epoch       Orbit state vector epoch, UTC
   * @param  xeci        Cartesian GCRF state vector, DU and DU/TU
   * @param  ecfeciSys   ECF/ECI conversion resource
   */
  KeplerProp(const std::string& orbit_name,
             const JulianDate& epoch,
             const Eigen::Matrix<double, 6, 1>& xeci,
             const std::shared_ptr<const EcfEciSys>& ecfeciSys);

  /**
   * @return  Unique ephemeris identifier
   */
  std::string getName() const override
  {
    return name;
  }

  /**
   * @return  Orbit epoch
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
    return ecfeci->getBeginTime();
  }

  /**
   * @return  Latest time for which ephemeris can be retrieved
   */
  JulianDate getEndTime() const override
  {
    return ecfeci->getEndTime();
  }

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
  std::string name {""};
  std::shared_ptr<const EcfEciSys> ecfeci {nullptr};
  JulianDate jd0;
  std::array<double, 6> x0;
};


}

#endif
