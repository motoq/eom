/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_VINTI_H
#define ASTRO_VINTI_H

#include <array>
#include <memory>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>


namespace eom {

/*
 * Wrapper class to the Gim J. Der, Albert T. Monuki, & Herb Reynolds
 * Vinti (Vinti6) orbit propagator.  See the Vinti6.c source code.
 *
 * Vinti6.c and Vinti.h source:
 * Supplemental Material
 * https://arc.aiaa.org/doi/suppl/10.2514/4.866487
 *
 * Orbital and Celestial Mechanics
 * Nino L. Bonavito, Gim J. Der, and John P. Vinti
 * AIAA, 1998
 */
class Vinti : public Ephemeris {
public:
  /**
   * Initialize Vinti6
   *
   * @param  epoch      Orbit state vector epoch, UTC
   * @param  xeci       Cartesian GCRF state vector, DU and DU/TU
   * @param  ecfeciSys  ECF/ECI conversion resource
   */
  Vinti(const JulianDate& epoch,
        const Eigen::Matrix<double, 6, 1>& xeci,
        const std::shared_ptr<const EcfEciSys>& ecfeciSys);

  /**
   * Compute state vector given a time
   *
   * @param  jd     Time of desired state vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian state vector at requested time in the requested
   *          reference frame, DU and DU/TU
   */
  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate& jd,
                                             EphemFrame frame) const override;

private:
  std::shared_ptr<const EcfEciSys> ecfeci {nullptr};
  std::array<double, 4> planet = {phy_const::km_per_du,
                                  phy_const::gm_km3_sec2,
                                  phy_const::j2, phy_const::j3};
  JulianDate jd0;
  std::array<double, 6> x0;
};


}

#endif
