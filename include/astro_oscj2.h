/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_OSCJ2_H
#define ASTRO_OSCJ2_H

#include <array>
#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>


namespace eom {

/*
 * Legacy osculating J2 propagator (GENPL option)
 */
class OscJ2 : public Ephemeris {
public:
  ~OscJ2() = default;
  OscJ2(const OscJ2&) = default;
  OscJ2& operator=(const OscJ2&) = default;
  OscJ2(OscJ2&&) = default;
  OscJ2& operator=(OscJ2&&) = default;

  /**
   * Initialize OscJ2
   *
   * @param  epoch      Orbit state vector epoch, UTC
   * @param  xeci       Cartesian GCRF state vector, DU and DU/TU
   * @param  ecfeciSys  ECF/ECI conversion resource
   */
  OscJ2(const JulianDate& epoch,
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
  JulianDate jd0;
  std::array<double, 63> elmn;
};


}

#endif
