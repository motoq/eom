/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_SRP_SPHERICAL_H
#define ASTRO_SRP_SPHERICAL_H

#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_force_model.h>

namespace eom {

/**
 * Spherical solar radiation pressure model.
 *
 * @author  Kurt Motekew
 * @date    2025/06/09
 */
class SrpSpherical : public ForceModel {
public:
  /**
   * Initialize with reflectivity and area over mass ratio
   *
   * @param  cr   Reflectivity.  Valid values range from 0 (translucent)
   *              to 2 (completely reflective) with 1 being black body
   * @param  aom  Area to mass ratio
   */
  SrpSpherical(double cr, double aom, std::unique_ptr<Ephemeris> sun_eph);

  /**
   * Compute SRP acceleration
   *
   * @param  state  ECI state vector at the point for which third
   *                body gravity is to be computed, DU, DU/TU
   *
   * @return  Cartesian acceleration at state vector, ECI, DU/TU^2
   */
  Eigen::Matrix<double, 3, 1>
      getAcceleration(const JulianDate& jd,
                      const Eigen::Matrix<double, 6, 1>& state) override;

private:
  double m_cr {};
  double m_aom {};
  std::unique_ptr<Ephemeris> m_sun_eph {nullptr};
};


}

#endif
