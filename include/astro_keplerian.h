/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_KEPLERIAN_H
#define ASTRO_KEPLERIAN_H

#include <ostream>
#include <array>

#include <Eigen/Dense>

#include <astro_orbital_elements.h>

namespace eom {

/**
 * Supports conversions between Cartesian state vectors and Keplerian
 * element sets.
 */
class Keplerian : public OrbitalElements {
public:
  /**
   * Initialize 
   */
  Keplerian(const std::array<double, 6>& kep);

  Keplerian(const Eigen::Matrix<double, 6, 1>& cart);

  /**
   * @return  Keplerian element set
   *            [0] a, semimajor axis, DU
   *            [1] e, eccentricity
   *            [2] i, inclination, radians
   *            [3] o, RAAN, radians
   *            [4] w, argument of perigee, radians
   *            [5] v, true anomaly, radians
   */
  std::array<double, 6> getOrbitalElements() const override
  {
    return m_oe;
  }

  /**
   * 
   */
  Eigen::Matrix<double, 6, 1> getCartesian() const override
  {
    return m_xyz;
  }

  void print(std::ostream& stream) const override;

private:
  std::array<double, 6> m_oe;
  Eigen::Matrix<double, 6, 1> m_xyz;
};


/**
 * Transforms mean anomaly to true anomaly
 *
 * @param  mean_anomaly  Mean anomaly, radians
 * @param  eccentricity  Orbit eccentricity
 *
 * @return  True anomaly, radians
 */
//double trueAnomaly(double mean_anomaly, double eccentricity);

/**
 * Transforms true anomaly to mean anomaly
 *
 * @param  true_anomaly  True anomaly, radians
 * @param  eccentricity  Orbit eccentricity
 *
 * @return  Mean anomaly, radians
 */
//double meanAnomaly(double true_anomaly, double eccentricity);

}

#endif
