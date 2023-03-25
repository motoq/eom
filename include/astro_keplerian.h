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
 * element sets.  The fast variable is true anomaly.  Mean and eccentric
 * anomaly are available via function calls.  Only elliptical orbits are
 * considered.
 *
 *
 * @author  Kurt Motekew
 * @date    2022/01/01
 * @date    2023/03/25  Added error checking
 */
class Keplerian : public OrbitalElements {
public:
  ~Keplerian() = default;
  Keplerian(const Keplerian&) = default;
  Keplerian& operator=(const Keplerian&) = default;
  Keplerian(Keplerian&&) = default;
  Keplerian& operator=(Keplerian&&) = default;

  /**
   * Initialize with Keplerian orbital elements.
   *
   * @param  oe  Keplarian orbital elements
   *              [0] a, semimajor axis, DU
   *              [1] e, eccentricity
   *              [2] i, inclination, radians
   *              [3] o, RAAN, radians
   *              [4] w, argument of perigee, radians
   *              [5] v, true anomaly, radians
   *
   * @throws  invalid_argument if either the eccentricity or inclination
   *          are too close to zero.  Use a different orbital element
   *          type for zero eccentricity, possibly just a different ECI
   *          reference frame for inclination.  Also thrown if the
   *          perigee radius is less than 1 DU since most gravity models
   *          are not valid below the gravitational scaling factor.
   */
  Keplerian(const std::array<double, 6>& oe);

  /**
   * Initialize with an inertial Cartesian state vector
   *
   * @param  cart  Position and velocity, DU and DU/TU
   *
   * @throws  invalid_argument if the state vector transforms to orbital
   *          elements with the same characteristics as the above
   *          constructor.
   */
  Keplerian(const Eigen::Matrix<double, 6, 1>& cart);

  /**
   * @return  Orbital elements as described in constructor, DU and
   *          radians
   */
  std::array<double, 6> getOrbitalElements() const override
  {
    return m_oe;
  }

  /**
   * @return  Cartesian state vector, DU and DU/TU
   */
  Eigen::Matrix<double, 6, 1> getCartesian() const override
  {
    return m_cart;
  }

  /**
   * @return  Orbit energy, negative for elliptical orbits, DU^2/TU^2
   */
  double getEnergy() const override
  {
    return m_sme;
  }

  /**
   * @return  Orbit angular momentum, DU^2/TU
   */
  double getAngularMomentum() const override
  {
    return m_hmag;
  }

  /**
   * Output orbital elements and state vector to the supplied stream.
   * Units of kilometers, seconds, and degrees are used.
   *
   * @param  stream  The stream to output the scenario setting to
   */
  void print(std::ostream& stream) const override;

  /*
   * @return the eccentric anomaly, radians
   */
  double getEccentricAnomaly() const;

  /*
   * @return the mean anomaly, radians
   */
  double getMeanAnomaly() const;

  /*
   * @return the orbit period, TU
   */
  double getPeriod() const;

  /**
   * @return  Perigee radial distance, DU
   */
  double getPerigeeRadius() const
  {
    return m_rp;
  }

  /*
   * Updates the true anomaly given the input mean anomaly.  Only the
   * true anomaly is modified.
   *
   * @param  ma  Mean anomaly, radians
   *
   * @throws  NonconvergenceException if the conversion fails to
   *          converge to the required tolerance within the required
   *          number of iterations.
   */
  void setWithMeanAnomaly(double ma);

private:
  void set(const std::array<double, 6>& oe);

  double m_sme {};
  double m_hmag {};
  double m_rp {};
  std::array<double, 6> m_oe;
  Eigen::Matrix<double, 6, 1> m_cart;
};


}

#endif
