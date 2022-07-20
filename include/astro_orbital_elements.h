/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_ORBITAL_ELEMENTS_H
#define ASTRO_ORBITAL_ELEMENTS_H

#include <ostream>
#include <array>

#include <Eigen/Dense>

#include <utl_printable.h>

namespace eom {

/**
 * Interface defining as set of orbital elements with conversions to and
 * from Cartesian state vectors.
 *
 * @author  Kurt Motekew
 * @date    2022/01/01
 */
class OrbitalElements : public Printable {
public:
  virtual ~OrbitalElements() = default;
  OrbitalElements() = default;
  OrbitalElements(const OrbitalElements&) = delete;
  OrbitalElements& operator=(const OrbitalElements&) = delete;
  OrbitalElements(OrbitalElements&&) = delete;
  OrbitalElements& operator=(OrbitalElements&&) = delete;

  /**
   * @return  Orbital elements, various systems of units depending on
   *          the theory, but typically DU and radians.
   */
  virtual std::array<double, 6> getOrbitalElements() const = 0;

  /**
   * @return  Equivalent Cartesian position and velocity state vector,
   *          DU and DU/TU
   */
  virtual Eigen::Matrix<double, 6, 1> getCartesian() const = 0;

  /**
   * @return  Orbit angular momentum, DU^2/TU
   */
  virtual double getAngularMomentum() const = 0;

};


}

#endif
