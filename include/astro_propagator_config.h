/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_PROPAGATOR_CONFIG_H
#define ASTRO_PROPAGATOR_CONFIG_H

namespace eom {

/**
 * Supported propagators
 */
enum class PropagatorType {
  Kepler1,                        ///< Gim J. Der version
  Vinti6,                         ///< Gim J. Der version
  VintiJ2,                        ///< Gim J. Der version with J3 = 0
  OscJ2
};


/**
 * Contains propagator configuration parameters
 */
class PropagatorConfig {
public:
  /**
   * Initilize with the Kepler1 version of the two-body propagator
   */
  PropagatorConfig() : prop {PropagatorType::Kepler1}
  {
  }

  /**
   * Initialize with only the propagator type.  Most general
   * perturbation propagators should use this version.  Special
   * perturbation propagators will be initialized with default values.
   */
  PropagatorConfig(PropagatorType prop_type) : prop {prop_type}
  {
  }

  /*
   * @return  The propagator type
   */
  PropagatorType getPropagatorType() const noexcept { return prop; }

private:
  PropagatorType prop;
};


}

#endif
