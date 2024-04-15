/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef OBS_RNG_AZ_EL_H
#define OBS_RNG_AZ_EL_H

namespace eom {

/**
 * Range, azimuth, and sine of elevation
 */
template<typename T>
struct rng_az_sinel {
  T range {static_cast<T>(0.0)};
  T azimuth {static_cast<T>(0.0)};
  T sinel {static_cast<T>(0.0)};
};

}

#endif
