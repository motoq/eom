/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef AXS_INTERVAL_H
#define AXS_INTERVAL_H

#include <cal_julian_date.h>

namespace eom {

/**
 * Access start and stop times
 */
struct axs_interval {
  JulianDate start;
  JulianDate stop;
  double sinel_start {};
  double sinel_stop {};
};

}

#endif
