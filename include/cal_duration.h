/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef CAL_DURATION_H
#define CAL_DURATION_H

#include <phy_const.h>

namespace eom {

/**
 * Represents a duration in time that will provides time in appropriate
 * units for functions that need to increment/decrement time.  E.G.,
 * when added to a JulianDate, the getDays() method can be used.  Most
 * other astrodynamics uses would request time units.
 *
 * @author  Kurt Motekew
 * @date    20210923
 */
class Duration {
public:
  /**
   * Default to 1 TU
   */
  Duration() : tu {1.0}, days {phy_const::day_per_tu}
  {
  }

  /**
   * Initialize with a scalar duration and scale factor to convert to
   * time units.
   *
   * @param  dur    Duration
   * @param  to_tu  Converts dur to TU.  Note utl_units.h has
   *                useful operators (e.g., 5 minutes as 5.0_min).
   */
  Duration(double dur, double to_tu) : tu {dur*to_tu},
                                       days {dur*to_tu*phy_const::day_per_tu}
  {
  }

  /**
   * @return  Duration in days.
   */
  double getDays() const noexcept
  {
    return days;
  }

  /**
   * @return  Duration in TU, the standard computational time units
   *          established by the application.
   */
  double getTu() const noexcept
  {
    return tu;
  }
  
private:
  double tu;
  double days;
};


}

#endif
