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
   * Default to 0
   */
  Duration()
  {
  }

  /**
   * Initialize with a scalar duration and scale factor to convert to
   * time units.
   *
   * @param  dur    Duration
   * @param  to_tu  Converts dur to TU.  Note utl_units.h has
   *                useful operators, e.g., 5 minutes: Duration(5.0, 1.0_min)
   */
  Duration(double dur, double to_tu) : tu {dur*to_tu},
                                       days {dur*to_tu*phy_const::day_per_tu}
  {
  }

  /**
   * Set with a scalar duration and scale factor to convert to
   * time units.
   *
   * @param  dur    Duration
   * @param  to_tu  Converts dur to TU.  Note utl_units.h has
   *                useful operators, e.g., 5 minutes: Duration(5.0, 1.0_min)
   */
  void set(double dur, double to_tu)
  {
    tu = dur*to_tu;
    days =  tu*phy_const::day_per_tu;
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

  /**
   * Multiply this duration by some value
   *
   * @param  sf  Scale factor
   *
   * @return  Copy of this Duration scaled by sf
   */
  Duration operator*(double sf) const noexcept
  {
    Duration dur {*this};
    dur.tu *= sf;
    dur.days *= sf;

    return dur;
  }
  
private:
  double tu {0.0};
  double days {0.0};
};


}

#endif
