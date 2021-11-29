/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef CAL_LEAP_SECONDS_H
#define CAL_LEAP_SECONDS_H

#include <cal_const.h>
#include <cal_julian_date.h>

namespace eom {

/**
 * Singleton class containing a manually set leap second.  The default
 * number of leap seconds is set to zero.  The proper number should be
 * set at program startup.
 *
 * Example use:
 *   eom::LeapSeconds& ls = eom::LeapSeconds::getInstance();
 *   ls.setTai_Utc(100.0);
 *   std::cout << ls.getTai_Utc();
 *
 * @author  Kurt Motekew
 * @date    202109  Initial
 */
class LeapSeconds
{
public:
  LeapSeconds(LeapSeconds const&) = delete;
  LeapSeconds(LeapSeconds const&&) = delete;
  void operator=(LeapSeconds const&)  = delete;
  void operator=(LeapSeconds const&&)  = delete;

  /**
   * @return  The static reference to the leap second utility
   */
  static LeapSeconds& getInstance()
  {
    static LeapSeconds instance;
    return instance;
  }

  /**
   * @param  dat  Set TAI - UTC, the number of leap seconds,
   *              in seconds
   */
  void setTai_Utc(double dat) noexcept
  {
    taimutc = dat;
  }

  /**
   * @return  TAI - UTC, the number of LeapSecondss, seconds
   */
  double getTai_Utc() const noexcept
  {
    return taimutc;
  }

  JulianDate utc2tt(const JulianDate& utc) const noexcept {
    return utc + (taimutc + cal_const::ttmtai)*cal_const::day_per_sec;
  }

private:
    // Disallow construction
  LeapSeconds() {}

    // TAI - UTC, seconds
  double taimutc {0.0};
};


}

#endif
