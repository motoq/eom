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
 *
 * @author  Kurt Motekew
 * @date    20210923
 */
class Duration {
public:

  /**
   */
  Duration(double dur, double to_tu) : tu {dur*to_tu},
                                       days {dur*to_tu*phy_const::DAY_PER_TU}
  {
  }

  double getDays() const noexcept
  {
    return days;
  }

  double getTU() const noexcept
  {
    return tu;
  }
  
private:
  double tu;
  double days;
};


}

#endif
