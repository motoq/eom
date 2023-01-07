/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef CAL_TIME_INDEXER_H
#define CAL_TIME_INDEXER_H

#include <vector>
#include <utility>

#include <cal_julian_date.h>

namespace eom {

/**
 * Utility
 *
 * @author  Kurt Motekew
 * @date    2023/01/02
 */
class TimeIndexer {
public:
  /**
   * Initialize with a vector of start and stop times that will be moved
   * here and used to form the index mapper.
   *
   * @param
   * @param  dt_eps  Optional endpoint tolerance.
   */
  TimeIndexer(std::vector<std::pair<JulianDate, JulianDate>> times);

  unsigned long getIndex(const JulianDate& jd) const;

private:
  double m_dt {};
  double m_dur {};
  unsigned long m_n {};
  JulianDate m_jd0;

  std::vector<std::pair<JulianDate, JulianDate>> m_times;
  std::vector<unsigned long> imap;
};


}

#endif
