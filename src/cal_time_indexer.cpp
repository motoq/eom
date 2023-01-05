/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <cal_time_indexer.h>

#include <iostream>

#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>

#include <cal_julian_date.h>

#include <utl_const.h>

namespace eom {

// Check second > first
TimeIndexer::TimeIndexer(std::vector<std::pair<JulianDate, JulianDate>> times)
{
  m_times = std::move(times);
  m_dur = m_times.back().second - m_times.front().first;
  
    // Find maximum dt and check second is always > first
  m_dt = m_times.front().second - m_times.front().first;
  for (const auto& interval : m_times) {
    if (interval.second <= interval.first) {
      throw std::invalid_argument("TimeIndexer::TimeIndexer(): Invalid times");
    }
    m_dt = std::max(m_dt, interval.second - interval.first);
  }
  unsigned long m_n {static_cast<unsigned long>(m_dur/m_dt)};

  jd0 = m_times.front().first;
  auto jd = jd0 + m_dt;
  unsigned long ii {0};
  while (jd <= m_times.back().second) {
    bool outside {false};
    while (!outside) {
      if (jd <= m_times[ii].second) {
        outside = true;
        imap.push_back(ii);
      }
      ii++;
    }
    jd += m_dt;
  }

/*
  jd = jd0;
  while (jd < m_times.back().second) {
    auto dur = jd - jd0;
    unsigned long ndx = static_cast<unsigned long>(m_n*(dur/m_dur));
    if (imap.size() <= ndx) {
      ndx = imap.size() - 1UL;
    }
    if (m_times[imap[ndx]].second < jd) {
      std::cout << '\n' << jd.to_str() << "  " <<
                           m_times[imap[ndx]].second.to_str();
    }
    jd += m_dt;
  }
*/

  std::cout << "\nMax Dt (sec): " << m_dt*utl_const::sec_per_day;
}


unsigned long TimeIndexer::getIndex(const JulianDate& jd) const
{
  auto dur = jd - jd0;
  unsigned long ndx = static_cast<unsigned long>(m_n*(dur/m_dur));
  if (imap.size() <= ndx) {
    ndx = imap.size() - 1UL;
  }

  bool found {false};
  while (ndx >= 0  &&  !found) {
    if (m_times[imap[ndx]].first <= jd  &&  jd <= m_times[imap[ndx]].second) {
      found = true;
    } else {
      ndx--;
    }
  }

  return ndx;
}

}
