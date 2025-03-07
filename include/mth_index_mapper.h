/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef MTH_INDEX_MAPPER_H
#define MTH_INDEX_MAPPER_H

#include <vector>
#include <utility>
#include <algorithm>
#include <stdexcept>

namespace eom {

/**
 * Given a list of contiguous blocks (e.g., time intervals), creates a
 * mapping scheme where a value covered by the blocks can be mapped to
 * the specific block containing the value.  For example, if
 * interpolating a set of values over unevenly spaced time values, this
 * utility can be initialized with pairs of time intervals and then
 * called upon to locate the index to the interval.
 *
 * Currently base on mapping using the largest block size.  Efficiency is
 * lost if many very small intervals occur with large intervals.  An
 * option to use the smallest block size may be added for these cases.
 *
 * @tparam  T  Type that can be used to mark the beginning and end
 *             of an interval.  It must support {-, <, <=} operations
 *             with subtraction returning a value that can be implicitly
 *             cast to a double.
 *
 * @author  Kurt Motekew
 * @date    2023/01/07
 */
template<typename T>
class IndexMapper {
public:
  /**
   * Initialize with a set of blocks
   *
   * @param  blocks  Contiguous intervals in increasing order.  Each
   *                 interval also represents a value where the first
   *                 boundary is < the second.  Possession of this
   *                 vector is taken via a move operation.
   *
   * @throws  invalid_argument if the second value defining an interval
   *          is less than or equal to the first.
   */
  IndexMapper(std::vector<std::pair<T, T>> blocks);

  /**
   * @param  Value for which a the index to the interval containing that
   *         value is to be found.
   *
   * @return  The index of the interval containing the input value.
   *
   * @throws  out_of_range if value is not covered
   */
  unsigned long getIndex(const T& val) const;

private:
  double m_bsize {};              // Maximum block size
  double m_range {};              // Size of all intervals
  unsigned long m_n {};
  T m_val0;

  std::vector<std::pair<T, T>> m_blocks;
  std::vector<unsigned long> imap;
};


template<typename T>
IndexMapper<T>::IndexMapper(std::vector<std::pair<T, T>> blocks)
{
  m_blocks = std::move(blocks);
  m_range = m_blocks.back().second - m_blocks.front().first;

    // Find (maximum) block size and check second is always > first
  m_bsize = m_blocks.front().second - m_blocks.front().first;
  for (const auto& interval : m_blocks) {
    if (interval.second <= interval.first) {
      throw std::invalid_argument("IndexMapper::IndexMapper(): Invalid blocks");
    }
    m_bsize = std::max(m_bsize, interval.second - interval.first);
  }
  m_n = static_cast<unsigned long>(m_range/m_bsize);

    // Ensure there are not gaps
  for (unsigned int ii=1; ii<m_blocks.size(); ++ii) {
    if (m_blocks[ii-1].second <m_blocks[ii].first) {
      throw std::invalid_argument(
          "IndexMapper::IndexMapper(): Disjoint blocks");
    }
  }

  m_val0 = m_blocks.front().first;
  auto val = m_val0 + m_bsize;
  unsigned long ii {0};
  while (val <= m_blocks.back().second) {
    bool outside {false};
    while (!outside) {
      if (val <= m_blocks[ii].second) {
        outside = true;
        imap.push_back(ii);
      }
      ii++;
    }
    val += m_bsize;
  }
}
  

template<typename T>
unsigned long IndexMapper<T>::getIndex(const T& val) const
{
    // Index mapping is evenly spaced - get approximate mapping index based
    // on percentage of the total number of values
  auto range = val - m_val0;
  unsigned long ndx = static_cast<unsigned long>(m_n*(range/m_range));
  if (imap.size() <= ndx) {
    ndx = imap.back();
  } else {
    ndx = imap[ndx];
  }

    // Now retrieve interval index based on initial estimate
    // The search should usually be backward - check for forward search
  if (m_blocks.size() <= ndx) {
    ndx = m_blocks.size() - 1UL;
  }
  auto ndx0 = ndx;
  bool found {false};
  while (!found) {
    if (m_blocks[ndx].second < val) {
      break;
    }
    if (m_blocks[ndx].first <= val  &&  val <= m_blocks[ndx].second) {
      found = true;
    } else if (ndx == 0) {
      break;
    } else {
      ndx--;
    }
  }
    // Forward search
  if (!found) {
    ndx = ndx0;
    while (ndx < m_blocks.size()  &&  !found) {
      if (m_blocks[ndx].first <= val  &&  val <= m_blocks[ndx].second) {
        found = true;
      } else {
        ndx++;
      }
    }
  }

  if (!found) {
    throw std::out_of_range("IndexMapper::getIndex() - bad value");
  }

  return ndx;
}


}

#endif
