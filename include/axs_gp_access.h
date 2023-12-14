/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef AXS_GP_ACCESS_H
#define AXS_GP_ACCESS_H

#include <string>
#include <iterator>
#include <vector>

#include <axs_interval.h>

namespace eom {

/**
 * Interface defining ground point to orbit access analysis
 * functionality.
 *
 * @author  Kurt Motekew
 * @date    2023/12/05
 */
class GpAccess {
public:
  virtual ~GpAccess() = default;
  GpAccess() = default;
  GpAccess(const GpAccess&) = default;
  GpAccess& operator=(const GpAccess&) = delete;
  GpAccess(GpAccess&&) = delete;
  GpAccess& operator=(GpAccess&&) = delete;

  /**
   * Locates and stores the next access interval
   *
   * @return  true if an interval was located.  False if no other
   *          intervals are present over the simulation time.
   */
  virtual bool findNextAccess() = 0;

  /**
   * Computes all remaining access over the entire simulation time and
   * stores them
   */
  virtual void findAllAccesses() = 0;

  /**
   * @return  Name (string identifier) associated with ground point
   */
  virtual std::string getGpName() const = 0;

  /**
   * @return  Name (string identifier) associated with orbit
   */
  virtual std::string getOrbitName() const = 0;

  /**
   * Provides constant iterator access to access interval structures
   *
   * @return  Beginning of range iterator for axs_interval structure
   */
  virtual std::vector<axs_interval>::const_iterator cbegin() const = 0;

  /**
   * Provides constant iterator access to access interval structures.
   *
   * @return  End of range iterator for axs_interval structure
   */
  virtual std::vector<axs_interval>::const_iterator cend() const = 0;

  /**
   * Provides constant iterator access to access interval structures.
   * Note constant implementation.
   *
   * @return  Beginning of range iterator for axs_interval structure
   */
  virtual std::vector<axs_interval>::const_iterator begin() const = 0;

  /**
   * Provides constant iterator access to access interval structures.
   * Note constant implementation.
   *
   * @return  End of range iterator for axs_interval structure
   */
  virtual std::vector<axs_interval>::const_iterator end() const = 0;

};


}

#endif
