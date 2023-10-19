/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef AXS_GP_CONSTRAINTS_H
#define AXS_GP_CONSTRAINTS_H

#include <cmath>

namespace eom {

/**
 * Ground Point access constraint class.
 *
 * @author  Kurt Motekew
 * @date    20230915
 */
class GpConstraints {
public:

  /**
   * Initialize with a minimum elevation angle, measured from the ground
   * point up from the plane tangent to the surface of the central body
   * Default to zero if not supplied.
   *
   * @param  min_el  Minimum elevation constraint, radians
   */
  GpConstraints(double min_el = 0) : m_sin_min_el {std::sin(min_el)}
  {
  }

  /**
   * @param  Update the minimum elevation, radians.
   */
  void setMinEl(double min_el)
  {
    m_sin_min_el = std::sin(min_el);
  }

  /**
   * @return  Minimum elevation, radians
   */
  double getMinEl() const
  {
    return std::asin(m_sin_min_el);
  }

  /**
   * @return  Sine of minimum elevation constraint
   */
  double getSineMinEl() const noexcept
  {
    return m_sin_min_el;
  }

  
private:
  double m_sin_min_el;
};


}

#endif
