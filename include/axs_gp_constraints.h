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
#include <memory>
#include <vector>

#include <Eigen/Dense>

#include <utl_const.h>
#include <utl_constraint_function.h>
#include <cal_julian_date.h>
#include <astro_ground_point.h>

namespace eom {

/**
 * Ground Point access constraint class.  Azimuth is defined to be
 * clockwise from north, from 0 to 360 degrees.
 *
 * @author  Kurt Motekew 2023/09/15  initial, min elevation only
 * @date    Kurt Motekew 2024/01/07  Added max el and min/max az
 */
class GpConstraints {
public:

  /**
   * Initialize with a minimum elevation angle of zero and no other
   * constraints.
   */
  GpConstraints()
  {
  }

  /**
   * Set minimum elevation angle, measured from the ground point up from
   * the plane tangent to the surface of the central body.
   *
   * @param  min_el  Minimum elevation constraint, radians
   *
   * @throws  invalid_argument if not: -PI/2 <= el <= PI/2
   */
  void setMinEl(double min_el);

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

  /**
   * Set maximum elevation angle, measured from the ground point up from
   * the plane tangent to the surface of the central body.
   *
   * @param  max_el  Maximum elevation constraint, radians
   *
   * @throws  invalid_argument if not: -PI/2 <= el <= PI/2
   */
  void setMaxEl(double max_el);

  /**
   * @return  Maximum elevation, radians
   */
  double getMaxEl() const
  {
    return std::asin(m_sin_max_el);
  }

  /**
   * @return  Sine of maximum elevation constraint
   */
  double getSineMaxEl() const noexcept
  {
    return m_sin_max_el;
  }

  /**
   * Set the minimum and maximum azimuth angles, measured clockwise
   * from north
   *
   * @param  min_az  Minimum azimuth constraint, radians
   * @param  max_az  Maximum azimuth constraint, radians
   *
   * @throws  invalid_argument if not:  0 <= az <= 2*pi
   */
  void setMinMaxAz(double min_az, double max_az);

  /**
   * @return  Minimum azimuth, radians
   */
  double getMinAz() const noexcept
  {
    return m_min_az;
  }

  /**
   * @return  Maximum azimuth, radians
   */
  double getMaxAz() const noexcept
  {
    return m_max_az;
  }

  /**
   * Add a time based constraint to be evaluated during access analysis.
   *
   * @param  constraint  Time dependent constraint function
   */
  void addConstraint(
      std::shared_ptr<const ConstraintFunction<JulianDate>> constraint);

  /**
   * Check for satisfaction of all geometric and time dependent constraints.
   *
   * @param  jd   Time of interest
   * @param  gp   Ground poin
   * @param  pos  Satellite position
   *
   * @return  If visible, return true
   */
  bool isVisible(const JulianDate& jd,
                 const GroundPoint& gp,
                 const Eigen::Matrix<double, 3, 1>& pos) const;
  
private:
  double m_sin_min_el {0.0};
  double m_sin_max_el {1.0};
  double m_min_az {0.0};
  double m_max_az {utl_const::tpi};

  bool m_check_az {false};

  std::vector<
      std::shared_ptr<const ConstraintFunction<JulianDate>>> m_constraints;
};


}

#endif
