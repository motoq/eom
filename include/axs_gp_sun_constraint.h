/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef AXS_GP_SUN_CONSTRAINT_H
#define AXS_GP_SUN_CONSTRAINT_H

#include <memory>

#include <utl_constraint_function.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <astro_ground_point.h>
#include <astro_sun_meeus.h>

namespace eom {

/**
 * Determines if sun elevation constraints are satisfied.
 *
 * @author  Kurt Motekew
 * @date    2024/08/20
 */
class GpSunConstraint : public ConstraintFunction<JulianDate> {
public:
  ~GpSunConstraint() = default;
  GpSunConstraint(const GpSunConstraint&) = delete;
  GpSunConstraint& operator=(const GpSunConstraint&) = delete;
  GpSunConstraint(GpSunConstraint&&) = default;
  GpSunConstraint& operator=(GpSunConstraint&&) = default;

  /**
   * Initialize parameters needed to compute constraints based on a
   * ground point and the sun.
   *
   * @param  gp         Ground point location
   * @param  ecfeciSys  ECF/ECI transformation service
   */
  GpSunConstraint(const GroundPoint& gp,
                  std::shared_ptr<const EcfEciSys> ecfeciSys);

  /**
   * Evaluates all enabled sun constraints w.r.t. the ground point.
   *
   * @param  utc  UTC time at which to evaluate if the sun constraint is
   *              satisfied
   *
   * @return  If true, all sun constraints are satisfied
   */
  bool isSatisfied(JulianDate utc) const override;

  /**
   * Enables and sets a maximum elevation of the sun w.r.t. the ground
   * point horizon
   *
   * @param  max_el  Maximum elevation above the plane tangent to the
   *                 ground point, radians
   */
  void setMaxElevation(double max_el);

private:
  GroundPoint m_gp;
  std::shared_ptr<const Ephemeris> m_sunEph;

  bool m_use_max_el {false};
  double m_max_sin_el {0.0};
};


}

#endif
