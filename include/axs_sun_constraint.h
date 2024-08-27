/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef AXS_SUN_CONSTRAINT_H
#define AXS_SUN_CONSTRAINT_H

#include <memory>

#include <utl_constraint_function.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>

namespace eom {

/**
 * Determines if sun elevation constraints are satisfied.
 *
 * @author  Kurt Motekew
 * @date    2024/08/20
 */
class SunConstraint : public ConstraintFunction<JulianDate> {
public:
  ~SunConstraint() = default;
  SunConstraint(const SunConstraint&) = delete;
  SunConstraint& operator=(const SunConstraint&) = delete;
  SunConstraint(SunConstraint&&) = default;
  SunConstraint& operator=(SunConstraint&&) = default;

  SunConstraint(std::shared_ptr<const Ephemeris> satEph,
                std::shared_ptr<const Ephemeris> sunEph,
                std::shared_ptr<const EcfEciSys> ecfeciSys);

  /**
   * @return  Time associated with current state vector and derivative, UTC
   */
  bool isSatisfied(JulianDate) const override
  {
    return true;
  }

private:
  std::shared_ptr<const Ephemeris> m_satEph;
  std::shared_ptr<const Ephemeris> m_sunEph;
  std::shared_ptr<const EcfEciSys> m_ecfeciSys;
};


}

#endif
