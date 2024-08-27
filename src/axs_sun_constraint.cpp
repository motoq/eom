/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_sun_constraint.h>

#include <memory>
#include <utility>

#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>

namespace eom {

SunConstraint::SunConstraint(std::shared_ptr<const Ephemeris> satEph,
                             std::shared_ptr<const Ephemeris> sunEph,
                             std::shared_ptr<const EcfEciSys> ecfeciSys) :
                             m_satEph {std::move(satEph)},
                             m_sunEph {std::move(sunEph)},
                             m_ecfeciSys {std::move(ecfeciSys)}
{
}


}
