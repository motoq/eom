/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_gp_access.h>

#include <memory>

#include <astro_ground_point.h>
#include <astro_ephemeris.h>
#include <axs_gp_constraints.h>

namespace eom {


GpAccess::GpAccess(const GroundPoint& gp,
                   const std::shared_ptr<Ephemeris>& eph,
                   const GpConstraints& xcs)
{
  m_gp = std::make_unique<GroundPoint>(gp);
  m_eph = eph;
  m_xcs = xcs;
}


}
