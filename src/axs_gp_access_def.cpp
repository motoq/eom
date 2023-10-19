/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <axs_gp_access_def.h>

#include <memory>

#include <astro_ground_point.h>
#include <astro_ephemeris.h>

namespace eom {

void GpAccessDef::setResources(const std::shared_ptr<GroundPoint>& gp,
                               const std::shared_ptr<eom::Ephemeris>& eph)
{
  m_gp = gp;
  m_eph = eph;
}


}
