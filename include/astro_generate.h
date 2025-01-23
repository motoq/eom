/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_GENERATE_H
#define ASTRO_GENERATE_H

#include <memory>
#include <string>

#include <cal_julian_date.h>
#include <cal_duration.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_propagator_config.h>

namespace eom {

/**
 * Creates an ephemeris "service" for a transfer orbit defined
 * by an initial orbit, destination orbit, transfer time, and duration.
 * Transfer time is currently limited to half a rev of the initial
 * orbit (direct transfers vs. multi-rev drift).  Current implementation
 * is primarily meant as a simple illustration of Wiesel's Boundary
 * Value Problem example (see function implementation for reference).
 *
 * @param  orbit_name     Name assigned to new orbit
 * @param  startOrbit     Ephemeris source of initial orbit
 * @param  endOrbit       Ephemeris source of destination orbit
 * @param  xferStartTime  Start time of transfer, UTC
 * @param  xferDur        Duration of transfer.  For the current
 *                        implementation, the transfer time should be
 *                        less than half the period of the shorter
 *                        orbit.
 * @param  propCfg        Configuration to use for transfer orbit.  Note,
 *                        celestial bodies will be limited to sun/moon
 *                        and Meeus' algorithms will be used if external
 *                        ephemerides are requested.
 * @param  ecfeciSys      Ecf/Eci utility service pointer that will be
 *                        copied into the Ephemeris object.
 *
 * @throws  std::invalid_argument  When incompatible arguments are
 *                                 provided, such as ephemeris services
 *                                 that are not valid over required time
 *                                 spans.
 *
 * @return  Transfer orbit
 */
std::unique_ptr<Ephemeris>
generate_xfer_orbit(const std::string& orbit_name,
                    const Ephemeris& startOrbit,
                    const Ephemeris& endOrbit,
                    const JulianDate& xferStartTime,
                    const Duration& xferDur,
                    const PropagatorConfig& propCfg,
                    std::shared_ptr<const EcfEciSys> ecfeciSys);


}

#endif
