/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_BUILD_H
#define ASTRO_BUILD_H

#include <memory>

#include <astro_orbit_def.h>
#include <astro_rel_orbit_def.h>
#include <astro_ephemeris_file.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>

namespace eom {

/**
 * EOM build utility function declarations
 */


/**
 * Creates an ephemeris "service" given an orbit definition.
 * The orbit definition is similar to a declaration, while the result of
 * this function is the full implementation behind the eom::Ephemeris
 * interface.
 *
 * @param  orbitParams  Orbit definition
 * @param  ecfeciSys    Ecf/Eci utility service pointer that will be
 *                      copied into the Ephemeris object.
 *
 * @return  Orbit implementation
 */
std::unique_ptr<Ephemeris>
build_orbit(const OrbitDef& orbitParams,
            const std::shared_ptr<const EcfEciSys>& ecfeciSys);

/**
 * Creates an ephemeris "service" based on a reference orbit and a
 * relative orbit definition.
 *
 * @param  relOrbi    Relative orbit parameters
 * @param  refOrbit   Reference/template orbit
 * @param  refEph     Reference orbit ephemeris
 * @param  ecfeciSys  Ecf/Eci utility service pointer that will be
 *                    copied into the Ephemeris object.
 *
 * @return  Orbit implementation
 */
std::unique_ptr<Ephemeris>
build_orbit(const RelOrbitDef& relOrbit,
            const OrbitDef& refOrbit,
            const std::shared_ptr<eom::Ephemeris>& refEph,
            const std::shared_ptr<const EcfEciSys>& ecfeciSys);

/**
 * Create and ephemeris "service" based on externally generated
 * ephemeris data.
 *
 * @param  efd        Ephemeris file definition
 * @param  startTime  Earliest time for which ephemeris needs to
 *                    be present in file.
 * @param  stopTime   Latest time for which ephemeris needs to be
 *                    present in file.
 */
std::unique_ptr<Ephemeris>
build_ephemeris(const EphemerisFile& efd,
                const JulianDate& startTime,
                const JulianDate& stopTime,
                const std::shared_ptr<const EcfEciSys>& ecfeciSys);
}

#endif
