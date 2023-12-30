/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_BUILD_H
#define ASTRO_BUILD_H

#include <string>
#include <vector>
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
 * Create an ephemeris "service" based on externally generated
 * ephemeris data.
 *
 * @param  efd        Ephemeris file definition
 * @param  startTime  Earliest time for which ephemeris needs to
 *                    be present in file.
 * @param  stopTime   Latest time for which ephemeris needs to be
 *                    present in file.
 * @param  ecfeciSys  Ecf/Eci utility service pointer that will be
 *                    copied into the Ephemeris object.
 */
std::unique_ptr<Ephemeris>
build_ephemeris(const EphemerisFile& efd,
                const JulianDate& startTime,
                const JulianDate& stopTime,
                const std::shared_ptr<const EcfEciSys>& ecfeciSys);

/**
 * Create an ephemeris "service" for a celestial object based on
 * an .emb (eom binary/unformatted) ephemeris file.
 *
 * @param  name_prefix  Name of the celestial body (Moon, Sun, Mercury,
 *                      Venus, Mars, Jupiter, Saturn, Uranus, Neptune,
 *                      and yes, the planet Pluto.  The filename to be
 *                      loaded is 'name_prefix.emb'.
 * @param  startTime    Earliest time for which ephemeris needs to
 *                      be present in file.
 * @param  stopTime     Latest time for which ephemeris needs to be
 *                      present in file.
 * @param  ecfeciSys    Ecf/Eci utility service pointer that will be
 *                      copied into the Ephemeris object.
 *
 * @throws  runtime_error if 'name.emb' can't be opened or the format is
 *          invalid.
 */
std::unique_ptr<Ephemeris>
build_celestial(const std::string& name_prefix,
                const JulianDate& startTime,
                const JulianDate& stopTime,
                const std::shared_ptr<const EcfEciSys>& ecfeciSys);

/**
 * Parse NGS SP3-c compatible ephemeris.  'V' format ECF position and
 * velocity format is expected - position only will thrown an exception.
 * "EP" and "EV" fields are skipped.  Each "ID" must be the same
 * throughout the file or an exception will be trown.
 *
 * @param  file_name  Filename with SP3-c compatible ephemeris
 * @param  jdStart    Start time for which to store ephemeris records
 * @param  jdStop     End time for which to store ephemeris records
 *
 * @return  Position and velocity records, ECF, DU and DU/TU
 *
 * @throws  runtime_error if file_name can't be opened or the format is
 *          invalid.
 */
std::vector<state_vector_rec> parse_sp3_file(const std::string& file_name,
                                             const JulianDate& jdStart,
                                             const JulianDate& jdStop);

}

#endif
