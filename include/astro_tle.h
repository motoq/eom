/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_TLE_H
#define ASTRO_TLE_H

#include <array>
#include <ostream>
#include <string>

#include <cal_julian_date.h>

namespace eom {

/**
 * Two Line Element set class.  Note accessor functions to individual
 * TLE elements return in published units vs. eom internal computational
 * units.
 *
 * @author  Kurt Motekew @date    2024/07/23
 */
class Tle {
public:

  /**
   * Default orbit
   */
  Tle() { }

  /**
   * Initialize with a NORAD two line element set.
   *
   * @param  tle1  Line 1 of element set.  Must be at least 63 characters
   * @param  tle2  Line 2 of element set.  Must be at least 63 characters
   *
   * @throw invalid_argument If an error is encounterd during parsing
   */
  Tle(const std::string& tle1, const std::string& tle2);

  /**
   * @return  Minimum length of each TLE line required by parser
   */
  static unsigned int minimumLineLength() { return 63; }

  /**
   * Same as constructor, but with an existing TLE object.
   */
  void set(const std::string& tle1, const std::string& tle2);

  /**
   * The SGP4 epoch is used to compute the independent variable used for
   * orbit propagation given a Julian date.
   *
   * @return  SGP4 epoch, Jan 1 1950 00:00:00, UTC
   */
  JulianDate getSgpEpoch() const;

  /**
   * @return  TLE line 1, as was used for initialization.
   */
  std::string getLineOne() const;

  /**
   * @return  TLE line 2, as was used for initialization.
   */
  std::string getLineTwo() const;

  /**
   * @return  Alphanumeric satellite designator
   */
  std::string getName() const;

  /**
   * @return  Alphanumeric satellite designator in form more compatible
   *          with legacy SGP related codebases.
   */
  std::array<char, 6> getSatN() const;

  /**
   * @return  TLE epoch
   */
  JulianDate getEpoch() const;

  /**
   * @return  Two digit TLE epoch year, as it appears in the TLE
   */
  int getEpochYear() const;

  /**
   * @return  TLE epoch day of year and fraction
   */
  double getEpochDayOfYear() const;

  /**
   * @return  TLE epoch, days from SGP4 epoch
   */
  double getTleSgpEpoch() const;

  /**
   * @return  First time derivative of mean motion fudge factor,
   *          1/2 rev/day^2
   */
  double getMeanMotionRate() const;

  /**
   * @return  Second time derivative of mean motion fudge factor,
   *          1/6 rev/day^3
   */
  double getMeanMotionSecondRate() const;

  /**
   * @return  Ballistic coefficient fudge factor, 1/ER
   */
  double getBstar() const;

  /**
   * @return  "Ephemeris" type.  0 indicates the Kozai mean elements
   *          expected by most SGP/SGP4 initialization stages.  Brouwer
   *          mean elements are indicated by a 2.  4 has been chosen to
   *          designate the new USSF SGP4-XP propagator.
   */
  int getElementType() const;

  /**
   * @return  Inclination, degrees
   */
  double getInclination() const;

  /**
   * @return  Right ascension of the ascending node, degrees
   */
  double getRaan() const;

  /**
   * @return  Eccentricity
   */
  double getEccentricity() const;

  /**
   * @return  Argument of perigee, deg
   */
  double getArgumentOfPerigee() const;

  /**
   * @return  Mean anomaly, deg
   */
  double getMeanAnomaly() const;

  /**
   * @return  Mean motion, rev/day
   */
  double getMeanMotion() const;

private:
  std::string m_tle1 {""};   // TLE line 1
  std::string m_tle2 {""};   // TLE line 2
  std::string m_satn {""};   // Satellite alphanumeric designator
  int m_epoch_yr {0};        // Two digit year
  double m_epoch_doy {1.0};  // Days
  double m_ndot {0.0};       // d(no)/dt
  double m_nddot {0.0};      // d2(no)/(dtdt)
  double m_bstar {0.0};      // B*
  int m_eph_type {0};        //
  double m_inclo {63.4};     // Incliniation, deg
  double m_nodeo {0.0};      // RAAN, deg
  double m_ecco {0.0001};    // Eccentricity
  double m_argpo {0.0};      // Argument of perigee, deg
  double m_m0 {0.0};         // Mean anomaly, deg
  double m_n0 {2.0};         // Mean motion, rev/day
    //
  JulianDate m_epoch;
  
};

/**
 * Override output stream as non-member function.
 */
std::ostream& operator<<(std::ostream& out, const Tle& kep);

}

#endif
