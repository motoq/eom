/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_ECFECI_SYS_H
#define ASTRO_ECFECI_SYS_H

#include <vector>
#include <ostream>

#include <Eigen/Dense>
#include <Eigen/Geometry>

#include <cal_julian_date.h>
#include <cal_duration.h>

namespace eom {

/* Enable after adding parsing of IERS EOP
struct raw_eop {
  long mjd {0L};
  double xp {0.0};
  double yp {0.0};
  double ut1mutc {0.0};
  double lod {0.0};
  double dx {0.0};
  double dy {0.0};
};
*/

  // The transformation directoin is ECF to ECI
/**
 * ECF to ECI transformations
 */
struct ecf_eci {
  double mjd2000 {0.0};      ///< Modified Julian Date w.r.t. 2000, days
  double ut1mutc {0.0};      ///< UT1-UTC, seconds
  double lod {0.0};          ///< Length of day, seconds (for TIRF to CIRF)
  Eigen::Quaterniond pm;     ///< Polar motion; ITRF to TIRF
  Eigen::Quaterniond bpn;    ///< Frame bias, precession, nutation; CIRF to GCRF
};

/**
 * System resource utility for ECF/ECI conversions.
 *
 * @author  Kurt Motekew
 * @date    20211020
 */
class EcfEciSys {
public:
  /**
   * This constructor creates an ECF to ECI conversion utility that
   * makes use of a list of generated precession and nutation data.  It
   * allows for accurate conversions with some upfront computational burden.
   * Full benefit of this option comes with enabling interpolation.
   * Interpolation also comes at a price since it requires interpolating
   * ECF to ECI data vs. directly referencing an existing structure of
   * data in a list.
   *
   * @param  startTime    Earliest UTC time for which ECF to ECI data is
   *                      to be generated or available.
   * @param  startTime    Latest UTC time for which ECF to ECI data is
   *                      to be generated or available.
   * @param  dt           Rate at which to generate ECF to ECI data
   *                      If zero, one set of data is generated at the
   *                      center of the timeframe.
   * @param  interpolate  If true, ECF to ECI data will be interpolated.
   *                      Defaults to true.
   */
  EcfEciSys(const JulianDate& startTime, const JulianDate& stopTime,
            const Duration& dt, bool interpolate = true);

  /**
   * Returns an ecf_eci struucture.  Primarily intended for internal use
   * but public given the potential usefulness.
   *
   * @param  utc  Time for which to return and ecf_eci structure
   *
   * @return  Structure of ECF to ECI parameters.  All transformations
   *          are in the direction of ECF to ECI.
   */
  ecf_eci getEcfEciData(const JulianDate& utc) const;

  /**
   * Convert an ECF position vector to ECI.
   *
   * @param  utc   UTC time of position vector
   * @param  posf  Cartesian ECF position vector
   *
   * @return  Cartesian ECI position vector of same units as innput
   *          vector
   */
  Eigen::Matrix<double, 3, 1>
  ecf2eci(const JulianDate& utc, const Eigen::Matrix<double, 3, 1>& posf) const;

  /**
   * Convert an ECF position and velocity state vector to ECI.
   *
   * @param  utc   UTC time of state vector
   * @param  posf  Cartesian ECF position vector, DU
   * @param  posf  Cartesian ECF velocity vector, DU/TU
   *
   * @return  Cartesian ECI state vector, DU and DU/TU
   */
  Eigen::Matrix<double, 6, 1>
  ecf2eci(const JulianDate& utc, const Eigen::Matrix<double, 3, 1>& posf,
                                 const Eigen::Matrix<double, 3, 1>& velf) const;

  /**
   * Convert an ECI position vector to ECF.
   *
   * @param  utc   UTC time of position vector
   * @param  posf  Cartesian ECI position vector
   *
   * @return  Cartesian ECF position vector of same units as innput
   *          vector
   */
  Eigen::Matrix<double, 3, 1>
  eci2ecf(const JulianDate& utc, const Eigen::Matrix<double, 3, 1>& posi) const;

  /**
   * Convert an ECI position and velocity state vector to ECF.
   *
   * @param  utc   UTC time of state vector
   * @param  posf  Cartesian ECI position vector, DU
   * @param  posf  Cartesian ECI velocity vector, DU/TU
   *
   * @return  Cartesian ECF state vector, DU and DU/TU
   */
  Eigen::Matrix<double, 6, 1>
  eci2ecf(const JulianDate& utc, const Eigen::Matrix<double, 3, 1>& posi,
                                 const Eigen::Matrix<double, 3, 1>& veli) const;

  /**
   * Convert an ECF position and velocity state vector to true equator
   * true equinox (TEME) using the IAU 1982 gmst angular rotation.
   *
   * @param  utc   UTC time of state vector
   * @param  posf  Cartesian ECF position vector, DU
   * @param  posf  Cartesian ECF velocity vector, DU/TU
   *
   * @return  Cartesian TEME state vector, DU and DU/TU
   */
  Eigen::Matrix<double, 6, 1>
  ecf2teme(const JulianDate& utc,
           const Eigen::Matrix<double, 3, 1>& posf,
           const Eigen::Matrix<double, 3, 1>& velf) const;

  /**
   * Convert a TEME position and velocity state vector to ECF.
   *
   * @param  utc   UTC time of state vector
   * @param  posf  Cartesian TEME position vector, DU
   * @param  posf  Cartesian TEME velocity vector, DU/TU
   *
   * @return  Cartesian ECF state vector, DU and DU/TU
   */
  Eigen::Matrix<double, 6, 1>
  teme2ecf(const JulianDate& utc,
           const Eigen::Matrix<double, 3, 1>& posi,
           const Eigen::Matrix<double, 3, 1>& veli) const;

  /**
   * Print stored ECF2ECI data to the supplied stream.
   *
   * @out  Output stream
   */
  void print(std::ostream& out);

private:
  JulianDate jdStart;
  JulianDate jdStop;
  double rate_days {0.0};
  unsigned long nfi {0UL};
  bool interpolate_bpnpm {true};
  std::vector<ecf_eci> f2iData;
};


}

#endif