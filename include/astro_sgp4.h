/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_SGP4_H
#define ASTRO_SGP4_H

#include <array>
#include <memory>
#include <string>

#include <Eigen/Dense>

#include <SGP4.h>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_tle.h>

namespace eom {

enum class Sgp4EarthModel {
  wgs72,                          ///< Legacy compatibility mode?
  wgs84                           ///< Updated constants
};

enum class Sgp4Mode {
  afspc,
  improved
};

/*
 * Wrapper class to the David Vallado's implementation hosted by
 * CelesTrak (https://celestrak.org/).
 */
class Sgp4 : public Ephemeris {
public:
  ~Sgp4() = default;
  Sgp4(const Sgp4&) = default;
  Sgp4& operator=(const Sgp4&) = default;
  Sgp4(Sgp4&&) = default;
  Sgp4& operator=(Sgp4&&) = default;

  /**
   * Initialize SGP4 propagator.  Kozai mean elements sets (Eph Type = 0)
   * expected.
   *
   * @param  orbit_name   Name (string identifier) associated with orbit
   *                      Not necessarily a NORAD designator
   * @param  tle1         NORAD formatted TLE line 1
   * @param  tle2         NORAD formatted TLE line 2
   * @param  ecfeciSys    ECF/ECI conversion resource
   * @param  earth_model  Inidcates source of earth model constants
   * @param  mode         Legacy AFSPC or more recent mode
   *
   * @throw invalid_argument  If an error is encountered during
   *                          initialization.
   */
  Sgp4(const std::string& orbit_name,
       const Tle& tle,
       std::shared_ptr<const EcfEciSys> ecfeciSys,
       Sgp4EarthModel earth_model = Sgp4EarthModel::wgs72,
       Sgp4Mode mode = Sgp4Mode::afspc);

  /**
   * @return  Unique ephemeris identifier
   */
  std::string getName() const override;

  /**
   * @return  Orbit epoch
   */
  JulianDate getEpoch() const override
  {
    return m_jd0;
  }

  /**
   * @return  Earliest time for which ephemeris can be retrieved
   */
  JulianDate getBeginTime() const override
  {
    return m_ecfeci->getBeginTime();
  }

  /**
   * @return  Latest time for which ephemeris can be retrieved
   */
  JulianDate getEndTime() const override
  {
    return m_ecfeci->getEndTime();
  }

  /**
   * Compute state vector given a time
   *
   * @param  jd     Time of desired state vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian state vector at requested time in the requested
   *          reference frame, DU and DU/TU
   *
   * @throws  out_of_range if the requested time is out of range.  This
   *          would be due to a time for which ECF/ECI transformation
   *          data is not available vs. a limitation of the analytic
   *          propagator itself.
   */
  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate& jd,
                                             EphemFrame frame) const override;

  /**
   * @param  jd     Time of desired position vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian position vector, DU
   *
   * @throws  out_of_range if the requested time is out of range.  This
   *          would be due to a time for which ECF/ECI transformation
   *          data is not available vs. a limitation of the analytic
   *          propagator itself.
   */
  Eigen::Matrix<double, 3, 1> getPosition(const JulianDate& jd,
                                          EphemFrame frame) const override;
private:
  std::string m_name;
  Tle m_tle;
  JulianDate m_jd0;
  std::shared_ptr<const EcfEciSys> m_ecfeci {nullptr};
  elsetrec m_satrec;
};


}

#endif
