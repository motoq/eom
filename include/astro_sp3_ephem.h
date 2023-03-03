/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_SP3_EPHEM_H
#define ASTRO_SP3_EPHEM_H

#include <string>
#include <vector>
#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_granule.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <mth_index_mapper.h>

namespace eom {

namespace sp3 {
  constexpr int order {5};
  constexpr int np {9};
}

/**
 * Interpolation granules
 */
struct sp3_granule {
  JulianDate jd1;                           ///< Interpolator start time
  JulianDate jd2;                           ///< Interpolator stop time
  Granule<sp3::order, sp3::np> tItp;        ///< Interpolator

  sp3_granule(const JulianDate& jdStart,
              const JulianDate& jdEnd,
              const Granule<sp3::order, sp3::np>& interp) : jd1(jdStart),
                                                            jd2(jdEnd),
                                                            tItp(interp)
  {
  }
};

/**
 * Parse NGS SP3-c compatible ephemeris.  'V' format ECF position and
 * velocity format is expected - position only will thrown an exception.
 * "EP" and "EV" fields are skipped.  Each "ID" must be the same
 * throughout the file or an exception will be trown.
 *
 * Chebyshev interpolation using
 *
 * @author  Kurt Motekew  2023/02/28
 */
class Sp3Ephem : public Ephemeris {
public:
  ~Sp3Ephem() = default;
  Sp3Ephem(const Sp3Ephem&) = default;
  Sp3Ephem& operator=(const Sp3Ephem&) = default;
  Sp3Ephem(Sp3Ephem&&) = default;
  Sp3Ephem& operator=(Sp3Ephem&&) = default;

  /**
   * Initialize with SP3 compatible format ephemeris.
   *
   * @param  name       Unique ephemeris identifier
   * @param  file_name  Filename with SP3-c compatible ephemeris.
   * @param  jdStart    Start time for which ephemeris must be available
   * @param  jdStop     End time for which ephemeris must be available
   * @param  ecfeciSys  ECF/ECI conversion resource
   *
   * @throws  runtime_error for parsing and processing errors
   */
  Sp3Ephem(const std::string& name,
           const std::string& file_name,
           const JulianDate& jdStart,
           const JulianDate& jdStop,
           const std::shared_ptr<const EcfEciSys>& ecfeciSys);

  /**
   * @return  Unique ephemeris identifier
   */
  std::string getName() const override
  {
    return m_name;
  }

  /**
   * @return  Get initializing state vector time
   */
  JulianDate getEpoch() const override
  {
    return m_jdEpoch;
  }

  /**
   * @return  Earliest time for which ephemeris can be retrieved
   */
  JulianDate getBeginTime() const override
  {
    return m_jdStart;
  }

  /**
   * @return  Latest time for which ephemeris can be retrieved
   */
  JulianDate getEndTime() const override
  {
    return m_jdStop;
  }

  /**
   * Interpolate state vector from stored ephemeris for given time
   *
   * @param  jd     Time of desired state vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian state vector at requested time in the requested
   *          reference frame, DU and DU/TU
   *
   * @throws  out_of_range if the requested time is out of range
   */
  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate&,
                                             EphemFrame frame) const override;

  /**
   * @param  jd     Time of desired position vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian position vector, DU
   *
   * @throws  out_of_range if the requested time is out of range
   */
  Eigen::Matrix<double, 3, 1> getPosition(const JulianDate& jd,
                                          EphemFrame frame) const override;

private:
  std::string m_name {""};
  JulianDate m_jdStart;
  JulianDate m_jdStop;
  JulianDate m_jdEpoch;
  std::shared_ptr<const EcfEciSys> m_ecfeciSys {nullptr};

  std::unique_ptr<IndexMapper<JulianDate>> m_ndxr {nullptr};
  std::vector<sp3_granule> m_eph_interpolators;
};


}

#endif
