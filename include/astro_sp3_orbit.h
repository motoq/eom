/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_SP3_ORBIT_H
#define ASTRO_SP3_ORBIT_H

#include <string>
#include <vector>
#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <mth_hermite1.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <mth_index_mapper.h>

namespace eom {

/**
 * SP3 records used for generating interpolators
 */
struct sp3_rec {
  JulianDate t;
  Eigen::Matrix<double, 3, 1> p;            ///< Position
  Eigen::Matrix<double, 3, 1> v;            ///< Velocity

  sp3_rec(const JulianDate& jd,
          const Eigen::Matrix<double, 3, 1>& x,
          const Eigen::Matrix<double, 3, 1>& dx) : t(jd),
                                                   p(x),
                                                   v(dx)
  {
  }
};

/**
 * Interpolation records generated from ephemeris
 */
struct sp3_interp_rec {
  JulianDate jd1;                           ///< Interpolator start time
  JulianDate jd2;                           ///< Interpolator stop time
  Hermite1<double, 3> hItp;                 ///< Interpolator

  sp3_interp_rec(const JulianDate& jdStart,
                 const JulianDate& jdEnd,
                 const Hermite1<double, 3>& hInterp) : jd1(jdStart),
                                                       jd2(jdEnd),
                                                       hItp(hInterp)
  {
  }
};

/**
 * Parse NGS SP3-c compatible ephemeris.  'V' format ECF position and
 * velocity format is expected - position only will thrown an exception.
 * "EP" and "EV" fields are skipped.  Each "ID" must be the same
 * throughout the file or an exception will be trown.
 *
 * @author  Kurt Motekew  2023/01/10
 */
class Sp3Orbit : public Ephemeris {
public:
  ~Sp3Orbit() = default;
  Sp3Orbit(const Sp3Orbit&) = default;
  Sp3Orbit& operator=(const Sp3Orbit&) = default;
  Sp3Orbit(Sp3Orbit&&) = default;
  Sp3Orbit& operator=(Sp3Orbit&&) = default;

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
  Sp3Orbit(const std::string& name,
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
   * Interpolate state vector from stored ephemeris for given time
   *
   * @param  jd     Time of desired state vector, UTC
   * @param  frame  Desired output reference frame
   *
   * @return  Cartesian state vector at requested time in the requested
   *          reference frame, DU and DU/TU
   *
   * @throws  invalid_argument if the requested time is out of range
   */
  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate&,
                                             EphemFrame frame) const override;

private:
  std::string m_name {""};
  JulianDate m_jdStart;
  JulianDate m_jdStop;
  JulianDate m_jdEpoch;
  std::shared_ptr<const EcfEciSys> m_ecfeciSys {nullptr};

  std::unique_ptr<IndexMapper<JulianDate>> m_ndxr {nullptr};
  std::vector<sp3_interp_rec> m_eph_interpolators;
};


}

#endif