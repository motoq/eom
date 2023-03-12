/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_SP3_HERMITE_H
#define ASTRO_SP3_HERMITE_H

#include <string>
#include <vector>
#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <mth_hermite2.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <mth_index_mapper.h>

namespace eom {

/**
 * Interpolation records generated from ephemeris
 */
struct sp3_hermite {
  JulianDate jd1;                           ///< Interpolator start time
  JulianDate jd2;                           ///< Interpolator stop time
  Hermite2<double, 3> hItp;                 ///< Interpolator

  sp3_hermite(const JulianDate& jdStart,
              const JulianDate& jdEnd,
              const Hermite2<double, 3>& hInterp) : jd1(jdStart),
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
 * Hermite interpolation using two position and velocity pairs is
 * employed, augmented with a J4 gravity model to create acceleration.
 * This method was chosen because the spacing for a 2 rev/day orbit was
 * 15 minutes, causing significant error with position-velocity only
 * Hermite interpolation.  Hermite interpolation does limit this class
 * to ephemeris where velocity is the derivative of position, and the J4
 * acceleration model is a sufficient representation of acceleration.
 * For a 2 rev/day orbit with 15 minute state vector spacing, results in
 * an oscillation of about 5 cm.
 *
 * @author  Kurt Motekew  2023/01/10
 */
class Sp3Hermite : public Ephemeris {
public:
  ~Sp3Hermite() = default;
  Sp3Hermite(const Sp3Hermite&) = default;
  Sp3Hermite& operator=(const Sp3Hermite&) = default;
  Sp3Hermite(Sp3Hermite&&) = default;
  Sp3Hermite& operator=(Sp3Hermite&&) = default;

  /**
   * Initialize with SP3 compatible format ephemeris.
   *
   * @param  name         Unique ephemeris identifier
   * @param  sp3_records  Position and velocity records to form Hermite
   *                      interpolation polynomials.  At least two must
   *                      be present and must cover jdStart and jdStop.
   *                      ECF, DU and DU/TU
   * @param  jdStart      Start time for which ephemeris must be available
   * @param  jdStop       End time for which ephemeris must be available
   * @param  ecfeciSys    ECF/ECI conversion resource
   *
   * @throws  runtime_error for parsing and processing errors
   */
  Sp3Hermite(const std::string& name,
             const std::vector<state_vector_rec>& sp3_records,
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
  std::vector<sp3_hermite> m_eph_interpolators;
};


}

#endif
