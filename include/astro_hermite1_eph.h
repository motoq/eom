/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_HERMITE1_EPH_H
#define ASTRO_HERMITE1_EPH_H

#include <memory>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include <mth_hermite1.h>
#include <mth_index_mapper.h>
#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>

namespace eom {

/**
 * Interpolation records generated from ephemeris
 */
struct hermite1_eph_rec {
  JulianDate jd1;                           ///< Interpolator start time
  JulianDate jd2;                           ///< Interpolator stop time
  Hermite1<double, 3> hItp;                 ///< Interpolator

  hermite1_eph_rec(const JulianDate& jdStart,
                   const JulianDate& jdEnd,
                   const Hermite1<double, 3>& hInterp) : jd1(jdStart),
                                                         jd2(jdEnd),
                                                         hItp(hInterp)
  {
  }
};

/**
 * Hermite interpolation using two position and velocity pairs
 * is employed to generate a state vector from a set of ECI
 * ephemeris given a time.  No augmentation via the use of an
 * external gravity model is employed making this class useful
 * for non-earth centered objects provided the state vectors are
 * appropriately spaced.
 *
 * @author  Kurt Motekew  2023/12/30
 */
class Hermite1Eph : public Ephemeris {
public:
  /**
   * Initialize with ECI position/velocity based ephemeris records.
   *
   * @param  name         Unique ephemeris identifier
   * @param  eph_records  Position and velocity records to form Hermite
   *                      interpolation polynomials.  At least two must
   *                      be present and must cover jdStart and jdStop.
   *                      ECI, DU and DU/TU
   * @param  jdStart      Start time for which ephemeris must be available
   * @param  jdStop       End time for which ephemeris must be available
   * @param  ecfeciSys    ECF/ECI conversion resource
   *
   * @throws  runtime_error for initialization error
   */
  Hermite1Eph(const std::string& name,
              const std::vector<state_vector_rec>& pv_records,
              const JulianDate& jdStart,
              const JulianDate& jdStop,
              std::shared_ptr<const EcfEciSys> ecfeciSys);

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
  std::string m_name;
  JulianDate m_jdStart;
  JulianDate m_jdStop;
  JulianDate m_jdEpoch;
  std::shared_ptr<const EcfEciSys> m_ecfeciSys {nullptr};

  std::unique_ptr<IndexMapper<JulianDate>> m_ndxr {nullptr};
  std::vector<hermite1_eph_rec> m_eph_interpolators;
};


}

#endif
