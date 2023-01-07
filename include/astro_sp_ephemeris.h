/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_SP_EPHEMERIS_H
#define ASTRO_SP_EPHEMERIS_H

#include <string>
#include <vector>
#include <memory>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <mth_hermite2.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <mth_ode_solver.h>
#include <mth_index_mapper.h>

namespace eom {

/**
 * Ephemeris records used for generating interpolators
 */
struct eph_record {
  JulianDate t;
  Eigen::Matrix<double, 3, 1> p;            ///< Position
  Eigen::Matrix<double, 3, 1> v;            ///< Velocity
  Eigen::Matrix<double, 3, 1> a;            ///< Acceleration

  eph_record(const JulianDate& jd,
             const Eigen::Matrix<double, 3, 1>& x,
             const Eigen::Matrix<double, 3, 1>& dx,
             const Eigen::Matrix<double, 3, 1>& ddx) : t(jd),
                                                       p(x),
                                                       v(dx),
                                                       a(ddx)
  {
  }
};

/**
 * Interpolation records generated from ephemeris
 */
struct interp_record {
  JulianDate jd1;                           ///< Interpolator start time
  JulianDate jd2;                           ///< Interpolator stop time
  Hermite2<double, 3> hItp;                 ///< Interpolator

  interp_record(const JulianDate& jdStart,
                const JulianDate& jdEnd,
                const Hermite2<double, 3>& hInterp) : jd1(jdStart),
                                                      jd2(jdEnd),
                                                      hItp(hInterp)
  {
  }
};

/**
 * Generates ephemeris through special perturbations methods and stores
 * as interpolators for retrieval.  Position, velocity, and acceleration
 * are used to form Hermite interpolators.
 *
 * @author  Kurt Motekew  2022/12/26
 */
class SpEphemeris : public Ephemeris {
public:
  ~SpEphemeris() = default;
  SpEphemeris(const SpEphemeris&) = default;             // copy constructor
  SpEphemeris& operator=(const SpEphemeris&) = default;  // copy assignment
  SpEphemeris(SpEphemeris&&) = default;                  // move constructor
  SpEphemeris& operator=(SpEphemeris&&) = default;       // move assignment

  /**
   * Initialize with orbital state and model/integrator.  Generate
   * ephemeris from jdStart to jdStop.
   *
   * @param  name       Unique ephemeris identifier
   * @param  jdStart    Start time for which ephemeris should be created
   * @param  jdStop     End time for which ephemeris should be created
   * @param  ecfeciSys  ECF/ECI conversion resource
   * @param  sp         Integrator with force model (EOM) used to
   *                    generate ephemeris.  SpEphemeris takes
   *                    ownership.
   */
  SpEphemeris(const std::string& name,
              const JulianDate& jdStart,
              const JulianDate& jdStop,
              const std::shared_ptr<const EcfEciSys>& ecfeciSys,
              std::unique_ptr<OdeSolver<JulianDate, double, 6>> sp);

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
  JulianDate m_jdEpoch;
  JulianDate m_jdStart;
  JulianDate m_jdStop;
  std::shared_ptr<const EcfEciSys> m_ecfeciSys {nullptr};

  std::unique_ptr<IndexMapper<JulianDate>> m_ndxr {nullptr};
  std::vector<interp_record> m_eph_interpolators;
};


}

#endif
