/*
 * Copyright 2021, 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_PROPAGATOR_CONFIG_H
#define ASTRO_PROPAGATOR_CONFIG_H

#include <cal_duration.h>
#include <cal_julian_date.h>

namespace eom {

/**
 * Propagator options
 */
enum class PropagatorType {
#ifdef GENPL
  sec_j2,
  osc_j2,
#endif
  fandg,                          ///< BMW f and g method
  secular_j2,                     ///< J2 Secular effects
  gpx,                            ///< Testing ground for GP propagation
  kepler1,                        ///< Gim J. Der version
  vinti6,                         ///< Gim J. Der version
  vinti_j2,                       ///< Gim J. Der version with J3 = 0
  kepler1mod,                     ///< Experimental, modified Gim J. Der version
  vinti6mod,                      ///< Experimental, modified Gim J. Der version
  sgp4,                           ///< Based on Vallado's version
  sp                              ///< Special Pert with more config options
};

/**
 * Numerical integration options 
 */
enum class Propagator {
#ifdef GENPL
  gj,                             ///< Gauss-Jackson
  gjs,                            ///< G-J with time regularization
#endif
  rk4,                            ///< Basic RK4 integration
  rk4s,                           ///< RK4 with time regularization
  adams4                          ///< Adams-Bashforth-Moulton
};

/**
 * Gravity model options
 */
enum class GravityModel {
#ifdef GENPL
  gravt,
#endif
  jn,                             ///< Simple zonal-only gravity model
  std,                            ///< Degree, Order gravity model
};

/**
 * Sun gravity model options
 */
enum class SunGravityModel {
  none,
  meeus,                          ///< Analytic Astronomical Algorithms
  eph                             ///< sun.emb file
};

/**
 * Moon gravity model options
 */
enum class MoonGravityModel {
  none,
  meeus,                          ///< Analytic Astronomical Algorithms
  eph                             ///< moon.emb file
};

/**
 * Solar radiation pressure model options
 */
enum class SrpModel {
  none,
  spherical                       ///< Two parameter spherical
};


/**
 * Contains propagator configuration parameters
 */
class PropagatorConfig {
public:
  /**
   * Initilize with default propagator
   */
  PropagatorConfig();

  /**
   * Initialize with only the propagator type.  Most general
   * perturbation propagators should use this version.  Special
   * perturbation propagators will be initialized with default values.
   */
  PropagatorConfig(PropagatorType prop_type);

  /*
   * @return  The propagator type
   */
  PropagatorType getPropagatorType() const noexcept;

  /**
   * Set time span over which ephemeris should be valid - typically
   * only applies to SP based ephemeris, indicating integration limits.
   *
   * @param  jdStart  Start time, UTC
   * @param  jdStop   Stop time, UTC
   */
  void setStartStopTime(const JulianDate& jdStart,
                        const JulianDate& jdStop);

  /**
   * @return  Required ephemeris start time when generated via SP
   *          methods, UTC
   */
  JulianDate getStartTime() const noexcept;

  /**
   * @return  Required ephemeris stop time when generated via SP
   *          methods, UTC
   */
  JulianDate getStopTime() const noexcept;

  /**
   * @param  Set the integration method to use
   */
  void setPropagator(Propagator integration_method);

  /**
   * @return  Integration method to use
   */
  Propagator getPropagator() const noexcept;

  /**
   * @return  Suggessted (starting) integration step size.  Defaults
   *          to zero if not explicitly set.
   */
  void setStepSize(const Duration& dt);

  /**
   * @return  Suggested integration step size.  Default value of
   *          zero typically interpreted as indicator to use
   *          default or self determined integration step size.
   */
  Duration getStepSize() const noexcept;

  /**
   * @param  Set the gravity model to use
   */
  void setGravityModel(GravityModel gravity_model);

  /**
   * @return  Gravity model to use
   */
  GravityModel getGravityModel() const noexcept;

  /**
   * @param  Set the solar gravity model to use
   */
  void setSunGravityModel(SunGravityModel sun_gravity);

  /**
   * @return  Solar gravity model to use
   */
  SunGravityModel getSunGravityModel() const noexcept;

  /**
   * @param  Set the lunar gravity model to use
   */
  void setMoonGravityModel(MoonGravityModel moon_gravity);

  /**
   * @return  Lunar gravity model to use
   */
  MoonGravityModel getMoonGravityModel() const noexcept;

  /**
   * When called, enables other gravity models based on celestial bodies
   * initialized via external ephemerides
   */
  void enableOtherGravityModels() noexcept;

  /**
   * When called, disables other gravity models based on celestial bodies
   * initialized via external ephemerides
   */
  void disableOtherGravityModels() noexcept;

  /**
   * @return  true if planetary, etc., gravity models are enabled
   */
  bool otherGravityModelsEnabled() const noexcept;

  /**
   * Order <= Degree
   *
   * @param  Degree of gravitational coefficients to consider
   * @param  Order of gravitational coefficients to consider
   */
  void setDegreeOrder(int degree, int order);

  /**
   * @return  Degree of gravity model
   */
  int getDegree() const noexcept;

  /**
   * @return  Order of gravity model
   */
  int getOrder() const noexcept;

  /**
   * @param  Set the SRP model to use
   */
  void setSrpModel(SrpModel srp_model);

  /**
   * @return  Selected SRP model
   */
  SrpModel getSrpModel();

  /**
   * @param  Set reflectivity for fixed value SRP models
   */
  void setReflectivity(double cr);

  /**
   * @return  Reflectivity for fixed value SRP models
   */
  double getReflectivity();

  /**
   * @param  Set fixed plate area over mass for SRP models.
   */
  void setAreaOverMass(double aom);

  /**
   * @param  Fixed plate area over mass for SRP models.
   */
  double getAreaOverMass();

private:
    // Required for all propagators
  PropagatorType m_prop_type {PropagatorType::kepler1};
    // Typically required only for SP methods to set integration limits
  JulianDate m_jdStart;
  JulianDate m_jdStop;
    // Integration method and step size
  Propagator m_propagator {Propagator::rk4};
  Duration m_dt;
    // Gravity model
  GravityModel m_gravity_model {GravityModel::jn};
  SunGravityModel m_sun_gravity {SunGravityModel::none};
  MoonGravityModel m_moon_gravity {MoonGravityModel::none};
  bool m_other_gravity {false};
    // Non-conservative forces
  SrpModel m_srp_model {SrpModel::none};
  double m_cr {0.0};                        // Reflectivity
  double m_aom {0.0};                       // Area over Mass

  int m_degree {0};
  int m_order {0};
};


}

#endif
