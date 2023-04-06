/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_VINTI_PROP_H
#define ASTRO_VINTI_PROP_H

#include <string>
#include <array>
#include <memory>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>
#include <astro_kepler_prop.h>


namespace eom {

/*
 * A C++ version of the the Gim J. Der & Herbert B. Reynolds Vinti6
 * analytic propagator.  See the Vinti6.c source code for the original.
 * This version differs from astro_vinti.h/cpp in that the Vinti6 code
 * from has been separated into an initialization portion executed once
 * in the constructor along with the Vinti6.c code being brought into
 * the class definition.  J3 effects were also removed.  The goal with
 * astro_vinti.cpp was to be a wrapper about the original implementation
 * while no effort is made with this version to honor the original
 * structure of the code.
 *
 * Vinti6.c and Vinti.h source:
 * Supplemental Material
 * https://arc.aiaa.org/doi/suppl/10.2514/4.866487
 *
 * Orbital and Celestial Mechanics
 * Nino L. Bonavito, Gim J. Der, and John P. Vinti
 * AIAA, 1998
 */
class VintiProp : public Ephemeris {
public:
  ~VintiProp() = default;
  VintiProp(const VintiProp&) = default;
  VintiProp& operator=(const VintiProp&) = default;
  VintiProp(VintiProp&&) = default;
  VintiProp& operator=(VintiProp&&) = default;

  /**
   * Initialize VintiProp
   *
   * @param  orbit_name      Name (string identifieer) associated with orbit
   * @param  epoch           Orbit state vector epoch, UTC
   * @param  xeci            Cartesian GCRF state vector, DU and DU/TU
   * @param  ecfeciSys       ECF/ECI conversion resource
   *
   * @throws  std::invalid_argument if orbit definition falls within the
   *          Vinti forbidden zone.
   */
  VintiProp(const std::string& orbit_name,
            const JulianDate& epoch,
            const Eigen::Matrix<double, 6, 1>& xeci,
            const std::shared_ptr<const EcfEciSys>& ecfeciSys);

  /**
   * @return  Unique ephemeris identifier
   */
  std::string getName() const override
  {
    return name;
  }

  /**
   * @return  Orbit epoch
   */
  JulianDate getEpoch() const override
  {
    return jd0;
  }

  /**
   * @return  Earliest time for which ephemeris can be retrieved
   */
  JulianDate getBeginTime() const override
  {
    return ecfeci->getBeginTime();
  }

  /**
   * @return  Latest time for which ephemeris can be retrieved
   */
  JulianDate getEndTime() const override
  {
    return ecfeci->getEndTime();
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

  /**
   * @return  Vinti mean elements at orbit epoch
   */
  std::array<double, 6> getVintiMean() const
  {
    return oe0;
  }
  
private:
  void vinti_local(const JulianDate& jd, std::array<double, 6>& x1) const;

  std::string name {""};
  std::shared_ptr<const EcfEciSys> ecfeci {nullptr};
  std::array<double, 4> planet = {phy_const::re,
                                  phy_const::gm,
                                  phy_const::j2, 0.0};
  JulianDate jd0;
  std::array<double, 3> pin;
  std::array<double, 3> vin;

  std::unique_ptr<KeplerProp> m_kep {nullptr};

    // t0 at jd0
  double t0 {0.0};
    // Vinti mean elements at epoch
  std::array<double, 6> oe0;
    // Initialization variables - constant after constructor
  double xhat0 {0.0};
    // Steps 1-3
  double alph3 {};
  double alph2 {};
  double b1 {};
  double a1 {};
  double gamma {};
  double gam1 {};
  double p {};
  double smgam {};
  double s1 {};
  double q1 {};
  double p1 {};
    // Step 4
  double gams3 {};
  //double s {};
  double q {0.0};
  double q2 {};
  double g {};
  double a {};
  double b {};
  double g2 {};
  double ecc2 {};
  double ecc {};                            // Vinti mean oe
  double rho1 {};
    // Wi for the R - integrals
  double x21 {};
  double x22 {};
  double x33 {};
  double x32 {};
  double x31 {};
  double x44 {};
  double x43 {};
  double x42 {};
  double x41 {};
  double x55 {};
  double x54 {};
  double x53 {};
  double x52 {};
  double x51 {};
  double x66 {};
  double x65 {};
  double x64 {};
  double x63 {};
  double x62 {};
  double x61 {};
  double x77 {};
  double x76 {};
  double x75 {};
  double x74 {};
  double x73 {};
  double x72 {};
  double x71 {};
    // R1 Coefficients
  double cr11 {};
  double cr12 {};
  double cr13 {};
  double cr14 {};
  double cr15 {};
  double cr16 {};
  double cr17 {};
    // R2 Coefficients
  double cr21 {};
  double cr22 {};
  double cr23 {};
  double cr24 {};
  double cr25 {};
  double cr26 {};
  double cr27 {};
    // R3 Coefficients
  double cr31 {};
  double cr32 {};
  double cr33 {};
  double cr34 {};
  double cr35 {};
    // N3 Coefficients
  double xmm1 {};
  double xmm2 {};
  double d10 {};
  double d20 {};
  double b1q {};
  double b2q {};
    //
  double ucf1 {};
  double ucf2 {};
  double ucf3 {};
  double denyst {};
    // N1 Coefficients
  double cn11 {};
  double cn12 {};
  double cn13 {};
  double cn14 {};
  double cn15 {};
  double cn16 {};
  double cn17 {};
    // N2 Coefficients
  double cn31 {};
  double cn32 {};
  double cn33 {};
  double cn34 {};
  double cn35 {};
  double cn36 {};
    // Step 5
  double somega {};
  double capt {};
  double comega {};
};


}

#endif
