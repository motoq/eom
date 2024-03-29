/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_VINTI_H
#define ASTRO_VINTI_H

#include <string>
#include <array>
#include <memory>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>


namespace eom {

enum class VintiPertModel {
  J2_ONLY,                             ///< Include only J2 effects
  J2_J3                                ///< Vinti J2 with J3 approx
};

/*
 * Wrapper class to the Gim J. Der, Albert T. Monuki, & Herb Reynolds
 * Vinti (Vinti6) orbit propagator.  See the Vinti6.c source code.
 *
 * Vinti6.c and Vinti.h source:
 * Supplemental Material
 * https://arc.aiaa.org/doi/suppl/10.2514/4.866487
 *
 * Orbital and Celestial Mechanics
 * Nino L. Bonavito, Gim J. Der, and John P. Vinti
 * AIAA, 1998
 */
class Vinti : public Ephemeris {
public:
  ~Vinti() = default;
  Vinti(const Vinti&) = default;
  Vinti& operator=(const Vinti&) = default;
  Vinti(Vinti&&) = default;
  Vinti& operator=(Vinti&&) = default;

  /**
   * Initialize Vinti6
   *
   * @param  orbit_name      Name (string identifieer) associated with orbit
   * @param  epoch           Orbit state vector epoch, UTC
   * @param  xeci            Cartesian GCRF state vector, DU and DU/TU
   * @param  ecfeciSys       ECF/ECI conversion resource
   * @param  VintiPertModel  Inidcates if J2 only or J2 with the J3
   *                         approximation should be included in the
   *                         gravity model.
   */
  Vinti(const std::string& orbit_name,
        const JulianDate& epoch,
        const Eigen::Matrix<double, 6, 1>& xeci,
        std::shared_ptr<const EcfEciSys> ecfeciSys,
        VintiPertModel pertModel = VintiPertModel::J2_J3);

  /**
   * @return  Unique ephemeris identifier
   */
  std::string getName() const override
  {
    return m_name;
  }

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
  std::shared_ptr<const EcfEciSys> m_ecfeci {nullptr};
  std::array<double, 4> m_planet = {phy_const::km_per_du,
                                    phy_const::gm_km3_sec2,
                                    phy_const::j2, phy_const::j3};
  JulianDate m_jd0;
  std::array<double, 6> m_x0;
};


}

#endif
