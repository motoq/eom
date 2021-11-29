#ifndef ASTRO_VINTI_H
#define ASTRO_VINTI_H

#include <array>
#include <memory>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_ecfeci_sys.h>


namespace eom {

/*
 * Supplemental Material
 * https://arc.aiaa.org/doi/suppl/10.2514/4.866487
 *
 * Orbital and Celestial Mechanics
 * Nino L. Bonavito, Gim J. Der, and John P. Vinti
 * AIAA, 1998
 */
class Vinti : public Ephemeris {
public:
  Vinti(const JulianDate& epoch,
        const Eigen::Matrix<double, 6, 1>& xeci,
        const std::shared_ptr<const EcfEciSys>& ecfeciSys);

  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate& jd,
                                             EphemFrame frame) const override;

private:
  std::shared_ptr<const EcfEciSys> ecfeci {nullptr};
  std::array<double, 4> planet = {phy_const::km_per_du,
                                  phy_const::gm_km3_sec2,
                                  phy_const::j2, phy_const::j3};
  JulianDate jd0;
  std::array<double, 6> x0;
};


}

#endif
