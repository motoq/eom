#ifndef ASTRO_KEPLER_H
#define ASTRO_KEPLER_H

#include <array>

#include <Eigen/Dense>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>

namespace eom {

class Kepler : public Ephemeris {
public:
  Kepler(const JulianDate& epoch,
         const Eigen::Matrix<double, 6, 1>& stateVector);

  Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate& jd,
                                             EphemFrame frame) const override;

private:
  std::array<double, 4> planet = {0.0, 0.0, 0.0, 0.0};
  JulianDate jd0;
  std::array<double, 6> x0;
};


}

#endif
