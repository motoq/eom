#include <astro_kepler.h>

#include <array>
#include <cmath>

#include <Vinti.h>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>

namespace eom {


Kepler::Kepler(const JulianDate& epoch,
               const Eigen::Matrix<double, 6, 1>& stateVector)
{
  planet[0] = 1.0;
  planet[1] = 1.0;
  jd0 = epoch;
  x0[0] = stateVector(0,0);
  x0[1] = stateVector(1,0);
  x0[2] = stateVector(2,0);
  x0[3] = stateVector(3,0);
  x0[4] = stateVector(4,0);
  x0[5] = stateVector(5,0);
}


Eigen::Matrix<double, 6, 1> Kepler::getStateVector(const JulianDate& jd,
                                                   EphemFrame frame) const 
{
  std::array<double, 6> x1;
  double x;
  double t1 {phy_const::tu_per_day*(jd - jd0)};
  Kepler1(planet.data(), 0.0, x0.data(), t1, x1.data(), &x);

  Eigen::Matrix<double, 6, 1> stateVector;
  stateVector(0,0) = x1[0];
  stateVector(1,0) = x1[1];
  stateVector(2,0) = x1[2];
  stateVector(3,0) = x1[3];
  stateVector(4,0) = x1[4];
  stateVector(5,0) = x1[5];

  if (frame == EphemFrame::ecf) {
    stateVector *= 1.0;
  }

  return stateVector;
}



}
