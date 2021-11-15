#ifndef ASTRO_EPHEMERIS_H
#define ASTRO_EPHEMERIS_H

#include <Eigen/Dense>

#include <cal_julian_date.h>

namespace eom {

enum class EphemFrame {
  eci,                            ///< GCRF (IAU 2000A/2006)
  ecf                             ///< ITRF (~WGS 84)
};

class Ephemeris {
public:
  virtual ~Ephemeris() {}

  virtual Eigen::Matrix<double, 6, 1> getStateVector(const JulianDate& jd,
                                                     EphemFrame frame) const=0;
};


}

#endif
