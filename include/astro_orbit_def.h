#ifndef ASTRO_ORBIT_DEF_H
#define ASTRO_ORBIT_DEF_H

#include <array>
#include <string>

#include <cal_julian_date.h>
#include <astro_propagator_config.h>

namespace eom {


enum class CoordType {
  cartesian
};

enum class FrameType {
  gcrf                            ///< GCRF (IAU 2000A/2006)
};



class OrbitDef {
public:
  OrbitDef(std::string orbit_name, const PropagatorConfig& propConfig,
           const JulianDate& epoch, const std::array<double, 6>& state,
                            CoordType coord_type, FrameType frame_type);

  std::string getOrbitName() const noexcept { return name; }

  PropagatorConfig getPropagatorConfig() const noexcept { return propCfg; }

  JulianDate getEpoch() const noexcept { return jd0; }

  std::array<double, 6> getInitialState() const noexcept { return x0; }

  CoordType getCoordinateType() const noexcept { return coord; }

  FrameType getReferenceFrameType() const noexcept { return frame; }

private:
  std::string name;
  PropagatorConfig propCfg;
  CoordType coord;
  FrameType frame;
  JulianDate jd0;
  std::array<double, 6> x0;
};


}

#endif
