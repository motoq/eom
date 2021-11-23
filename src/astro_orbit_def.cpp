#include <astro_orbit_def.h>

#include <array>
#include <string>

#include <cal_julian_date.h>
#include <astro_propagator_config.h>

namespace eom {

OrbitDef::OrbitDef(std::string orbit_name, const PropagatorConfig& propConfig,
                   const JulianDate& epoch, const std::array<double, 6>& state,
                                    CoordType coord_type, FrameType frame_type)
{
  name = orbit_name;
  propCfg = propConfig,
  jd0 = epoch;
  x0 = state;
  coord = coord_type;
  frame = frame_type;
}


}
