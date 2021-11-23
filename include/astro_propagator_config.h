#ifndef ASTRO_PROPAGATOR_CONFIG_H
#define ASTRO_PROPAGATOR_CONFIG_H

namespace eom {


enum class PropagatorType {
  Kepler1
};


class PropagatorConfig {
public:
  PropagatorConfig() : prop {PropagatorType::Kepler1}
  {
  }

  PropagatorConfig(PropagatorType prop_type) : prop {prop_type}
  {
  }

  PropagatorType getPropagatorType() const noexcept { return prop; }

private:
  PropagatorType prop;
};


}

#endif
