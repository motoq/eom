/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_XFER_ORBIT_H
#define EOM_XFER_ORBIT_H

#include <deque>
#include <memory>
#include <string>
#include <unordered_map>

#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_orbit_def.h>
#include <astro_propagator_config.h>

#include <eom_command.h>
#include <eom_config.h>

namespace eom_app {

/**
 * EOM Command type that creates a Matlab/Octave function that plots the
 *
 * @author  Kurt Motekew
 * @date    2025/01/05
 */
class EomXferOrbit : public EomCommand {
public:
  /**
   * Converts string tokens into a command computing generating a
   * transfer orbit based on initial and destination orbits.
   *
   * @param  tokens      Tokenized parameters with the orbit name, output
   *                     reference frame type (ITRF or GCRF), and output
   *                     filename.  Tokens are consumed as they are
   *                     used.
   * @param  cfg         Scenario configuration
   *
   * @throws  invalid_argument if exactly 3 tokens are not present.
   *          Orbit names will be checked during the validate step.
   */
  EomXferOrbit(std::deque<std::string>& tokens, const EomConfig& cfg);

  /**
   * Checks that listed ephemeris sources are valid.
   *
   * @param  ephemerides  All  available ephemeris resources
   *
   * @throws  CmdValidateException if validation fails (invalid orbit
   *          name encountered).
   */
  void validate(const std::unordered_map<
      std::string, std::shared_ptr<eom::Ephemeris>>& ephemerides) override;

  void validate(const std::unordered_map<
      std::string, std::shared_ptr<eom::Ephemeris>>& ephemerides,
      const std::vector<eom::OrbitDef>& orbits,
      std::shared_ptr<const eom::EcfEciSys> ecfeciSys) override;

  /**
   * Writes .m function plotting the range between two ephemeris sources
   */
  void execute() const override;

private:
    // Initialization
  std::string m_start_orbit_name {""};
  std::string m_end_orbit_name {""};
  eom::JulianDate m_xfer_start;
  eom::Duration m_xfer_dur;
  std::string m_func_name;
  eom::Duration m_dtOut;
  std::string m_distanceUnitsLbl;
  double m_to_time_units;
  double m_to_distance_units;
    // Set during validate
  std::shared_ptr<const eom::EcfEciSys> m_f2i {nullptr};
  std::shared_ptr<eom::Ephemeris> m_start_eph;
  std::shared_ptr<eom::Ephemeris> m_end_eph;
  eom::PropagatorConfig m_propCfg;
};


}

#endif

