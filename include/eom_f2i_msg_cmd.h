/*
 * Copyright 2025 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_F2I_MSG_CMD_H
#define EOM_F2I_MSG_CMD_H

#include <deque>
#include <memory>
#include <string>
#include <unordered_map>

#include <cal_julian_date.h>
#include <astro_ecfeci_sys.h>
#include <astro_ephemeris.h>
#include <astro_orbit_def.h>

#include <eom_command.h>
#include <eom_config.h>

namespace eom_app {

/**
 * EOM Command type that writes ECF to ECI data to a requested date and
 * time.
 *
 * @author  Kurt Motekew
 * @date    2026/03/10
 */
class EomF2iMsgCmd : public EomCommand {
public:
  /**
   * Converts string tokens into a command indicating at what time
   * ECFECI data should be generated.
   *
   * @param  tokens      Tokenized parameters with 
   * @param  cfg         Scenario configuration
   *
   * @throws  invalid_argument if exactly n tokens are not present.
   */
  EomF2iMsgCmd(std::deque<std::string>& tokens, const EomConfig& cfg);

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
   * Writes .lis file with requested ECF to ECI data
   */
  void execute() const override;

private:
    // Initialization
  eom::JulianDate m_time;
  std::string m_filename;
    // Set during validate
  std::shared_ptr<const eom::EcfEciSys> m_f2i {nullptr};
};


}

#endif

