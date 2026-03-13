/*
 * Copyright 2026 Kurt Motekew
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
 * EOM Command type that writes ECF to ECI transformation data to a
 * file for a requested date and time.
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
   * @param  tokens      Tokenized parameters consisting of a datetime
   *                     entry followed by the filename to for which
   *                     to which the information will be saved.
   * @param  cfg         Not used.
   *
   * @throws  invalid_argument if exactly n tokens are not present.
   */
  EomF2iMsgCmd(std::deque<std::string>& tokens, const EomConfig& cfg);

  /**
   * Not used for this command
   */
  void validate(const std::unordered_map<
      std::string, std::shared_ptr<eom::Ephemeris>>& ephemerides) override;

  /**
   * Moves a copy of the ECFECI service here
   *
   * @param  orbits       Not used
   * @param  ephemerides  Not used
   * @param  ecfeciSys    Used to generate ECFECI data
   *
   * @throws  CmdValidateException if the requested datetime of the
   *          ECFECI data is not covered by the EcfEciSys service.
   */
  void validate(const std::unordered_map<
      std::string, std::shared_ptr<eom::Ephemeris>>& ephemerides,
      const std::vector<eom::OrbitDef>& orbits,
      std::shared_ptr<const eom::EcfEciSys> ecfeciSys) override;

  /**
   * Writes file with requested ECF to ECI data
   */
  void execute() const override;

private:
    // During construction
  eom::JulianDate m_time;
  std::string m_filename;
    // Set during validate()
  std::shared_ptr<const eom::EcfEciSys> m_f2i {nullptr};
};


}

#endif

