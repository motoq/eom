/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_EPHEM_PRINTER_H
#define EOM_EPHEM_PRINTER_H

#include <memory>
#include <string>
#include <deque>
#include <unordered_map>

#include <cal_julian_date.h>
#include <astro_orbit_def.h>
#include <astro_ephemeris.h>

#include <eom_command.h>
#include <eom_config.h>

namespace eom_app {

/**
 * EOM Command type that prints ephemeris of an orbit to a file.
 * Cartesian ITRF or GCRF are supported as outputs.
 *
 * @author  Kurt Motekew
 * @date    202111
 */
class EomEphemPrinter : public EomCommand {
public:
  ~EomEphemPrinter() = default;
  EomEphemPrinter(const EomEphemPrinter&) = default;
  EomEphemPrinter& operator=(const EomEphemPrinter&) = default;
  EomEphemPrinter(EomEphemPrinter&&) = default;
  EomEphemPrinter& operator=(EomEphemPrinter&&) = default;

  /**
   * Converts string tokens into an ephemeris print command.
   *
   * @param  tokens      Tokenized parameters with the orbit name, output
   *                     reference frame type (ITRF or GCRF), and output
   *                     filename.  Tokens are consumed as they are
   *                     used.
   * @param  cfg         Scenario configuration
   * @param  orbit_defs  Current orbit definitions for token validation
   *
   * @throws  invalid_argument if exactly 3 tokens are not present or 
   *          the indicated reference frame is not valid.  Orbit names
   *          will be checked during the validate step.
   */
  EomEphemPrinter(std::deque<std::string>& tokens,
                  const EomConfig& cfg,
                  const std::vector<eom::OrbitDef>& orbit_defs);


  /**
   * Checks that input ephemeris source is valid
   *
   * @param  ephemerides  All  available ephemeris resources
   *
   * @throws  CmdValidateException if validation fails (desired orbit
   *          name is not valid).
   */
  void validate(const std::unordered_map<
      std::string, std::shared_ptr<eom::Ephemeris>>& ephemerides) override;

  /**
   * Writes .e format ephemeris to the previously specified file.
   */
  void execute() const override;

private:
  eom::EphemFrame m_frame;
  std::string m_orbit_name;
  std::string m_file_name;
  eom::JulianDate m_jdStart;
  eom::JulianDate m_jdStop;
  std::shared_ptr<eom::Ephemeris> m_eph;
};


}

#endif

