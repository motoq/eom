/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_ORBIT_PRINTER_H
#define EOM_ORBIT_PRINTER_H

#include <memory>
#include <string>
#include <deque>
#include <unordered_map>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>

#include <eom_command.h>
#include <eom_config.h>

namespace eom_app {

/**
 * EOM Command type that creates a Matlab/Octave function that plots the
 * 3D position and velocity of a satellite in DU and DU/TU.
 *
 * @author  Kurt Motekew
 * @date    2023/04/15
 */
class EomOrbitPrinter : public EomCommand {
public:
  /**
   * Converts string tokens into a command
   *
   * @param  tokens      Tokenized parameters with the orbit name,
   *                     reference frame to use, and the output filename
   *                     prefix.
   * @param  cfg         Scenario configuration
   *
   * @throws  invalid_argument if exactly 3 tokens are not present or 
   *          the indicated reference frame is not valid.  Orbit names
   *          will be checked during the validate step.
   */
  EomOrbitPrinter(std::deque<std::string>& tokens,
                  const EomConfig& cfg);

  /**
   * Checks that listed ephemeris sources are valid.
   *
   * @param  ephemerides  All  available ephemeris resources
   *
   * @throws  CmdValidateException if validation fails (desired orbit
   *          name is not valid).
   */
  void validate(const std::unordered_map<
      std::string, std::shared_ptr<eom::Ephemeris>>& ephemerides) override;

  /**
   * Writes .m function
   */
  void execute() const override;

private:
  std::string m_orbit_name;                 // Name of orbit to print
  std::string m_func_name;                  // Function name (file prefix)
  std::string m_file_name;                  // func_name.m
  eom::JulianDate m_jdStart;
  eom::JulianDate m_jdStop;
  eom::Duration m_dtOut;
  eom::EphemFrame m_frame;

  std::shared_ptr<eom::Ephemeris> m_eph;
};


}

#endif

