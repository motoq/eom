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

#include <cal_julian_date.h>
#include <astro_ephemeris.h>

#include <eom_command.h>
#include <eom_config.h>

namespace eom_app {

/**
 * EOM Command type that creates a Matlab/Octave function that plots the
 * position and velocity of a satellite, relative to ECI, in DU and
 * DU/TU.
 *
 * @author  Kurt Motekew
 * @date    2023/04/14
 */
class EomOrbitPrinter : public EomCommand {
public:
  ~EomOrbitPrinter() = default;
  EomOrbitPrinter(const EomOrbitPrinter&) = default;
  EomOrbitPrinter& operator=(const EomOrbitPrinter&) = default;
  EomOrbitPrinter(EomOrbitPrinter&&) = default;
  EomOrbitPrinter& operator=(EomOrbitPrinter&&) = default;

  /**
   * Converts string tokens into a command
   *
   * @param  tokens        Tokenized parameters with the orbit name and
   *                       output filename prefix.
   * @param  jdEphStart    Time of first 3D position and velocity output
   * @param  jdEphStop     Time of final 3D position and velocity output
   * @param  ephemerides   List of ephemeris sources
   *
   * @throws  invalid_argument if exactly 2 tokens are not present, or 
   *          the parsed orbit name is not a valid orbit.
   */
  EomOrbitPrinter(std::deque<std::string>& tokens, const EomConfig& cfg,
                  const std::shared_ptr<std::unordered_map<std::string,
                  std::shared_ptr<eom::Ephemeris>>>& ephemerides);

  /**
   * Checks that listed ephemeris sources are valid.
   *
   * @throws  CmdValidateException if validation fails
   */
  void validate() override;

  /**
   * Writes .m function
   */
  void execute() const override;

private:
  std::string orbit_name;                   // Name of orbit to print
  std::string func_name;                    // Function name (file prefix)
  std::string file_name;                    // func_name.m
  eom::JulianDate jdStart;
  eom::JulianDate jdStop;
  eom::Duration dtOut;
  std::shared_ptr<eom::Ephemeris> eph;
  std::shared_ptr<std::unordered_map<std::string,
                  std::shared_ptr<eom::Ephemeris>>> m_ephemerides;
};


}

#endif

