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

#include <cal_julian_date.h>
#include <astro_ephemeris.h>

#include <eom_command.h>

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
   * @param  tokens        Tokenized parameters with the orbit name, output
   *                       reference frame type (ITRF or GCRF), and output
   *                       filename.  Tokens are consumed as they are
   *                       used.
   * @param  jdEphStart    Time of first ephemeris output
   * @param  jdEphStop     Time of final ephemeris output
   * @param  ephemerides   Ordered list of ephemeris sources
   *
   * @throws  invalid_argument if exactly 3 tokens are not present or 
   *          the indicated reference frame is not valid.  Orbit names
   *          will be checked during the validate step.
   */
  EomEphemPrinter(std::deque<std::string>& tokens,
      const eom::JulianDate& jdEphStart, const eom::JulianDate& jdEphStop,
      const std::shared_ptr<std::unordered_map<std::string,
                            std::shared_ptr<eom::Ephemeris>>>& ephemerides);


  /**
   * Checks that input ephemeris source is valid
   *
   * @throws  CmdValidateException if validation ails (desired orbit
   *          name is not valid).
   */
  void validate() override;

  /**
   * Writes .e format ephemeris to the previously specified file.
   */
  void execute() const override;

private:
  int endx;                                 // Index into ephemerides
  eom::EphemFrame frame;
  std::string orbit_name;
  std::string file_name;
  eom::JulianDate jdStart;
  eom::JulianDate jdStop;
  std::shared_ptr<eom::Ephemeris> eph;
  std::shared_ptr<std::unordered_map<std::string,
                  std::shared_ptr<eom::Ephemeris>>> m_ephemerides;
};


}

#endif

