/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_RANGE_PRINTER_H
#define EOM_RANGE_PRINTER_H

#include <memory>
#include <string>
#include <deque>
#include <array>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>

#include <eom_command.h>
#include <eom_config.h>

namespace eom_app {

/**
 * EOM Command type that prints ephemeris of an orbit to a file.
 * Cartesian ITRF or GCRF are supported as outputs.
 *
 * @author  Kurt Motekew
 * @date    20211215
 */
class EomRangePrinter : public EomCommand {
public:
  /**
   * Converts string tokens into a command computing the range between
   * two ephemeris sources as a function of time.
   *
   * @param  tokens        Tokenized parameters with the orbit name, output
   *                       reference frame type (ITRF or GCRF), and output
   *                       filename.  Tokens are consumed as they are
   *                       used.
   * @param  jdEphStart    Time of first range output
   * @param  jdEphStop     Time of final range output
   * @param  ephemerides   List of ephemeris sources
   *
   * @throw  invalid_argument if exactly 3 tokens are not present, or an
   *         error is encountered (such as the input orbit name is not a
   *         valid orbit, or an invalid output reference frame is not
   *         selected.
   */
  EomRangePrinter(std::deque<std::string>& tokens, const EomConfig& cfg,
      const std::shared_ptr<std::unordered_map<std::string,
                            std::shared_ptr<eom::Ephemeris>>>& ephemerides);

  /**
   * Checks that listed ephemeris sources are valid.
   *
   * @throw  CmdValidateException if validation fails
   */
  void validate() override;

  /**
   * Writes .m function plotting the range between two ephemeris sources
   */
  void execute() const override;

private:
  std::array<std::string, 2> orbit_names;
  eom::EphemFrame frame;
  std::string func_name;
  std::string file_name;
  eom::JulianDate jdStart;
  eom::JulianDate jdStop;
  eom::Duration dtOut;
  std::string timeUnitsLbl;
  std::string distanceUnitsLbl;
  double to_time_units;
  double to_distance_units;
  std::array<std::shared_ptr<eom::Ephemeris>, 2> eph;
  std::shared_ptr<std::unordered_map<std::string,
                  std::shared_ptr<eom::Ephemeris>>> m_ephemerides;
};


}

#endif

