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
#include <vector>
#include <deque>
#include <array>

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
   * @param  ephem_ndxs    Locates orbit index given the orbit name
   * @param  ephem_list    Ordered list of ephemeris sources
   *
   * @throw  invalid_argument if exactly 3 tokens are not present, or an
   *         error is encountered (such as the input orbit name is not a
   *         valid orbit, or an invalid output reference frame is not
   *         selected.
   */
  EomRangePrinter(std::deque<std::string>& tokens,
      const eom::JulianDate& jdEphStart, const eom::JulianDate& jdEphStop,
      const std::shared_ptr<std::unordered_map<std::string, int>>& ephem_ndxs,
      const std::shared_ptr<std::vector<std::shared_ptr<eom::Ephemeris>>>&
                                                                    ephem_list);

  /**
   * Writes .e format ephemeris to the previously specified file.
   */
  void execute() const override;

private:
  std::array<double, 2> endxs;              // Index into ephemeris sources
  eom::EphemFrame frame;
  std::string file_name;
  eom::JulianDate jdStart;
  eom::JulianDate jdStop;
  eom::Duration dtout;
  std::string units;
  double to_units;
  std::shared_ptr<std::vector<std::shared_ptr<eom::Ephemeris>>> ephemerides;
};


}

#endif

