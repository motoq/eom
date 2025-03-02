/*
 * Copyright 2021, 2023 Kurt Motekew
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
#include <unordered_map>
#include <array>

#include <cal_julian_date.h>
#include <astro_ephemeris.h>

#include <eom_command.h>
#include <eom_config.h>

namespace eom_app {

/**
 * EOM Command type that creates a Matlab/Octave function that plots the
 * range between two satellites over time.  Optionally, the amplitude
 * spectrum of the range is plotted for perturbation analysis.
 *
 * @author  Kurt Motekew
 * @date    2021/12/15
 * @date    2023/07/08  Added amplitude spectrum option
 */
class EomRangePrinter : public EomCommand {
public:
  /**
   * Converts string tokens into a command computing the range between
   * two ephemeris sources as a function of time.
   *
   * @param  tokens      Tokenized parameters with the orbit name, output
   *                     reference frame type (ITRF or GCRF), and output
   *                     filename.  Tokens are consumed as they are
   *                     used.
   * @param  cfg         Scenario configuration
   * @param  do_spectrum  Also create range amplitude spectrum plot.
   *                      This option is geared towards comparison
   *                      of the same orbit propagated via methods
   *                      capturing different perturbation effects.
   *
   * @throws  invalid_argument if exactly 3 tokens are not present.
   *          Orbit names will be checked during the validate step.
   */
  EomRangePrinter(std::deque<std::string>& tokens,
                  const EomConfig& cfg,
                  bool do_spectrum = false);

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

  /**
   * Writes .m function plotting the range between two ephemeris sources
   */
  void execute() const override;

private:
  std::array<std::string, 2> m_orbit_names;
  std::string m_func_name;
  std::string m_file_name;
  eom::JulianDate m_jdStart;
  eom::JulianDate m_jdStop;
  eom::Duration m_dtOut;
  std::string m_timeUnitsLbl;
  std::string m_distanceUnitsLbl;
  bool m_spectrum;
  double m_to_time_units;
  double m_to_distance_units;

  std::array<std::shared_ptr<eom::Ephemeris>, 2> m_eph;
};


}

#endif

