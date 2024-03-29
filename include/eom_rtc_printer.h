/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef EOM_RTC_PRINTER_H
#define EOM_RTC_PRINTER_H

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
 * position of one satellite relative to another in Cartesian radial,
 * transverse, cross-track coordinates (RTC, RSW, etc.).  The
 * orientation of the RTC reference frame is based on the inertial
 * velocity vector of the chief satellite.
 *
 * @author  Kurt Motekew
 * @date    2022/01/09
 */
class EomRtcPrinter : public EomCommand {
public:
  ~EomRtcPrinter() = default;
  EomRtcPrinter(const EomRtcPrinter&) = default;
  EomRtcPrinter& operator=(const EomRtcPrinter&) = default;
  EomRtcPrinter(EomRtcPrinter&&) = default;
  EomRtcPrinter& operator=(EomRtcPrinter&&) = default;

  /**
   * Converts string tokens into a command computing the position of a
   * deputy satellite w.r.t. a chief satellite using the RTC reference
   * frame.
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
  EomRtcPrinter(std::deque<std::string>& tokens, const EomConfig& cfg);

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
  std::string m_distanceUnitsLbl;
  double m_to_time_units;
  double m_to_distance_units;

  std::array<std::shared_ptr<eom::Ephemeris>, 2> m_eph;
};


}

#endif

