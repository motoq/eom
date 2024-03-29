/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_rtc_printer.h>

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <deque>
#include <unordered_map>
#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_attitude.h>

#include <eom_command.h>
#include <eom_config.h>

namespace eom_app {


EomRtcPrinter::EomRtcPrinter(std::deque<std::string>& tokens,
                             const EomConfig& cfg)
{
    // Read orbit name, etc.
  using namespace std::string_literals;
  if (tokens.size() != 3) {
    throw std::invalid_argument("EomRtcPrinter::EomRtcPrinter() "s +
                                "PrintRtc requires 3 arguments "s +
                                "vs. input "s +
                                std::to_string(tokens.size()));
  }
  for (int ii=0; ii<2; ++ii) {
    m_orbit_names[ii] = tokens[0];
    tokens.pop_front();
    if (!cfg.pendingOrbit(m_orbit_names[ii])) {
      throw std::invalid_argument("EomRtcPrinter::EomRtcPrinter() "s +
                                  "Invalid orbit name " + m_orbit_names[ii]);
    }
  }
  m_func_name = tokens[0];
  tokens.pop_front();
  m_file_name = m_func_name + ".m"s;
  m_jdStart = cfg.getStartTime();
  m_jdStop = cfg.getStopTime();
  m_dtOut = cfg.getOutputRate();
  m_distanceUnitsLbl = cfg.getIoDistansUnits();
  m_to_distance_units = cfg.getIoPerDu();
  m_to_time_units = cfg.getIoPerTu();
}


/*
 * Set ephemeris pointers using orbit_names from initialization
 */
void EomRtcPrinter::validate(const std::unordered_map<
    std::string, std::shared_ptr<eom::Ephemeris>>& ephemerides)
{
  using namespace std::string_literals;
  for (int ii=0; ii<2; ++ii) {
    try {
      m_eph[ii] = ephemerides.at(m_orbit_names[ii]);
    } catch (const std::out_of_range& oor) {
      throw CmdValidateException("EomRtcPrinter::validate() "s +
                                 "Invalid orbit name in PrintRange: "s +
                                 m_orbit_names[ii]);
    }
  }
}

void EomRtcPrinter::execute() const
{
  std::ofstream fout(m_file_name);

  if (fout.is_open()) {
    double tot_time {m_to_time_units*phy_const::tu_per_day*(m_jdStop -
                                                            m_jdStart)};
    double dt {m_to_time_units*m_dtOut.getTu()};
    unsigned long int nrec {static_cast<unsigned long int>(tot_time/dt)};
    nrec++;

      // Function header
    fout << "function [gxh, time_rtc] = " << m_func_name;
    fout << "\n% RTC is an EOM generated Matlab/Octave function that";
    fout << "\n% plots relative position as a function of time";
    fout << "\n%";
    fout << "\n% Outputs:";
    fout << "\n%   gxh       Graphics handle to new image";
    fout << "\n%   time_rtc  Nx4 matrix of time and RTC coordinates";
    fout << '\n';

      // Create time and range data
    fout << "\ntime_rtc = [";
    fout << std::scientific;
    fout.precision(16);
    for (unsigned long int ii=0UL; ii<nrec; ++ii) {
      if (ii > 0UL) {
        fout << ';';
      }
      double dtnow {ii*dt};
      fout << "\n  " << dtnow << " ";
      eom::JulianDate jdNow {m_jdStart +
                             phy_const::day_per_tu*(dtnow/m_to_time_units)};
      Eigen::Matrix<double, 6, 1> pv1 =
          m_eph[0]->getStateVector(jdNow, eom::EphemFrame::eci);
      Eigen::Matrix<double, 3, 1> r1 {pv1.block<3,1>(0,0)};
      Eigen::Matrix<double, 3, 1> r2 =
          m_eph[1]->getPosition(jdNow, eom::EphemFrame::eci);
      Eigen::Matrix<double, 3, 1> dr = r1 - r2;
      Eigen::Matrix<double, 3, 1> v1 {pv1.block<3,1>(3,0)};
      Eigen::Matrix<double, 3, 3> i2rtcDcm {eom::AttitudeRtc<double>(r1, v1)};
      dr = i2rtcDcm*dr;
      for (int jj=0; jj<3; ++jj) {
        fout << " " << m_to_distance_units*dr(jj);
      }
    }
      // Make the plot and annotate
    fout << "\n];";
    fout << "\nn = size(time_rtc,1);";
    fout << "\ngxh = figure; hold on;";
    fout << "\nplot3(time_rtc(:,2), time_rtc(:,3), time_rtc(:,4));";
    fout << "\nscatter3(time_rtc(1,2), time_rtc(1,3), time_rtc(1,4), 'g');";
    fout << "\nscatter3(time_rtc(n,2), time_rtc(n,3), time_rtc(n,4), 'r');";
    fout << "\nscatter3(0, 0, 0, 'b');";
    fout << "\nxlabel('Radial (" << m_distanceUnitsLbl << ")');";
    fout << "\nylabel('Transverse (" << m_distanceUnitsLbl << ")');";
    fout << "\nzlabel('Cross-Track (" << m_distanceUnitsLbl << ")');";
    fout << "\ntitle('" << m_orbit_names[0] << '-' << m_orbit_names[1] <<
            " RTC on " << m_jdStart.to_dmy_str() << "');";
    fout << "\naxis equal;";
    fout << "\nend\n";
    fout.close();
  } else {
    std::cerr << "\nCan't open " << m_file_name << '\n';
  }


}

}
