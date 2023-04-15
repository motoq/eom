/*
 * Copyright 2023 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_orbit_printer.h>

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <deque>
#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>

#include <eom_command.h>
#include <eom_config.h>

namespace eom_app {


EomOrbitPrinter::EomOrbitPrinter(
    std::deque<std::string>& tokens, const EomConfig& cfg,
    const std::shared_ptr<std::unordered_map<std::string,
                          std::shared_ptr<eom::Ephemeris>>>& ephemerides)
{
  m_ephemerides = ephemerides;

    // Read orbit name, output frame, and output filename
  using namespace std::string_literals;
  if (tokens.size() != 2) {
    throw std::invalid_argument("EomOrbitPrinter::EomOrbitPrinter() "s +
                                "PrintOrbit requires 2 arguments "s +
                                "vs. input "s +
                                std::to_string(tokens.size()));
  }
  orbit_name = tokens[0];
  tokens.pop_front();
  func_name = tokens[0];
  tokens.pop_front();
  file_name = func_name + ".m"s;
  jdStart = cfg.getStartTime();
  jdStop = cfg.getStopTime();
  dtOut = cfg.getOutputRate();
}


/*
 * Set ephemeris pointer using orbit_name from initialization
 */
void EomOrbitPrinter::validate()
{
  using namespace std::string_literals;
  try {
    eph = m_ephemerides->at(orbit_name);
  } catch (const std::out_of_range& oor) {
    throw CmdValidateException("EomOrbitPrinter::validate() "s +
                               "Invalid orbit name in PrintRange: "s +
                               orbit_name);
  }
}

void EomOrbitPrinter::execute() const
{
  std::ofstream fout(file_name);

  if (fout.is_open()) {
    double tot_time {phy_const::tu_per_day*(jdStop - jdStart)};
    double dt {dtOut.getTu()};
    unsigned long int nrec {static_cast<unsigned long int>(tot_time/dt)};
    nrec++;

      // Function header
    fout << "function [gxh, time_xyz] = " << func_name;
    fout << "\n% Orbit is an EOM generated Matlab/Octave function that";
    fout << "\n% plots a 3D orbit trace";
    fout << "\n%";
    fout << "\n% Outputs:";
    fout << "\n%   gxh       Graphics handle to new image";
    fout << "\n%   time_xyz  Nx4 matrix of time and ECI coordinates, DU and TU";
    fout << '\n';

      // Create time and range data
    fout << "\ntime_xyz = [";
    fout << std::scientific;
    fout.precision(16);
    for (unsigned long int ii=0UL; ii<nrec; ++ii) {
      if (ii > 0UL) {
        fout << ';';
      }
      double dtnow {ii*dt};
      fout << "\n  " << dtnow << " ";
      eom::JulianDate jdNow {jdStart + phy_const::day_per_tu*dtnow};
      Eigen::Matrix<double, 6, 1> pv =
          eph->getStateVector(jdNow, eom::EphemFrame::eci);
      for (int jj=0; jj<6; ++jj) {
        fout << " " << pv(jj);
      }
    }
      // Make the plot and annotate
    fout << "\n];";
    fout << "\nn = size(time_xyz,1);";
    fout << "\ngxh = figure; hold on;";
    fout << "\nplot3(time_xyz(:,2), time_xyz(:,3), time_xyz(:,4));";
    fout << "\nscatter3(time_xyz(1,2), time_xyz(1,3), time_xyz(1,4), 'g');";
    fout << "\nscatter3(time_xyz(n,2), time_xyz(n,3), time_xyz(n,4), 'r');";
    fout << "\nscatter3(0, 0, 0, 'b');";
    fout << "\nplot3(time_xyz(:,5), time_xyz(:,6), time_xyz(:,7));";
    fout << "\nxlabel('X');";
    fout << "\nylabel('Y');";
    fout << "\nzlabel('Z');";
    fout << "\ntitle('" << orbit_name <<
            " on " << jdStart.to_dmy_str() << "');";
    fout << "\naxis equal;";
    fout << "\nend\n";
    fout.close();
  } else {
    std::cerr << "\nCan't open " << file_name << '\n';
  }


}

}
