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
  if (tokens.size() != 3) {
    throw std::invalid_argument("EomOrbitPrinter::EomOrbitPrinter() "s +
                                "PrintOrbit requires 3 arguments "s +
                                "vs. input "s +
                                std::to_string(tokens.size()));
  }
  orbit_name = tokens[0];
  tokens.pop_front();
  if (tokens[0] == "ECI") {
    frame = eom::EphemFrame::eci;
    tokens.pop_front();
  } else if (tokens[0] == "ECF") {
    frame = eom::EphemFrame::ecf;
    tokens.pop_front();
  } else {
    throw std::invalid_argument("EomOrbitPrinter::EomOrbitPrinter() "s +
                                "Invalid reference frame "s + tokens[0]);
  }
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
    fout << "function [gxh, tpv] = " << func_name;
    fout << "\n% Orbit is an EOM generated Matlab/Octave function that";
    fout << "\n% plots a 3D orbit trace in ECI or ECF coordinates";
    fout << "\n%";
    fout << "\n% Outputs:";
    fout << "\n%   gxh  Graphics handle to new image";
    fout << "\n%   tpv  Nx7 matrix of time, pos, vel, in ECI or ECF";
    fout << "\n%        coordinates, units of DU and DU/TU";
    fout << '\n';

      // Create time and range data
    fout << "\ntpv = [";
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
          eph->getStateVector(jdNow, frame);
      for (int jj=0; jj<6; ++jj) {
        fout << " " << pv(jj);
      }
    }

      // Make the plot and annotate
    std::string coords = "ECI";
    if (frame == eom::EphemFrame::ecf) {
      coords = "ECF";
    }
    fout << "\n];";
    fout << "\nn = size(tpv,1);";
    fout << "\ngxh = figure; hold on;";
    fout << "\nplot3(tpv(:,2), tpv(:,3), tpv(:,4), '.b');";
    fout << "\nscatter3(tpv(1,2), tpv(1,3), tpv(1,4), 'g');";
    fout << "\nscatter3(tpv(n,2), tpv(n,3), tpv(n,4), 'r');";
    fout << "\nscatter3(0, 0, 0, 'b');";
    fout << "\nplot3(tpv(:,5), tpv(:,6), tpv(:,7), '.m');";
    fout << "\nscatter3(tpv(1,5), tpv(1,6), tpv(1,7), 'g');";
    fout << "\nscatter3(tpv(n,5), tpv(n,6), tpv(n,7), 'r');";
    fout << "\nxlabel('X, dX/dT');";
    fout << "\nylabel('Y, dY/dT');";
    fout << "\nzlabel('Z, dZ/dT');";
    fout << "\ntitle('" << coords << " " << orbit_name <<
            " on " << jdStart.to_dmy_str() << "');";
    fout << "\naxis equal;";
    fout << "\nend\n";
    fout.close();
  } else {
    std::cerr << "\nCan't open " << file_name << '\n';
  }


}

}
