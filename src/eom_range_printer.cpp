/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <eom_range_printer.h>

#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <vector>
#include <deque>
#include <stdexcept>

#include <Eigen/Dense>

#include <cal_const.h>
#include <phy_const.h>
#include <utl_units.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_print.h>

#include <eom_parse.h>

namespace eom_app {


EomRangePrinter::EomRangePrinter(std::deque<std::string>& tokens,
       const eom::JulianDate& jdEphStart, const eom::JulianDate& jdEphStop,
       const std::shared_ptr<std::unordered_map<std::string, int>>& ephem_ndxs,
       const std::shared_ptr<std::vector<std::shared_ptr<eom::Ephemeris>>>&
                                                                  ephem_list)
{
    // Read orbit name, output frame, and output filename
  using namespace std::string_literals;
  if (tokens.size() != 6) {
    throw std::invalid_argument("EomRangePrinter::EomRangePrinter:"s +
                                " PrintRange requires 6 arguments"s +
                                " vs. input "s +
                                std::to_string(tokens.size()));
  }
  for (int ii=0; ii<2; ++ii) {
    orbit_names[ii] = tokens[0];
    tokens.pop_front();
    try {
      endxs[ii] = ephem_ndxs->at(orbit_names[ii]);
    } catch (std::out_of_range& oor) {
      throw std::invalid_argument("EomRangePrinter::EomRangePrinter:"s +
                                  " Invalid orbit name in PrintRange: "s +
                                    orbit_names[ii]);
    }
  }
  dtout = parse_duration(tokens);
  file_name = tokens[0] + ".m"s;
  tokens.pop_front();
  units = tokens[0];
  tokens.pop_front();
  try {
    to_units = utl_units::per_du.at(units);
  } catch (std::out_of_range& oor) {
    throw std::invalid_argument("EomRangePrinter::EomRangePrinter:"s +
                                " Invalid units type: "s + units);
  }
  jdStart = jdEphStart;
  jdStop = jdEphStop;
  ephemerides = ephem_list;
}

void EomRangePrinter::execute() const
{
  std::ofstream fout(file_name.c_str());

  if (fout.is_open()) {
    double seconds {cal_const::sec_per_day*(jdStop - jdStart)};
    double dtsec   {phy_const::sec_per_tu*dtout.getTu()};
    unsigned long int nrec {static_cast<unsigned long int>(seconds/dtsec)};
    nrec++;

      // Function header
    fout << "function [gxh, time_range] = rng";
    fout << "\n% RNG is an EOM generated Matlab/Octave function that";
    fout << "\n% plots range as a function of time";
    fout << "\n%";
    fout << "\n% Outputs:";
    fout << "\n%   gxh         Graphics handle to new image";
    fout << "\n%   time_range  Nx2 matrix of time and range values";
    fout << '\n';

      // Create time and range data
    fout << "\ntime_range = [";
    fout << std::scientific;
    fout.precision(16);
    for (unsigned long int ii=0UL; ii<nrec; ++ii) {
      if (ii > 0UL) {
        fout << ';';
      }
      double tsec {ii*dtsec};
      fout << "\n  " << tsec << " ";
      eom::JulianDate jdNow {jdStart + cal_const::day_per_sec*tsec};
      Eigen::Matrix<double, 6, 1> pv1 =
        (*ephemerides)[endxs[0]]->getStateVector(jdNow, eom::EphemFrame::eci);
      Eigen::Matrix<double, 6, 1> pv2 =
        (*ephemerides)[endxs[1]]->getStateVector(jdNow, eom::EphemFrame::eci);
      Eigen::Matrix<double, 3, 1> r1 = pv1.block<3,1>(0,0);
      Eigen::Matrix<double, 3, 1> r2 = pv2.block<3,1>(0,0);
      Eigen::Matrix<double, 3, 1> dr = r1 - r2;
      double range = dr.norm();
      fout << " " << to_units*range;
    }
      // Make the plot and annotate
    fout << "\n];";
    fout << "\nfigure; hold on;";
    fout << "\nplot(time_range(:,1), time_range(:,2));";
    fout << "\nxlabel('seconds');";
    fout << "\nylabel('" << units << "');";
    fout << "\ntitle('" << orbit_names[0] << '-' << orbit_names[1] <<
            " Range on " << jdStart.to_dmy_str() << "');";
    fout << "\nend\n";
    fout.close();
  } else {
    std::cerr << "\nCan't open " << file_name << '\n';
  }


}

}
