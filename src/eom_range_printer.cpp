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
#include <iomanip>
#include <memory>
#include <string>
#include <deque>
#include <stdexcept>

#include <Eigen/Dense>

#include <phy_const.h>
#include <cal_julian_date.h>
#include <astro_ephemeris.h>
#include <astro_keplerian.h>

#include <eom_command.h>
#include <eom_config.h>

namespace eom_app {

EomRangePrinter::
EomRangePrinter(std::deque<std::string>& tokens, 
                const EomConfig& cfg,
                const std::shared_ptr<std::unordered_map<std::string,
                std::shared_ptr<eom::Ephemeris>>>& ephemerides,
                bool do_spectrum)
{
  m_ephemerides = ephemerides;
  m_spectrum = do_spectrum;

    // Read orbit name, output frame, and output filename
  using namespace std::string_literals;
  if (tokens.size() != 3) {
    throw std::invalid_argument("EomRangePrinter::EomRangePrinter() "s +
                                "PrintRange requires 3 arguments "s +
                                "vs. input "s +
                                std::to_string(tokens.size()));
  }
  for (int ii=0; ii<2; ++ii) {
    orbit_names[ii] = tokens[0];
    tokens.pop_front();
  }
  func_name = tokens[0];
  tokens.pop_front();
  file_name = func_name + ".m"s;
  jdStart = cfg.getStartTime();
  jdStop = cfg.getStopTime();
  dtOut = cfg.getOutputRate();
  timeUnitsLbl = cfg.getIoTimeUnits();
  distanceUnitsLbl = cfg.getIoDistansUnits();
  to_distance_units = cfg.getIoPerDu();
  to_time_units = cfg.getIoPerTu();
}


/*
 * Set ephemeris pointers using orbit_names from initialization
 */
void EomRangePrinter::validate()
{
  using namespace std::string_literals;
  for (int ii=0; ii<2; ++ii) {
    try {
      eph[ii] = m_ephemerides->at(orbit_names[ii]);
    } catch (const std::out_of_range& oor) {
      throw CmdValidateException("EomRangePrinter::validate() "s +
                                 "Invalid orbit name in PrintRange: "s +
                                 orbit_names[ii]);
    }
  }
}

void EomRangePrinter::execute() const
{
  std::ofstream fout(file_name);

  if (fout.is_open()) {
    double tot_time {to_time_units*phy_const::tu_per_day*(jdStop - jdStart)};
    double dt {to_time_units*dtOut.getTu()};
    unsigned long int nrec {static_cast<unsigned long int>(tot_time/dt)};
    nrec++;

      // Function header
    fout << "function [gxh, time_range] = " << func_name;
    fout << "(detrend_opt)";
    fout << "\n% RNG is an EOM generated Matlab/Octave function that";
    fout << "\n% plots range as a function of time.";
    if (m_spectrum) {
      fout << " The amplitude";
      fout << "\n% spectrum is also plotted.";
    }
    fout << "\n%";
    fout << "\n% Input:";
    fout << "\n%   detrend_opt  If > 0, perform indicated detrend on fft";
    fout << "\n%                operation.  Optional input, no detrend by";
    fout << "\n%                default";
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
      double dtnow {ii*dt};
      fout << "\n  " << dtnow << " ";
      eom::JulianDate jdNow {jdStart +
                             phy_const::day_per_tu*(dtnow/to_time_units)};
      Eigen::Matrix<double, 3, 1> r1 =
              eph[0]->getPosition(jdNow, eom::EphemFrame::eci);
      Eigen::Matrix<double, 3, 1> r2 =
              eph[1]->getPosition(jdNow, eom::EphemFrame::eci);
      Eigen::Matrix<double, 3, 1> dr = r1 - r2;
      double range = dr.norm();
      fout << " " << to_distance_units*range;
    }
      // Make the plot and annotate
    fout << "\n];";
    fout << "\ngxh = figure; hold on;";
    fout << "\nplot(time_range(:,1), time_range(:,2));";
    fout << "\nxlabel('" << timeUnitsLbl << "');";
    fout << "\nylabel('" << distanceUnitsLbl << "');";
    fout << "\ntitle('" << orbit_names[0] << '-' << orbit_names[1] <<
            " Range on " << jdStart.to_dmy_str() << "');";
    fout << "\naxis tight";

      // Create and plot amplitude spectrum using Matlab/Octave fft
    if (m_spectrum) {
      eom::Keplerian kep(eph[0]->getStateVector(jdStart, eom::EphemFrame::eci));
      auto dt_revs = dtOut.getTu()/kep.getPeriod();
      fout << "\n\nrate = 1/" << dt_revs << ';';
      fout << "\nnsamp = size(time_range,1);";
      fout << "\nrng_freqs = rate*(0:(nsamp/2))/nsamp;";
      fout << std::fixed;
      fout << std::setprecision(2);
      fout << "\nif nargin > 0";
      fout << "\n  rngs = detrend(time_range(:,2)', detrend_opt);";
      fout << "\n  stitle = sprintf('Amplitude Spectrum, " <<
                   phy_const::tu_per_day/kep.getPeriod() << " rev/day,";
      fout << " detrend = %i', detrend_opt);";
      fout << "\nelse";
      fout << "\n  rngs = time_range(:,2)';";
      fout << "\n  stitle = 'Amplitude Spectrum, " <<
                   phy_const::tu_per_day/kep.getPeriod() << " rev/day';";
      fout << "\nend";
      fout << "\nfft_rngs = fft(rngs);";
      fout << "\nP2 = 2.0*abs(fft_rngs/nsamp);";
      fout << "\nP1 = P2(1:nsamp/2+1);";
        // Make the plot and annotate
      fout << "\nfigure; hold on;";
      fout << "\nplot(rng_freqs, P1);";
      fout << "\nxlabel('revs');";
      fout << "\nylabel('" << distanceUnitsLbl << "');";
      fout << "\ntitle(stitle);";
      fout << "\naxis tight";
      fout << "\ngrid on;";
      fout << "\ngrid minor on;";
    }

    fout << "\nend\n";

    fout.close();
  } else {
    std::cerr << "\nCan't open " << file_name << '\n';
  }


}

}
