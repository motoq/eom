/*
 * Copyright 2021 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <astro_print.h>

#include <memory>
#include <string>
#include <iostream>
#include <fstream>

#include <Eigen/Dense>

#include <cal_const.h>
#include <cal_duration.h>
#include <cal_julian_date.h>
#include <phy_const.h>
#include <astro_ephemeris.h>

namespace eom {

void print_ephemeris(std::string file_name,
                     const JulianDate& jdStart, const JulianDate& jdStop,
                     const Duration& dtout, EphemFrame frame,
                     const std::shared_ptr<const Ephemeris>& orbit)
{
  std::ofstream fout(file_name);
  if (fout.is_open()) {
    double seconds {cal_const::sec_per_day*(jdStop - jdStart)};
    double dtsec   {phy_const::sec_per_tu*dtout.getTu()};

    unsigned long int nrec {static_cast<unsigned long int>(seconds/dtsec)};
    nrec++;
    nrec++;

    fout << "stk.v.11.0";
    fout << "\n\n# WrittenBy    eom";
    fout << "\n\nBEGIN Ephemeris";
    fout << "\n\n    NumberOfEphemerisPoints\t\t " << nrec;
    fout << "\n\n    ScenarioEpoch\t\t " << jdStart.to_dmy_str();
    fout << "\n\n\n    InterpolationMethod\t\t Lagrange";
    fout << "\n\n    InterpolationSamplesM1\t\t 7";
    fout << "\n\n    CentralBody\t\t Earth";
    if (frame == EphemFrame::eci) {
      fout << "\n\n    CoordinateSystem\t\t ICRF";
    } else {
      fout << "\n\n    CoordinateSystem\t\t Fixed";
    }
    fout << "\n\n    EphemerisTimePosVel";


    fout << std::scientific;
    fout.precision(16);

    fout << '\n';
    for (unsigned long int ii=0UL; ii<nrec; ++ii) {
      double tsec {ii*dtsec};
      fout << "\n " << tsec << " ";
      JulianDate jdNow {jdStart + cal_const::day_per_sec*tsec};
      Eigen::Matrix<double, 6, 1> pv = orbit->getStateVector(jdNow, frame);
      for (int jj=0; jj<3; ++jj) {
        fout << phy_const::m_per_du*pv(jj) << "  ";
      }
      for (int jj=3; jj<6; ++jj) {
        fout << phy_const::m_per_du*pv(jj)*phy_const::tu_per_sec << "  ";
      }
    }
    fout << "\n\n\nEND Ephemeris\n";
    fout.close();
  } else {
    std::cerr << "\nCan't open " << file_name << '\n';
  }

}


}

