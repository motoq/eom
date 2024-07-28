/*
 * Copyright 2024 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#include <iostream>
#include <string>

/*
#include <fstream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
*/

#include <utl_const.h>
#include <cal_greg_date.h>
#include <cal_julian_date.h>
#include <astro_sgp4.h>
#include <astro_tle.h>

#include <SGP4.h>

/**
 * @author  Kurt Motekew
 * @date    2024/07/22  Initial
 */
int main(int argc, char* argv[])
{
  if (argc != 2) {
    std::cerr << "\nProper use is:  " << argv[0] << " <tle_filename>\n";
    return 0;
  }

    // External input...
  std::string tle1 =
      "1 24208U 96044A   06177.04061740 -.00000094 -00004-6  10002-3 0  1600";
  std::string tle2 =
      "2 24208   3.8536  80.0121 0026640 311.0977  48.3000  1.00778054 36119";
  
    // Ingest TLE
  eom::Tle tle(tle1, tle2);
  std::cout << '\n' << tle << '\n';

    // Config
  gravconsttype gctype {gravconsttype::wgs72};
  char opsmode {'a'};

  elsetrec satrec;

  SGP4Funcs::sgp4init(gctype,
                      opsmode,
                      tle.getSatN().data(),
                      tle.getTleSgpEpoch(),
                      0.0*tle.getBstar(),
                      0.0*tle.getMeanMotionRate(),
                      0.0*tle.getMeanMotionSecondRate(),
                      tle.getEccentricity(),
                      utl_const::deg_per_rad*tle.getArgumentOfPerigee(),
                      utl_const::deg_per_rad*tle.getInclination(),
                      utl_const::deg_per_rad*tle.getMeanAnomaly(),
                      (utl_const::tpi/1440.0)*tle.getMeanMotion(),
                      utl_const::deg_per_rad*tle.getRaan(),
                      satrec);

  double r_teme[3];
  double v_teme[3];
  double tsince {0.0};                 // Minutes since TLE epoch
  if (SGP4Funcs::sgp4(satrec, tsince, r_teme, v_teme)) {
    std::cout << "Position: " << r_teme[0] << '\n' <<
                 "          " << r_teme[1] << '\n' <<
                 "          " << r_teme[2] << " km\n";
    std::cout << "Velocity: " << v_teme[0] << '\n' <<
                 "          " << v_teme[1] << '\n' <<
                 "          " << v_teme[2] << " km/sec";
  } else {
    std::cout << "\n\nError with TLE Conversion\n";
  }

  std::cout << '\n';
}
