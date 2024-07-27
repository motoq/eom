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

  const eom::GregDate gd0 {1950, 1, 1};
  const eom::JulianDate jd0 {gd0};

  std::string tle1 =
      "1 24208U 96044A   06177.04061740 -.00000094 -00004-6  10002-3 0  1600";
  std::string tle2 =
      "2 24208   3.8536  80.0121 0026640 311.0977  48.3000  1.00778054 36119";
  
  eom::Tle tle(tle1, tle2);
  std::cout << '\n' << tle << '\n';

  eom::JulianDate jd;


  gravconsttype gctype {gravconsttype::wgs84};
  char opsmode {'a'};
  char satn[5];
  strcpy(satn, "666");
  double epoch {jd - jd0};
  double xbstar {0.0};
  double xndot {0.0};
  double xnddot {0.0};
  double xecco {0.0026640};
  double xargpo {utl_const::rad_per_deg*311.0977};
  double xinclo {utl_const::rad_per_deg*3.8536};
  double xmo {utl_const::rad_per_deg*48.3000};
  double xno_kozai {1.00778054};
  double xnodeo {utl_const::rad_per_deg*80.0121};

  elsetrec satrec;

  SGP4Funcs::sgp4init(gctype,
                      opsmode,
                      satn,
                      epoch,
                      xbstar,
                      xndot,
                      xnddot,
                      xecco,
                      xargpo,
                      xinclo,
                      xmo,
                      xno_kozai,
                      xnodeo,
                      satrec);

  std::cout << '\n';
}
