/*
 * Copyright 2022 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef ASTRO_EOP_SYS_H
#define ASTRO_EOP_SYS_H

#include <string>
#include <vector>

#include <cal_julian_date.h>

namespace eom {

// Raw daily IERS EOP data
struct eop_record {
  long mjd {0L};
  double xp {0.0};
  double yp {0.0};
  double ut1mutc {0.0};
  double lod {0.0};
  double dx {0.0};
  double dy {0.0};
};

/**
 * System resource utility for IERS EOP data.
 *
 * @author  Kurt Motekew
 * @date    20220702
 */
class EopSys {
public:
  /**
   * This constructor parses an IERS semicolon separated values (.csv)
   * formatted file of EOP data.  The first line of the file is expected
   * to include the column labels "MJD", "x_pole" (arcsec), "y_pole"
   * (arcsec), "UT1-UTC" (sec), * "LOD" (millisec), "dX" (MAS), and "dY"
   * (MAS).  The file can be found from the IERS EOP datacenter finals.all
   * (IAU2000) version metadata link to the finals2000A.all.csv file.
   *
   * The first column is expected to be an integer MJD value, with each
   * line of data separated by one day.
   *
   * @param  fname        Filename to open and close, from which to
   *                      parse IERS data.  The IERS "csv" semicolon
   *                      separated value format (uses ';' delimiters)
   * @param  startTime    Earliest UTC time for which to parse and store
   *                      IERS EOP data.
   * @param  stopTime     Latest UTC time for which to parse and store
   *                      IERS EOP data.
   */
  EopSys(std::string fname,
         const JulianDate& startTime, const JulianDate& stopTime);

private:
  unsigned long mjd_first {0UL};
  unsigned long mjd_last {0UL};
  std::vector<eop_record> eopData;
};


}

#endif