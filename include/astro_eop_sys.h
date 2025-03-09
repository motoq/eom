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

// Raw daily IERS EOP data - left in published units
struct eop_record {
  long mjd {0L};             ///< Modified Julian Date
  float mjdf {0};            ///< Fractional day - zero unless interpolated
  double xp {0.0};           ///< Polar motion, arcseconds (")
  double yp {0.0};
  double ut1mutc {0.0};      ///< UT1-UTC, seconds
  double lod {0.0};          ///< Length of day, milliseconds
  double dx {0.0};           ///<  Corrections to pole location,
  double dy {0.0};           ///<  miliiarcseconds (MAS)
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

  /**
   * @param  jd  Time for which to return EOP data, UTC
   *
   * @return  EOP for the requested time.  Daily published values are
   *          interpolated.  Continuity is guaranteed.
   */
  eop_record getEop(const JulianDate& jd) const;

  /**
   * @return iterator pointing to beginning of EOP data
   */
  auto cbegin() const
  {
    return eopData.cbegin(); 
  }

  /**
   * @return iterator pointing to just after the last EOP data
   *         structure
   */
  auto cend() const
  {
    return eopData.cend();
  }

  /**
   * @return  The number of eop records
   */
  auto size() const
  {  
    return eopData.size();
  }

private:
  unsigned long mjd_first {0UL};
  unsigned long mjd_last {0UL};
  std::vector<eop_record> eopData;
};


}

#endif
