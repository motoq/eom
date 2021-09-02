/*
 * Copyright 2016 Kurt Motekew
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef CAL_JULIAN_DATE_H
#define CAL_JULIAN_DATE_H

#include <string>

#include <cal_const.h>
#include <cal_greg_date.h>

namespace eom {

/**
 * A Julian Date class designed to preserve precision by splitting the date
 * into high and low values along with an optional seconds storage counter.
 *
 * @author  Kurt Motekew
 * @date    20160314
 */
class JulianDate {
public:
  /** Initialize with J2000 */
  JulianDate() {}

  /**
   * See equivalent set method for details.
   */
  JulianDate(double jdDays, double jdFrac = 0.0, double jdSeconds = 0.0);

  /**
   * See equivalent set method for details.
   */
  JulianDate(const JulianDate& jd);

  /**
   * See equivalent set method for details.
   */
  JulianDate(const GregDate& gd, int hr = 0, int min = 0, double sec = 0.0);

  /**
   * Initialize using Julian Date components
   *
   * @param   jdDays      Days portion of Julian Date, or full JD if other
   *                      parameters are zero
   * @param   jdFrac      Fractional days to add to jdDays
   * @param   jsSeconds   Seconds to add to jdDays and jdFrac
   */
  void set(double jdDays, double jdFrac = 0.0, double jdSeconds = 0.0);

  /**
   * Initialize using an existing JulianDate.
   *
   * @param   jd   JD to copy
   */
  void set(const JulianDate& jd);

  /**
   * Initialize with a Gregorian date and time of day.  The time parameters
   * are not bounded typical clock 12/24 and 60 base systems - as many
   * hours, minutes, and seconds from the gd, +/-, are allowed and
   * incorporated into the internal JD representation.
   *
   * @param   gd    Gregorian Date
   * @param   hr    Hours from gd
   * @param   min   Minutes from gd and hr
   * @param   sec   Seconds from gd, hr, and min
   */
  void set(const GregDate& gd, int hr = 0, int min = 0, double sec = 0.0);

  /**
   * @return   Julian Date, scalar form
   */
  double getJd() const;

  /**
   * @return   Modified Julian Date, scalar
   */
  double getMjd() const;

  /**
   * @return   Large portion of the Julian Date, on the order of
   *           2,400,000.  Ideally, this would be the Julian Date
   *           corresponding to noon, but this is not necessary,
   *           depending on the value of jdLowVal.  Units of days.
   */
  double getJdHigh()   const { return jdHi; }

  /**
   * @return   Small portion of the Julian Date, typically on the order
   *           of a fraction of a day.  Units of days.
   */
  double getJdLow() const;

  /**
   * @param   days   Days to add to this Julian Date (or subtracted, if
   *                 negative).
   */
  JulianDate& operator+=(double days);

  /**
   * @param   days   Days to add (or subtract, if negative).
   *
   * @return   Copy of this JulianDate, adjusted by days
   */
  JulianDate operator+(double days);

  /**
   * @return  The days and fraction difference between the two Julian
   *          dates:  LargerJD - SmallerJD = PositiveDays
   */
  double operator-(const JulianDate& jd);

    /**
     * @return   Gregorian Date and time as a string.  Time is
     *           in base 24:60:60 format.
     */
  std::string to_str();

private:
  double jdHi {cal_const::J2000};         // Days
  double jdLow {0.0};                     // Fraction of a day
  double jdSec {0.0};                     // Seconds

  double gd2jd(int year, int month, int day);
  void normalize();
  void jd2gd(int& year, int& month, int& day,
             int& hour, int& minutes, double& seconds);
};

}

#endif
